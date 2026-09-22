"""Independent verifier for one sealed trial (V0..V6 of the frozen grading rules).

    verify_run.py --run-id <id>

Works only on copies under control/verifier/<run_id>/; never writes into arena/ or control/runs/.
Deterministic checks only. Items that need a human crystallographer are reported as unknown.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
SHELXL = ROOT / "toolbox" / "shelx" / "shelxl.exe"
PLATON = ROOT / "toolbox" / "shelx" / "platon.exe"
GEMMI_PY = Path(r"H:\CrystalPilot\.venv\Scripts\python.exe")   # gemmi lives here; used as a parser only
R1_TOL = 1e-4


def sha256_of(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


# ------------------------------------------------------------------ V0 seal & locate the delivery
def locate_delivery(arena: Path, cond: str) -> dict:
    """The agent's declared delivery: CrystalPilot's results folder (H0/H1) or deliverables/ (C0/P1).
    Falls back to any structure file in the arena, flagged as 'collected, not formally delivered'."""
    cands = []
    exts = (".cif", ".res", ".ins", ".hkl", ".fcf", ".fab", ".lst")
    results = arena / "CrystalPilot Results"
    formal = True
    if cond in ("H0", "H1") and results.is_dir():
        finals = sorted(results.rglob("final.cif"), key=lambda p: p.stat().st_mtime)
        if finals:
            cands = [p for p in finals[-1].parent.iterdir() if p.is_file()]
    if not cands and (arena / "deliverables").is_dir():
        cands = [p for p in (arena / "deliverables").rglob("*") if p.is_file()]
    if not any(p.suffix.lower() in (".cif", ".res") for p in cands):
        formal = False
        cands = [p for p in arena.rglob("*") if p.is_file() and p.suffix.lower() in exts
                 and "inputs" not in p.relative_to(arena).parts[:1] and ".crystalpilot" not in p.parts]
    files = {}
    for p in cands:
        try:
            files[p.relative_to(arena).as_posix()] = {"sha256": sha256_of(p), "size": p.stat().st_size,
                                                        "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(p.stat().st_mtime))}
        except OSError:
            pass
    # primary CIF: the one the agent linked in its final message, else a top-level deliverable, else 'final', else newest
    pick = {}
    cifs = [p for p in cands if p.suffix.lower() == ".cif"]
    linked = None
    fm = CTRL / "runs" / arena.name / "final_message.txt"
    if fm.exists():
        txt = fm.read_text(encoding="utf-8", errors="replace").replace("\\", "/")
        for p in cifs:
            if p.relative_to(arena).as_posix() in txt or p.name in txt and txt.count(p.name) == 1 and len(cifs) == 1:
                linked = p
                break
        if linked is None:
            for p in cifs:
                if p.relative_to(arena).as_posix() in txt:
                    linked = p
                    break
    if cifs:
        def cif_key(p):
            depth = len(p.relative_to(arena).parts)
            return (p != linked, "final" not in p.name.lower(), depth, "raw" in p.name.lower() or "sensitivity" in p.relative_to(arena).as_posix().lower(), -p.stat().st_mtime)
        cifs.sort(key=cif_key)
        pick[".cif"] = cifs[0]
    anchor = pick.get(".cif")
    # companions: same directory and same stem as the primary CIF first, then same directory, then anywhere
    for ext in (".res", ".ins", ".hkl", ".fcf", ".fab"):
        best = [p for p in cands if p.suffix.lower() == ext]
        if not best:
            continue
        def comp_key(p):
            same_dir = anchor is not None and p.parent == anchor.parent
            same_stem = anchor is not None and p.stem.lower() == anchor.stem.lower()
            return (not (same_dir and same_stem), not same_dir, "final" not in p.name.lower(),
                    "sensitivity" in p.relative_to(arena).as_posix().lower(), -p.stat().st_mtime)
        best.sort(key=comp_key)
        pick[ext] = best[0]
    return {"formal_delivery": formal, "files": files, "picked": {k: str(v.relative_to(arena)) for k, v in pick.items()},
            "primary_cif_linked_in_final_message": bool(linked), "_paths": pick}


# ------------------------------------------------------------------ V1 parse
def parse_cif(path: Path) -> dict:
    code = r'''
import json, sys, math
import gemmi
p = sys.argv[1]
doc = gemmi.cif.read_file(p)
out = {"blocks": []}
for b in doc:
    d = {"name": b.name}
    def g(*tags):
        for t in tags:
            v = b.find_value(t)
            if v is not None and v not in ("?", "."):
                return gemmi.cif.as_string(v) if v.startswith(("'", '"')) or not any(ch.isdigit() for ch in v) else v
        return None
    for key, tags in {
        "space_group": ("_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"),
        "space_group_number": ("_space_group_IT_number", "_symmetry_Int_Tables_number"),
        "a": ("_cell_length_a",), "b": ("_cell_length_b",), "c": ("_cell_length_c",),
        "alpha": ("_cell_angle_alpha",), "beta": ("_cell_angle_beta",), "gamma": ("_cell_angle_gamma",), "volume": ("_cell_volume",),
        "formula_sum": ("_chemical_formula_sum",), "formula_moiety": ("_chemical_formula_moiety",), "Z": ("_cell_formula_units_Z",),
        "wavelength": ("_diffrn_radiation_wavelength",), "theta_max": ("_diffrn_reflns_theta_max",), "theta_full": ("_diffrn_reflns_theta_full",),
        "reflns_total": ("_diffrn_reflns_number",), "reflns_unique": ("_reflns_number_total",), "reflns_gt": ("_reflns_number_gt",),
        "R_int": ("_diffrn_reflns_av_R_equivalents",), "completeness": ("_diffrn_measured_fraction_theta_max",),
        "R1_gt": ("_refine_ls_R_factor_gt",), "R1_all": ("_refine_ls_R_factor_all",), "wR2": ("_refine_ls_wR_factor_ref",), "GooF": ("_refine_ls_goodness_of_fit_ref",),
        "n_reflns_ls": ("_refine_ls_number_reflns",), "n_params": ("_refine_ls_number_parameters",), "n_restraints": ("_refine_ls_number_restraints",),
        "H_treatment": ("_refine_ls_hydrogen_treatment",), "flack": ("_refine_ls_abs_structure_Flack",), "abs_details": ("_refine_ls_abs_structure_details",),
        "weighting": ("_refine_ls_weighting_details",), "max_peak": ("_refine_diff_density_max",), "min_peak": ("_refine_diff_density_min",),
        "squeeze_voids": ("_platon_squeeze_void_nr",), "smtbx_voids": ("_smtbx_masks_void_nr",), "special_details": ("_refine_special_details",),
        "absorpt": ("_exptl_absorpt_correction_type",), "computing_solution": ("_computing_structure_solution",), "computing_refinement": ("_computing_structure_refinement",),
    }.items():
        d[key] = g(*tags)
    atoms = []
    t = b.find(["_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"])
    occ = b.find_loop("_atom_site_occupancy"); adp = b.find_loop("_atom_site_adp_type"); uiso = b.find_loop("_atom_site_U_iso_or_equiv")
    for i, row in enumerate(t):
        atoms.append({"label": row[0], "type": row[1], "occ": (occ[i] if occ and i < len(occ) else None),
                      "adp": (adp[i] if adp and i < len(adp) else None), "Ueq": (uiso[i] if uiso and i < len(uiso) else None)})
    d["n_atoms"] = len(atoms); d["elements"] = sorted({a["type"].rstrip("+-0123456789") for a in atoms})
    d["n_anisotropic"] = sum(1 for a in atoms if (a["adp"] or "").lower() == "uani")
    d["n_hydrogen"] = sum(1 for a in atoms if a["type"].upper().startswith("H") and not a["type"].upper().startswith("HG") and not a["type"].upper().startswith("HF") and not a["type"].upper().startswith("HO"))
    def num(s):
        try: return float(str(s).split("(")[0])
        except Exception: return None
    d["occupancy_out_of_range"] = [a["label"] for a in atoms if a["occ"] is not None and (num(a["occ"]) is None or not (0 < num(a["occ"]) <= 1.0001))]
    aniso = b.find(["_atom_site_aniso_label", "_atom_site_aniso_U_11", "_atom_site_aniso_U_22", "_atom_site_aniso_U_33", "_atom_site_aniso_U_23", "_atom_site_aniso_U_13", "_atom_site_aniso_U_12"])
    npd = []
    for row in aniso:
        try:
            u11,u22,u33,u23,u13,u12 = [num(row[i]) for i in range(1,7)]
            m = [[u11,u12,u13],[u12,u22,u23],[u13,u23,u33]]
            # positive definiteness via leading principal minors (Sylvester)
            d1 = u11; d2 = u11*u22 - u12*u12
            d3 = (u11*(u22*u33-u23*u23) - u12*(u12*u33-u23*u13) + u13*(u12*u23-u22*u13))
            if not (d1 > 0 and d2 > 0 and d3 > 0): npd.append(row[0])
        except Exception: pass
    d["npd_adp"] = npd; d["n_aniso_rows"] = len(aniso)
    if d.get("wavelength") and d.get("theta_max"):
        try: d["d_min"] = round(num(d["wavelength"]) / (2*math.sin(math.radians(num(d["theta_max"])))), 3)
        except Exception: d["d_min"] = None
    # shortest non-H distance via gemmi structure
    try:
        # small-molecule CIF: all unit-cell sites (symmetry expanded), shortest non-H distances under 1 A
        ss = gemmi.read_small_structure(p)
        sites = [s for s in ss.get_all_unit_cell_sites() if s.element.name != "H"]
        cell = ss.cell
        short = []
        for i in range(len(sites)):
            pi = cell.orthogonalize(sites[i].fract)
            for j in range(i + 1, len(sites)):
                pj = cell.orthogonalize(sites[j].fract)
                dist = cell.find_nearest_image(pi, pj).dist()
                if dist < 1.0 and not (dist < 0.02 and sites[i].label == sites[j].label):
                    short.append([sites[i].label, sites[j].label, round(dist, 3)])
        uniq = {}
        for a, b_, dd in short:
            uniq.setdefault((tuple(sorted([a, b_])), round(dd, 2)), [a, b_, dd])
        d["contacts_below_1A_nonH"] = sorted(uniq.values(), key=lambda x: x[2])[:20]
        d["n_unit_cell_sites_nonH"] = len(sites)
    except Exception as e:
        d["contacts_below_1A_nonH"] = f"unavailable: {type(e).__name__}: {str(e)[:80]}"
    out["blocks"].append(d)
print(json.dumps(out, ensure_ascii=False))
'''
    r = subprocess.run([str(GEMMI_PY), "-X", "utf8", "-c", code, str(path)], capture_output=True, text=True, timeout=300, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        return {"error": r.stderr[-800:]}
    try:
        return json.loads(r.stdout)
    except ValueError:
        return {"error": "unparseable output", "stdout": r.stdout[-500:]}


def parse_shelx(path: Path) -> dict:
    d = {"cards": {}, "atoms": 0, "sfac": None, "cell": None, "latt": None, "symm": 0, "hklf": None, "ls": None, "flags": {}}
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as e:
        return {"error": str(e)}
    in_atoms = False
    sfac_syms = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("REM") or line.startswith("!"):
            continue
        tok = line.split()
        key = tok[0].upper()
        if key == "SFAC":
            sfac_syms += tok[1:] if not any(t.replace('.', '', 1).replace('-', '', 1).isdigit() for t in tok[1:2]) else [tok[1]]
            d["sfac"] = sfac_syms
        if key == "CELL":
            d["cell"] = tok[1:]
        if key == "LATT":
            d["latt"] = tok[1]
        if key == "SYMM":
            d["symm"] += 1
        if key == "HKLF":
            d["hklf"] = tok[1:]
        if key in ("L.S.", "CGLS"):
            d["ls"] = " ".join(tok)
        if key in ("TWIN", "BASF", "ABIN", "FAB", "OMIT", "SHEL", "WGHT", "FVAR", "PART", "EADP", "SIMU", "DELU", "RIGU", "ISOR", "DFIX", "SADI", "FLAT", "AFIX", "EXYZ", "ANIS", "MERG", "ACTA", "TEMP", "SIZE", "LIST", "EXTI", "SWAT", "HTAB", "BOND", "CONF", "PLAN"):
            d["cards"][key] = d["cards"].get(key, 0) + 1
        # an atom line: name, sfac index, x, y, z, ...
        if len(tok) >= 5 and key not in ("CELL", "ZERR", "UNIT", "SFAC", "DISP", "SYMM", "LATT", "HKLF", "WGHT", "FVAR", "OMIT", "SHEL", "TWIN", "BASF", "EXTI", "SWAT", "TEMP", "SIZE", "L.S.", "CGLS", "MERG", "FMAP", "PLAN", "BOND", "CONF", "HTAB", "LIST", "ACTA", "ANIS", "AFIX", "PART", "EADP", "SIMU", "DELU", "RIGU", "ISOR", "DFIX", "SADI", "FLAT", "EXYZ", "ABIN", "FAB", "TITL", "REM", "MORE", "END", "STIR", "DAMP", "BLOC", "FREE", "BIND", "MOVE", "RESI", "SAME", "CHIV", "DANG", "BUMP", "NCSY", "SUMP", "DEFS", "SPEC", "GRID", "MPLA", "RTAB", "WPDB", "HFIX", "EQIV", "LAUE", "HOPE", "XNPD", "PRIG", "NEUT", "SHEL", "ZERR"):
            try:
                int(tok[1]); float(tok[2]); float(tok[3]); float(tok[4])
                d["atoms"] += 1
            except ValueError:
                pass
    d["flags"] = {"twin": "TWIN" in d["cards"] or "BASF" in d["cards"], "mask_abin": "ABIN" in d["cards"], "omit": d["cards"].get("OMIT", 0),
                  "shel": "SHEL" in d["cards"], "restraints": sum(d["cards"].get(k, 0) for k in ("SIMU", "DELU", "RIGU", "ISOR", "DFIX", "SADI", "FLAT", "EADP")),
                  "parts": d["cards"].get("PART", 0)}
    return d


def hkl_stats(path: Path) -> dict:
    n = 0; nonzero = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if len(line) < 20:
                    continue
                try:
                    h, k, l = int(line[0:4]), int(line[4:8]), int(line[8:12])
                except ValueError:
                    continue
                n += 1
                if (h, k, l) != (0, 0, 0):
                    nonzero += 1
                else:
                    break
    except OSError as e:
        return {"error": str(e)}
    return {"reflection_lines": nonzero, "sha256": sha256_of(path), "size": path.stat().st_size}


# ------------------------------------------------------------------ V3 fixed-model recompute
LST_R1 = re.compile(r"R1\s*=\s*([0-9.]+)\s+for\s+(\d+)\s+Fo\s*>\s*4sig\(Fo\)\s+and\s+([0-9.]+)\s+for\s+all\s+(\d+)\s+data", re.I)
LST_WR2 = re.compile(r"wR2\s*=\s*([0-9.]+),\s*GooF\s*=\s*S\s*=\s*([0-9.]+)", re.I)


def shelxl_zero_cycle(work: Path, res: Path, hkl: Path, fab: Path | None, tag: str, drop_abin: bool = False, timeout: int = 1800) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    stem = f"verify_{tag}"
    text = res.read_text(encoding="utf-8", errors="replace").splitlines()
    out_lines = []
    for line in text:
        key = line.strip().split()[0].upper() if line.strip() else ""
        if key in ("L.S.", "CGLS"):
            out_lines.append("L.S. 0")
            continue
        if key == "ACTA":
            continue
        if drop_abin and key == "ABIN":
            continue
        if key == "END":
            out_lines.append("END")
            break
        out_lines.append(line)
    if not any(l.strip().upper().startswith("L.S. 0") for l in out_lines):
        # no refinement card at all: insert before the first atom-ish or before HKLF
        idx = next((i for i, l in enumerate(out_lines) if l.strip().upper().startswith("HKLF")), len(out_lines))
        out_lines.insert(idx, "L.S. 0")
    (work / f"{stem}.ins").write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    shutil.copy2(hkl, work / f"{stem}.hkl")
    if fab and fab.exists() and not drop_abin:
        shutil.copy2(fab, work / f"{stem}.fab")
    t0 = time.time()
    rc = run_engine([str(SHELXL), stem], work, work / f"{stem}.stdout.txt", timeout=timeout)
    lst = work / f"{stem}.lst"
    rec = {"tag": tag, "exit_code": rc, "seconds": round(time.time() - t0, 1), "ins_sha256": sha256_of(work / f"{stem}.ins"), "hkl_sha256": sha256_of(work / f"{stem}.hkl"),
           "drop_abin": drop_abin, "R1_gt": None, "n_gt": None, "R1_all": None, "n_all": None, "wR2": None, "GooF": None}
    if lst.exists():
        t = lst.read_text(encoding="utf-8", errors="replace")
        m = LST_R1.findall(t)
        if m:
            r1g, ng, r1a, na = m[-1]
            rec.update({"R1_gt": float(r1g), "n_gt": int(ng), "R1_all": float(r1a), "n_all": int(na)})
        w = LST_WR2.findall(t)
        if w:
            rec.update({"wR2": float(w[-1][0]), "GooF": float(w[-1][1])})
        err = re.findall(r"\*\*\s.*", t)
        rec["shelxl_warnings"] = [e.strip()[:120] for e in err[:10]]
    return rec


# ------------------------------------------------------------------ V4 PLATON checkCIF
def run_engine(argv: list, cwd: Path, log: Path, timeout: int, done_marker: str | None = None, marker_file: Path | None = None,
               env: dict | None = None):
    """Run a science program with stdout/stderr on a file (never a pipe: PLATON keeps running after it has
    written its report and would hold a pipe open forever). Ends when the process exits, when `done_marker`
    appears in `marker_file`, or at `timeout`; in the last two cases the process tree is killed."""
    import psutil
    t0 = time.time()
    with log.open("w", encoding="utf-8") as fh:
        proc = subprocess.Popen(argv, cwd=str(cwd), stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=env)
        rc = None
        while True:
            rc = proc.poll()
            if rc is not None:
                break
            if done_marker and marker_file and marker_file.exists():
                try:
                    if done_marker in marker_file.read_text(encoding="utf-8", errors="replace"):
                        rc = "done_marker"
                        break
                except OSError:
                    pass
            if time.time() - t0 > timeout:
                rc = "timeout"
                break
            time.sleep(1.0)
        if rc in ("done_marker", "timeout"):
            time.sleep(2.0)
            try:
                root = psutil.Process(proc.pid)
                for p in root.children(recursive=True) + [root]:
                    try:
                        p.kill()
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
            except psutil.NoSuchProcess:
                pass
    return rc


def platon_checkcif(work: Path, cif: Path, timeout: int = 900) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    target = work / "verify.cif"
    shutil.copy2(cif, target)
    env = dict(os.environ)
    # PLATON recreates the FCF with SHELXL and needs it on PATH or in SHLEXE
    env["SHLEXE"] = str(SHELXL)
    env["PATH"] = str(SHELXL.parent) + ";" + env.get("PATH", "")
    checkdef = ROOT / "toolbox" / "shelx" / "check.def"
    if checkdef.exists():
        env["CHECKDEF"] = str(checkdef)
    t0 = time.time()
    rec = {"exit_code": None, "seconds": None, "alerts": {"A": [], "B": [], "C": [], "G": []}, "chk_file": None}
    rec["exit_code"] = run_engine([str(PLATON), "-u", "verify.cif"], work, work / "platon.stdout.txt", timeout=timeout,
                                  done_marker="CheckCIF out on", marker_file=work / "platon.out", env=env)
    rec["seconds"] = round(time.time() - t0, 1)
    if (work / "platon.out").exists():
        rec["platon_out_tail"] = (work / "platon.out").read_text(encoding="utf-8", errors="replace").strip().splitlines()[-8:]
        rec["fcf_recreated"] = "PROBLEM to Recreate FCF" not in "\n".join(rec["platon_out_tail"])
    chk = work / "verify.chk"
    if chk.exists():
        rec["chk_file"] = str(chk)
        for line in chk.read_text(encoding="utf-8", errors="replace").splitlines():
            m = re.match(r"\s*(?:PLAT)?(\d{3})_ALERT_(\d)_([ABCG])\s+(.*)", line)
            if m:
                rec["alerts"][m.group(3)].append({"code": "PLAT" + m.group(1), "type": int(m.group(2)), "text": m.group(4).strip()[:120]})
    rec["counts"] = {k: len(v) for k, v in rec["alerts"].items()}
    return rec


# ------------------------------------------------------------------ V5 claims
def extract_claims(texts: list) -> dict:
    blob = "\n".join(t for t in texts if t)
    r1s = [float(x) for x in re.findall(r"R1\s*[=:：≈~]\s*0?\.(\d{3,4})", blob) for x in [("0." + x)]]
    sgs = re.findall(r"\b(P\s?-?1|P\s?2[1₁]\s?2[1₁]\s?2[1₁]|P\s?2[1₁]/[cn]|P\s?2[1₁]|C\s?2/c|P\s?6/m\s?m\s?m|P\s?6/mmm|P\s?-?3\w*|Pnma|Pbca|R\s?-?3\w*|I\s?4/m\w*|Fm-?3m|Im-?3m|P\s?4/m\w*)\b", blob)
    unresolved = any(k in blob for k in ("未解决", "未完成", "尚未", "unresolved", "not resolved", "remaining issues", "尚存", "未能"))
    # a claim of publication readiness counts only when it is not negated in the same sentence
    positive = re.findall(r"[^。.\n]*(?:发表质量|发表级|publication[- ]quality|publication[- ]ready)[^。.\n]*", blob, re.I)
    negations = ("尚未", "未达", "没有达到", "未能", "不到", "距离", "not yet", "not reach", "below", "short of", "not publication")
    completed_claim = any(not any(n in s for n in negations) for s in positive)
    return {"R1_claims": sorted(set(r1s)), "space_group_mentions": sorted(set(s.replace(" ", "") for s in sgs)), "mentions_unresolved": unresolved,
            "claims_completion": completed_claim, "publication_sentences": [s.strip()[:160] for s in positive][:6], "chars": len(blob)}


# ------------------------------------------------------------------ V6 status & rubric
def decide(rep: dict) -> dict:
    run = rep["run"]
    exit_reason = run.get("exit_reason")
    cif = (rep["parsed"].get("cif") or {}).get("blocks", [{}])[0] if rep["parsed"].get("cif") else {}
    res = rep["parsed"].get("res") or {}
    has_model = bool((cif.get("n_atoms") or 0) > 0 or (res.get("atoms") or 0) > 0)
    recomp = rep.get("recompute", {}).get("main") or {}
    recomputed = recomp.get("R1_gt") is not None
    reported = None
    try:
        reported = float(str(cif.get("R1_gt")).split("(")[0]) if cif.get("R1_gt") else None
    except ValueError:
        reported = None
    r1_diff = abs(recomp["R1_gt"] - reported) if (recomputed and reported is not None) else None
    platon = rep.get("checkcif") or {}
    # level-A alerts that only report missing CIF metadata (cell-measurement items, temperature, radiation
    # type, crystal description, transmission factors) cannot make a structure wrong; the others can
    metadata_codes = {"PLAT058", "PLAT059", "PLAT183", "PLAT184", "PLAT185", "PLAT186", "PLAT187", "PLAT188", "PLAT189", "PLAT190",
                      "PLAT191", "PLAT192", "PLAT193", "PLAT194", "PLAT195", "PLAT196", "PLAT197", "PLAT198", "PLAT199", "PLAT660",
                      "PLAT699", "PLAT700", "PLAT701", "PLAT702", "PLAT703", "PLAT704", "PLAT705", "PLAT706", "PLAT707", "PLAT708",
                      "PLAT709", "PLAT710", "PLAT711", "PLAT712", "PLAT713", "PLAT714", "PLAT715", "PLAT716", "PLAT717", "PLAT718",
                      "PLAT719", "PLAT720", "PLAT721", "PLAT722", "PLAT723", "PLAT724", "PLAT725", "PLAT726", "PLAT727", "PLAT728",
                      "PLAT729", "PLAT730", "PLAT731", "PLAT732", "PLAT733", "PLAT734", "PLAT735", "PLAT736", "PLAT737", "PLAT738",
                      "PLAT739", "PLAT740", "PLAT741", "PLAT742", "PLAT743", "PLAT744", "PLAT745", "PLAT746", "PLAT747", "PLAT748",
                      "PLAT749", "PLAT750", "PLAT751", "PLAT752", "PLAT753", "PLAT754", "PLAT755", "PLAT756", "PLAT757", "PLAT758",
                      "PLAT759", "PLAT760", "PLAT761", "PLAT762", "PLAT763", "PLAT764", "PLAT765", "PLAT766", "PLAT767", "PLAT768",
                      "PLAT769", "PLAT770", "PLAT771", "PLAT772", "PLAT773", "PLAT774", "PLAT775", "PLAT776", "PLAT777", "PLAT778",
                      "PLAT779", "PLAT780", "PLAT781", "PLAT782", "PLAT783", "PLAT784", "PLAT785", "PLAT786", "PLAT787", "PLAT788",
                      "PLAT789", "PLAT790", "PLAT791", "PLAT792", "PLAT793", "PLAT794", "PLAT795", "PLAT796", "PLAT797", "PLAT798",
                      "PLAT799", "PLAT911", "PLAT912", "PLAT913", "PLAT914", "PLAT915", "PLAT916", "PLAT917", "PLAT918", "PLAT919",
                      "PLAT960", "PLAT961", "PLAT962", "PLAT963", "PLAT964", "PLAT965", "PLAT966", "PLAT967", "PLAT968", "PLAT969"}
    a_list = (platon.get("alerts") or {}).get("A") or []
    a_meta = [a for a in a_list if a.get("code") in metadata_codes]
    a_subst = [a for a in a_list if a.get("code") not in metadata_codes]
    a_alerts = len(a_subst)
    rep["checkcif_a_classification"] = {"metadata": [a["code"] for a in a_meta], "substantive": [a["code"] + " " + a["text"][:60] for a in a_subst]}
    fatal = []
    if cif.get("npd_adp"):
        fatal.append(f"non-positive-definite ADPs: {cif['npd_adp'][:5]}")
    if cif.get("occupancy_out_of_range"):
        fatal.append(f"occupancies out of range: {cif['occupancy_out_of_range'][:5]}")
    if isinstance(cif.get("contacts_below_1A_nonH"), list) and cif["contacts_below_1A_nonH"]:
        fatal.append(f"non-H contacts below 1 A: {cif['contacts_below_1A_nonH'][:3]}")
    if r1_diff is not None and r1_diff > R1_TOL:
        fatal.append(f"fixed-model R1 differs from reported by {r1_diff:.4f}")
    if a_alerts:
        fatal.append(f"{a_alerts} substantive level-A checkCIF alerts: " + "; ".join(a["code"] + " " + a["text"][:45] for a in a_subst[:4]))
    if exit_reason in ("preflight_blocked", "infra_failed", "channel_closed"):
        status = "not_evaluable"
    elif not has_model:
        status = "no_usable_result"
    elif recomputed and platon.get("chk_file") and not fatal:
        status = "candidate_checks_passed"
    elif has_model and (rep["claims"].get("mentions_unresolved") or fatal):
        status = "diagnostic_delivery" if rep["claims"].get("mentions_unresolved") and rep["delivery"].get("formal_delivery") else "partial_result"
    else:
        status = "partial_result"
    # rubric: 5 dimensions x 4 items, 5 / 2.5 / 0 / null
    def item(cond_pass, cond_partial=False, unknown=False):
        if unknown:
            return None
        return 5.0 if cond_pass else (2.5 if cond_partial else 0.0)
    inputs_ok = not run.get("inputs_tampered")
    chain = rep.get("chain", {})
    rubric = {
        "数据与来源": {
            "输入来源完整且一致": item(inputs_ok),
            "处理链/依赖可执行": item(chain.get("dials_outputs_present") is True or rep["run"]["case_id"] == "nu1000" and has_model, cond_partial=bool(chain.get("any_logs"))),
            "数据质量与取舍有记录": item(bool(cif.get("R_int") or chain.get("scaling_log")), cond_partial=has_model, unknown=not has_model),
            "指标与反射集合定义清楚": item(bool(cif.get("n_reflns_ls") and cif.get("theta_max")), cond_partial=bool(recomp.get("n_all")), unknown=not has_model),
        },
        "结构与化学": {
            "元素及占有率有支持": item(False, unknown=True),
            "连接/几何合理": item(False, cond_partial=False, unknown=True) if not fatal else 0.0,
            "对称性/无序/ADP 判断有支持": item(bool(has_model and not cif.get("npd_adp") and (cif.get("n_anisotropic") or 0) > 0), cond_partial=has_model, unknown=not has_model),
            "重要缺失密度与模型不确定性得到检查": item(False, unknown=True),
        },
        "可复现性": {
            "产物可解析": item(bool(cif.get("n_atoms")) and "error" not in (rep["parsed"].get("cif") or {}), cond_partial=has_model),
            "最终结构与反射/参数匹配": item(bool(rep["delivery"]["picked"].get(".hkl") and rep["delivery"]["picked"].get(".res")), cond_partial=has_model),
            "固定模型复算成立": item(recomputed and r1_diff is not None and r1_diff <= R1_TOL, cond_partial=recomputed, unknown=not has_model),
            "关键处理链复跑成立": item(False, cond_partial=chain.get("dials_outputs_present") is True, unknown=rep["run"]["case_id"] == "nu1000" and has_model),
        },
        "证据与判断": {
            "关键模型变更有具体证据": item(False, unknown=True),
            "竞争假设/失败尝试如实留存": item(False, unknown=True),
            "报告数值与最终模型一致": item(bool(rep["claims"]["R1_claims"]) and reported is not None and any(abs(c - reported) < 0.0015 for c in rep["claims"]["R1_claims"]),
                                     cond_partial=(reported is not None and not rep["claims"]["R1_claims"]), unknown=not has_model),
            "对完成程度与未知内容的声明准确": item(bool(rep["claims"]["mentions_unresolved"]) and not (status == "no_usable_result" and rep["claims"]["claims_completion"]),
                                              cond_partial=not rep["claims"]["claims_completion"]),
        },
        "交付完整性": {
            "可用 CIF/相应精修文件": item(bool(rep["delivery"]["picked"].get(".cif") and rep["delivery"]["picked"].get(".res")), cond_partial=has_model),
            "完整校验文件": item(bool(rep["delivery"].get("has_validation")), cond_partial=bool(platon.get("chk_file"))),
            "重跑入口与版本": item(bool(rep["delivery"].get("has_readme_or_summary") and rep["delivery"].get("has_logs")), cond_partial=bool(rep["delivery"].get("has_logs"))),
            "清楚列出未解决问题和任务范围": item(bool(rep["claims"]["mentions_unresolved"]), cond_partial=bool(rep["claims"]["chars"] > 200)),
        },
    }
    confirmed = sum(v for dim in rubric.values() for v in dim.values() if v is not None)
    unknown_items = sum(1 for dim in rubric.values() for v in dim.values() if v is None)
    return {"scientific_status": status, "publication_readiness": "not_expert_reviewed", "ground_truth_correctness": "unknown",
            "fatal_findings": fatal, "r1_reported": reported, "r1_recomputed": recomp.get("R1_gt"), "r1_abs_difference": r1_diff,
            "r1_definition": "R1 for Fo > 4sig(Fo), SHELXL L.S. 0 on the delivered RES with the delivered HKL" if recomputed else None,
            "rubric": rubric, "score_lower_bound": confirmed, "score_upper_bound": confirmed + 5.0 * unknown_items, "unknown_items": unknown_items}


# ------------------------------------------------------------------ main
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cpucap  # noqa: E402

def main() -> int:
    print("cpu cap:", cpucap.limit_self(), flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--skip-platon", action="store_true")
    a = ap.parse_args()
    evid = CTRL / "runs" / a.run_id
    manifest = json.loads((evid / "manifest.json").read_text(encoding="utf-8"))
    arena = Path(manifest["arena"])
    vdir = CTRL / "verifier" / a.run_id
    if vdir.exists():
        raise SystemExit(f"verifier dir exists for {a.run_id}; refusing to overwrite")
    vdir.mkdir(parents=True)
    rep = {"run_id": a.run_id, "verified_at": now(), "run": manifest, "tolerance_R1": R1_TOL}
    cond, case = manifest["condition"], manifest["case_id"]
    # V0
    d = locate_delivery(arena, cond)
    paths = d.pop("_paths")
    sealed = vdir / "sealed_delivery"
    sealed.mkdir()
    for rel in d["files"]:
        src = arena / rel
        dst = sealed / Path(rel).name
        if dst.exists():
            dst = sealed / (Path(rel).parent.name + "__" + Path(rel).name)
        shutil.copy2(src, dst)
    d["has_validation"] = any(n.lower().endswith((".chk", "validation.md", "checkcif.html", ".fcf")) for n in d["files"])
    d["has_readme_or_summary"] = any(Path(n).name.lower() in ("summary.md", "readme.md", "report.json", "report.md", "notes.md") for n in d["files"]) or \
        any(p.name.lower() in ("summary.md", "readme.md", "report.md") for p in arena.rglob("*.md") if "inputs" not in p.parts)
    d["has_logs"] = any(n.lower().endswith((".lst", ".log", ".txt")) for n in d["files"]) or any(arena.rglob("*.lst")) or any(arena.rglob("*.log"))
    rep["delivery"] = d
    # V1
    parsed = {}
    if ".cif" in paths:
        parsed["cif"] = parse_cif(paths[".cif"])
    for ext in (".res", ".ins"):
        if ext in paths:
            parsed[ext[1:]] = parse_shelx(paths[ext])
    if ".hkl" in paths:
        parsed["hkl"] = hkl_stats(paths[".hkl"])
    rep["parsed"] = parsed
    # V2 chain
    chain = {}
    if case == "alanine":
        expts = [p for p in arena.rglob("*.expt") if "inputs" not in p.parts]
        refls = [p for p in arena.rglob("*.refl") if "inputs" not in p.parts]
        logs = [p for p in arena.rglob("dials.*.log") if "inputs" not in p.parts]
        chain = {"dials_outputs_present": bool(expts and refls), "n_expt": len(expts), "n_refl": len(refls), "any_logs": bool(logs),
                 "dials_logs": sorted(p.name for p in logs)[:30], "scaling_log": any(p.name in ("dials.scale.log",) for p in logs),
                 "export_log": any(p.name == "dials.export.log" for p in logs),
                 "note": "presence of the raw-frame processing chain is recorded; a full replay from frames was not executed (partial replay)"}
    else:
        start = json.loads((ROOT / "staging" / "nu1000" / "manifest.json").read_text(encoding="utf-8"))
        chain = {"start_hkl_sha256": start["start.hkl"]["sha256"], "final_hkl_sha256": parsed.get("hkl", {}).get("sha256"),
                 "final_hkl_identical_to_start": parsed.get("hkl", {}).get("sha256") == start["start.hkl"]["sha256"],
                 "final_hkl_reflection_lines": parsed.get("hkl", {}).get("reflection_lines"), "start_hkl_lines": 337816 - 1}
    rep["chain"] = chain
    # V3 recompute (needs res + hkl)
    recompute = {}
    if ".res" in paths and ".hkl" in paths and (parsed.get("res", {}).get("atoms") or 0) > 0:
        fab = paths.get(".fab")
        recompute["main"] = shelxl_zero_cycle(vdir / "recompute", paths[".res"], paths[".hkl"], fab, "main")
        if parsed.get("res", {}).get("flags", {}).get("mask_abin"):
            recompute["no_mask"] = shelxl_zero_cycle(vdir / "recompute", paths[".res"], paths[".hkl"], fab, "nomask", drop_abin=True)
    else:
        recompute["note"] = "no RES+HKL pair with atoms; fixed-model recompute not possible"
    rep["recompute"] = recompute
    # V4 checkCIF
    if ".cif" in paths and not a.skip_platon:
        rep["checkcif"] = platon_checkcif(vdir / "checkcif", paths[".cif"])
    else:
        rep["checkcif"] = {"note": "no CIF or PLATON skipped"}
    # V5 claims
    texts = []
    for name in ("final_message.txt",):
        p = evid / name
        if p.exists():
            texts.append(p.read_text(encoding="utf-8", errors="replace"))
    for p in list(arena.rglob("SUMMARY.md")) + list(arena.rglob("summary.md")) + list(arena.rglob("README.md")):
        if "inputs" not in p.parts:
            texts.append(p.read_text(encoding="utf-8", errors="replace")[:20000])
    rep["claims"] = extract_claims(texts)
    # V6
    rep["verdict"] = decide(rep)
    (vdir / "report.json").write_text(json.dumps(rep, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    v = rep["verdict"]
    print(json.dumps({"run_id": a.run_id, "status": v["scientific_status"], "r1_reported": v["r1_reported"], "r1_recomputed": v["r1_recomputed"],
                      "diff": v["r1_abs_difference"], "fatal": v["fatal_findings"], "score": [v["score_lower_bound"], v["score_upper_bound"]],
                      "checkcif_counts": (rep["checkcif"] or {}).get("counts")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
