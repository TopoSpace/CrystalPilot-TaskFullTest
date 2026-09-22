"""Independent structure-factor recompute for runs that delivered a CIF plus reflections but no SHELX RES/INS
(the pure baselines wrote their own refinement code). The model in the CIF and the agent's own reflection file
are converted to a SHELXL job in which every atomic parameter is fixed (SHELXL's +10 convention) and only the
overall scale factor is refined: R1 and wR2 of the delivered model on the delivered data, computed by an
independent program. (SHELXL does not optimise the scale in an L.S. 0 job; a wrong starting scale inflated R1
from 0.18 to 0.44 on a control model, so the scale must be refined while the model stays fixed.) Validated on the
r04 SHELXL model: 0.1816 vs the RES-based zero-cycle value 0.1817.

    recompute_from_cif.py --run-id <rid>            (uses the verifier's picked .cif / .hkl)
    recompute_from_cif.py --cif <file> --hkl <file> --out <dir>

Two variants: MERG 2 (SHELXL default: Friedel pairs kept for non-centrosymmetric structures) and MERG 4 (Friedel
pairs merged, f'' = 0), the latter matching self-written code that ignores anomalous scattering. Results go to
control/verifier/<rid>/recompute_from_cif/ and recompute_from_cif.json; report.json gets a `recompute_from_cif` key.
Weights: the CIF weighting scheme (a, b) when parseable, else WGHT 0.1; R1 does not depend on it."""
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
SHELXL = ROOT / "toolbox" / "shelx" / "shelxl.exe"
LATT = {"P": 1, "I": 2, "R": 3, "F": 4, "A": 5, "B": 6, "C": 7}


def read_cif(cif: Path) -> dict:
    import gemmi
    st = gemmi.read_small_structure(str(cif))
    doc = gemmi.cif.read(str(cif))
    block = doc.sole_block()
    wl = block.find_value("_diffrn_radiation_wavelength")
    wdet = block.find_value("_refine_ls_weighting_details") or ""
    sg = None
    hm = st.spacegroup_hm or block.find_value("_space_group_name_H-M_alt") or block.find_value("_symmetry_space_group_name_H-M")
    if hm:
        sg = gemmi.find_spacegroup_by_name(hm.strip("'\" "))
    if sg is None:
        ops = [v.strip("'\" ") for v in block.find_values("_space_group_symop_operation_xyz")] or \
              [v.strip("'\" ") for v in block.find_values("_symmetry_equiv_pos_as_xyz")]
        if ops:
            sg = gemmi.find_spacegroup_by_ops(gemmi.GroupOps([gemmi.Op(o) for o in ops]))
    if sg is None:
        raise SystemExit("space group not resolvable from the CIF")
    sites = []
    for s in st.sites:
        el = s.element.name if s.element else (s.type_symbol or s.label)[:2].rstrip("0123456789+-")
        an = None
        if s.aniso.nonzero():
            a = s.aniso
            an = (a.u11, a.u22, a.u33, a.u23, a.u13, a.u12)
        sites.append({"label": s.label, "element": el.capitalize(), "x": s.fract.x, "y": s.fract.y, "z": s.fract.z,
                      "occ": s.occ if s.occ else 1.0, "uiso": s.u_iso if s.u_iso else 0.05, "aniso": an})
    return {"cell": [st.cell.a, st.cell.b, st.cell.c, st.cell.alpha, st.cell.beta, st.cell.gamma], "sg": sg,
            "wavelength": float(str(wl).split("(")[0]) if wl and str(wl)[0].isdigit() else 0.71073, "sites": sites, "weighting": wdet,
            "hm": sg.hm, "n_ops_total": len(sg.operations()) if hasattr(sg.operations(), "__len__") else None}


def shelx_symm(sg) -> tuple:
    """LATT code and SYMM triplets: one representative per inversion pair when centrosymmetric."""
    import gemmi
    ops = sg.operations()
    centro = ops.is_centrosymmetric()
    letter = sg.hm.strip()[0].upper()
    latt = LATT.get(letter, 1) * (1 if centro else -1)
    keep, seen = [], set()
    ident = gemmi.Op("x,y,z")
    inv = gemmi.Op("-x,-y,-z")
    for op in ops.sym_ops:
        if op == ident:
            continue
        neg = (inv * op).wrap()
        key, nkey = op.wrap().triplet(), neg.triplet()
        if key in seen or nkey in seen:
            continue
        seen.add(key)
        if centro and neg == ident:
            continue
        keep.append(key.upper())
    return latt, keep


def site_order(sg, fract, tol: float = 0.002) -> int:
    """Number of space-group operations that map the site onto itself (site-symmetry order). SHELXL expects the
    site occupation factor of an atom on a special position to be 1/order (e.g. 10.5 on a mirror plane); the CIF
    convention is 1. Computed from the operations, so it does not depend on optional CIF tags."""
    x, y, z = fract
    n = 0
    for op in sg.operations():
        px, py, pz = op.apply_to_xyz([x, y, z])
        if all(abs(d - round(d)) < tol for d in (px - x, py - y, pz - z)):
            n += 1
    return max(1, n)


DIALS_PY = Path(r"C:\users\<user>\miniforge3\envs\dials\python.exe")


def disp_lines(elements: list, wavelength: float) -> list:
    """DISP cards (f', f'') for a non-standard wavelength from the cctbx Sasaki tables (Henke as fallback,
    hydrogen 0). SHELXL only knows Mo and Cu Kalpha by itself; at 0.689 A the Zr f' is about -9 e, without which
    a Zr framework model is evaluated wrongly."""
    if abs(wavelength - 0.71073) < 0.002 or abs(wavelength - 1.54184) < 0.002 or not DIALS_PY.exists():
        return []
    code = ("import sys\nfrom cctbx.eltbx import sasaki, henke\nwl=float(sys.argv[1])\n"
            "for el in sys.argv[2:]:\n"
            "    if el.upper()=='H':\n        print(el, 0.0, 0.0); continue\n"
            "    try:\n        t=sasaki.table(el).at_angstrom(wl); print(el, t.fp(), t.fdp())\n"
            "    except Exception:\n"
            "        try:\n            t=henke.table(el).at_angstrom(wl); print(el, t.fp(), t.fdp())\n"
            "        except Exception:\n            print(el, 'nan', 'nan')\n")
    import os
    env = dict(os.environ)
    d = DIALS_PY.parent
    env["PATH"] = ";".join([str(d / "Library" / "bin"), str(d / "Scripts"), str(d), env.get("PATH", "")])
    try:
        r = subprocess.run([str(DIALS_PY), "-c", code, f"{wavelength:.6f}", *elements], capture_output=True, text=True, timeout=120, env=env)
    except (subprocess.TimeoutExpired, OSError):
        return []
    out = []
    for line in r.stdout.splitlines():
        p = line.split()
        if len(p) == 3 and p[1] != "nan":
            out.append(f"DISP {p[0]} {float(p[1]):.5f} {float(p[2]):.5f} 0.0")   # mu placeholder: not used for Fc
    return out


def parse_weighting(text: str) -> str:
    m = re.search(r"\(\s*([0-9.]+)\s*P\s*\)\s*\^?2\^?\s*(?:\+\s*([0-9.]+)\s*P)?", text or "")
    if m:
        a = float(m.group(1)); b = float(m.group(2)) if m.group(2) else 0.0
        return f"WGHT {a:.4f} {b:.4f}"
    return "WGHT 0.1000"


def write_hkl(src: Path, dst: Path) -> dict:
    rows = []
    for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        parts = t.split()
        if len(parts) < 5:
            continue
        try:
            h, k, l = int(parts[0]), int(parts[1]), int(parts[2]); i, s = float(parts[3]), float(parts[4])
        except ValueError:
            continue
        if h == 0 and k == 0 and l == 0:
            continue
        rows.append((h, k, l, i, s))
    mx = max(max(abs(r[3]), abs(r[4])) for r in rows) if rows else 1.0
    scale = 1.0
    while mx * scale > 99999.0:
        scale /= 10.0
    with dst.open("w", encoding="ascii", newline="\n") as fh:
        for h, k, l, i, s in rows:
            fh.write(f"{h:4d}{k:4d}{l:4d}{i*scale:8.2f}{s*scale:8.2f}\n")
        fh.write("   0   0   0    0.00    0.00\n")
    return {"n_reflections": len(rows), "scale_applied": scale, "source": str(src)}


def build_ins(model: dict, stem: str, merg: int, shel: float | None = None) -> str:
    a, b, c, al, be, ga = model["cell"]
    latt, symm = shelx_symm(model["sg"])
    elements = []
    for s in model["sites"]:
        if s["element"] not in elements:
            elements.append(s["element"])
    nops = len(model["sg"].operations())
    unit = []
    for el in elements:
        unit.append(sum(s["occ"] for s in model["sites"] if s["element"] == el) * nops)
    lines = [f"TITL {stem} recompute of delivered CIF model, atoms fixed, scale refined",
             f"CELL {model['wavelength']:.5f} {a:.4f} {b:.4f} {c:.4f} {al:.3f} {be:.3f} {ga:.3f}",
             f"ZERR {nops} 0.001 0.001 0.001 0.01 0.01 0.01", f"LATT {latt}"]
    lines += [f"SYMM {s}" for s in symm]
    lines += ["SFAC " + " ".join(elements)]
    lines += disp_lines(elements, model["wavelength"])
    lines += ["UNIT " + " ".join(f"{u:.0f}" for u in unit),
              f"MERG {merg}", "L.S. 4", "LIST 4", "FMAP 2", "PLAN 5", parse_weighting(model["weighting"])]
    if shel:
        lines.append(f"SHEL 999 {shel:g}")
    lines.append("FVAR 1.0")
    used = set()

    def fx(v: float) -> float:      # SHELXL fixes a parameter when 10 is added to it
        return v + 10.0

    for s in model["sites"]:
        lab = re.sub(r"[^A-Za-z0-9]", "", s["label"])[:4] or s["element"]
        base, n = lab, 1
        while lab.upper() in used:
            lab = (base[:3] + str(n))[:4]; n += 1
        used.add(lab.upper())
        sfac = elements.index(s["element"]) + 1
        occ = 10.0 + s["occ"] / site_order(model["sg"], (s["x"], s["y"], s["z"]))
        if s["aniso"]:
            u11, u22, u33, u23, u13, u12 = s["aniso"]
            lines.append(f"{lab} {sfac} {fx(s['x']):.6f} {fx(s['y']):.6f} {fx(s['z']):.6f} {occ:.5f} {fx(u11):.5f} {fx(u22):.5f} =\n"
                         f"   {fx(u33):.5f} {fx(u23):.5f} {fx(u13):.5f} {fx(u12):.5f}")
        else:
            lines.append(f"{lab} {sfac} {fx(s['x']):.6f} {fx(s['y']):.6f} {fx(s['z']):.6f} {occ:.5f} {fx(s['uiso']):.5f}")
    lines += ["HKLF 4", "END", ""]
    return "\n".join(lines)


def _fvar(path: Path) -> float | None:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.upper().startswith("FVAR"):
            try:
                return float(line.split()[1])
            except (IndexError, ValueError):
                return None
    return None


def run_shelxl(work: Path, stem: str, timeout: int = 900, max_iter: int = 8) -> dict:
    """Scale-only refinement (all atomic parameters fixed with the +10 convention, L.S. 4). The job is repeated
    with the updated FVAR until the scale is stable, so a start far from the data scale still converges; atom
    parameters never change."""
    log = work / f"{stem}.stdout.txt"
    t0 = time.time()
    rc, iters, history = None, 0, []
    with log.open("w", encoding="utf-8") as fh:
        while iters < max_iter:
            iters += 1
            before = _fvar(work / f"{stem}.ins")
            try:
                rc = subprocess.run([str(SHELXL), stem], cwd=str(work), stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, timeout=timeout).returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
                break
            res = work / f"{stem}.res"
            if not res.exists():
                break
            after = _fvar(res)
            history.append((before, after))
            if before is None or after is None or after == 0:
                break
            if abs(after - before) / abs(before) < 1e-4:
                break
            shutil.copy2(res, work / f"{stem}.ins")
    rec = {"exit_code": rc, "seconds": round(time.time() - t0, 1), "scale_iterations": iters, "scale_history": history}
    lst = work / f"{stem}.lst"
    if lst.exists():
        t = lst.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"R1 =\s+([0-9.]+) for\s+(\d+) Fo > 4sig\(Fo\) and\s+([0-9.]+) for all\s+(\d+) data", t)
        if m:
            rec.update({"R1_gt": float(m.group(1)), "n_gt": int(m.group(2)), "R1_all": float(m.group(3)), "n_all": int(m.group(4))})
        m2 = re.search(r"wR2 =\s+([0-9.]+),\s*GooF = S =\s+([0-9.]+)", t)
        if m2:
            rec.update({"wR2": float(m2.group(1)), "GooF": float(m2.group(2))})
        m3 = re.search(r"(\d+)\s+Reflections read, of which\s+(\d+)\s+rejected", t)
        if m3:
            rec.update({"reflections_read": int(m3.group(1)), "rejected": int(m3.group(2))})
        m4 = re.search(r"(\d+)\s+Unique reflections, of which\s+(\d+)\s+suppressed", t)
        if m4:
            rec.update({"unique": int(m4.group(1)), "suppressed": int(m4.group(2))})
        warn = [ln.strip() for ln in t.splitlines() if "**" in ln][:8]
        rec["warnings"] = warn
    return rec


def update_verdict(rep: dict, result: dict) -> None:
    """When the RES-based zero-cycle recompute was impossible (no SHELX files delivered), the rubric item for the
    fixed-model recompute was scored 0. The CIF-based recompute establishes the same fact by another route, so the
    item and the R1 fields are filled from it: 5 if |dR1| <= 0.005 (the agent's own weighting and merging differ
    from SHELXL's), 2.5 if <= 0.02, else 0. Score bounds move with the item (lower = known points, upper = lower +
    5 x unknown items). Nothing else in the verdict changes."""
    v = rep.get("verdict") or {}
    best = result.get("closest_variant")
    if v.get("r1_recomputed") is not None or not best:
        return
    d = best.get("abs_difference_to_reported")
    v["r1_recomputed"] = best["R1_gt"]
    v["r1_abs_difference"] = d
    v["r1_definition"] = f"delivered CIF model on the delivered reflections, SHELXL with all atomic parameters fixed and the scale refined ({best['variant']}); the RES-based zero-cycle route was not available"
    v["recompute_source"] = "recompute_from_cif"
    new_item = 5.0 if d is not None and d <= 0.005 else (2.5 if d is not None and d <= 0.02 else 0.0)
    rub = v.get("rubric") or {}
    for group, items in rub.items():
        if isinstance(items, dict) and "固定模型复算成立" in items:
            old_item = items["固定模型复算成立"]
            items["固定模型复算成立"] = new_item
            if isinstance(old_item, (int, float)):
                delta = new_item - old_item
                v["score_lower_bound"] = round((v.get("score_lower_bound") or 0) + delta, 1)
                v["score_upper_bound"] = round((v.get("score_upper_bound") or 0) + delta, 1)
            else:
                v["score_lower_bound"] = round((v.get("score_lower_bound") or 0) + new_item, 1)
                v["unknown_items"] = max(0, (v.get("unknown_items") or 1) - 1)
    v["fatal_findings"] = [f for f in (v.get("fatal_findings") or []) if "recompute not possible" not in f]
    if d is not None and d > 0.005:
        v["fatal_findings"].append(f"fixed-model R1 (from CIF) differs from reported by {d:.4f}")
    rep["verdict"] = v


def main() -> int:
    args = sys.argv[1:]
    rid = cif = hkl = out = None
    shel = None
    if "--shel" in args:
        i = args.index("--shel"); shel = float(args[i + 1]); del args[i:i + 2]
    for flag, name in (("--run-id", "rid"), ("--cif", "cif"), ("--hkl", "hkl"), ("--out", "out")):
        if flag in args:
            i = args.index(flag); locals_val = args[i + 1]; del args[i:i + 2]
            if name == "rid": rid = locals_val
            elif name == "cif": cif = Path(locals_val)
            elif name == "hkl": hkl = Path(locals_val)
            else: out = Path(locals_val)
    if rid:
        rep_p = CTRL / "verifier" / rid / "report.json"
        rep = json.loads(rep_p.read_text(encoding="utf-8"))
        arena = Path(json.loads((CTRL / "runs" / rid / "manifest.json").read_text(encoding="utf-8"))["arena"])
        picked = rep["delivery"]["picked"]
        cif = cif or arena / picked[".cif"]
        hkl = hkl or arena / picked[".hkl"]
        out = out or CTRL / "verifier" / rid / "recompute_from_cif"
    if not (cif and hkl and out):
        print(__doc__); return 2
    out.mkdir(parents=True, exist_ok=True)
    model = read_cif(cif)
    result = {"cif": str(cif), "hkl": str(hkl), "space_group": model["hm"], "wavelength": model["wavelength"],
              "n_sites": len(model["sites"]), "n_aniso": sum(1 for s in model["sites"] if s["aniso"]), "elements": sorted({s["element"] for s in model["sites"]}),
              "resolution_cutoff_A": shel, "cif_reported": {}, "variants": {}, "note": "delivered CIF model + delivered reflections evaluated by SHELXL with every atomic parameter fixed and only the "
                                                        "overall scale refined (independent program; model unchanged); dispersion terms for the wavelength from cctbx tables; "
                                                        "MERG 2 = Friedel pairs kept (non-centrosymmetric), MERG 4 = Friedel pairs merged and f'' = 0"}
    import gemmi
    block = gemmi.cif.read(str(cif)).sole_block()
    for tag in ("_refine_ls_R_factor_gt", "_refine_ls_wR_factor_ref", "_refine_ls_number_reflns", "_reflns_number_gt", "_refine_ls_goodness_of_fit_ref"):
        v = block.find_value(tag)
        if v:
            result["cif_reported"][tag] = v
    hk = write_hkl(hkl, out / "model.hkl")
    result["hkl_written"] = hk
    for merg in (2, 4):
        stem = f"m{merg}"
        (out / f"{stem}.ins").write_text(build_ins(model, stem, merg, shel), encoding="ascii", errors="replace")
        shutil.copy2(out / "model.hkl", out / f"{stem}.hkl")
        result["variants"][f"MERG{merg}"] = run_shelxl(out, stem)
    try:
        rep_r1 = float(str(result["cif_reported"].get("_refine_ls_R_factor_gt", "nan")).split("(")[0])
    except ValueError:
        rep_r1 = float("nan")
    best = None
    for k, v in result["variants"].items():
        if "R1_gt" in v:
            d = abs(v["R1_gt"] - rep_r1)
            if best is None or d < best[1]:
                best = (k, d, v["R1_gt"])
    result["closest_variant"] = {"variant": best[0], "R1_gt": best[2], "abs_difference_to_reported": round(best[1], 4)} if best else None
    (out.parent / "recompute_from_cif.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    if rid:
        rep["recompute_from_cif"] = result
        update_verdict(rep, result)
        rep_p.write_text(json.dumps(rep, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print(json.dumps({"space_group": result["space_group"], "sites": result["n_sites"], "reported": result["cif_reported"].get("_refine_ls_R_factor_gt"),
                      "variants": {k: {kk: v.get(kk) for kk in ("R1_gt", "n_gt", "R1_all", "n_all", "wR2", "exit_code", "warnings")} for k, v in result["variants"].items()},
                      "closest": result["closest_variant"]}, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
