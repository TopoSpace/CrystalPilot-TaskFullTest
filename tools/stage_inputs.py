"""Seal the two cases' inputs into staging/ with a per-file exclusion record and SHA-256 manifests.

Alanine: the raw Rigaku frames and the acquisition metadata are copied; every file that carries a
data-reduction result or a database lookup (space-group determination, CSD cell matches, reduction
plots and intermediates) is excluded and listed with its reason. The source tree is never modified.
NU-1000: the three files named by the user are copied byte-for-byte under the names the protocol fixes.
"""
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
SRC_ALA = Path(r"H:\CrystalPilot-campaigns\usertest\test-alanine\alanine\CCDC_2360269")
SRC_NU = Path(r"H:\CrystalPilot-campaigns\usertest\test-NU1000-4")
STG = ROOT / "staging"
PROV = ROOT / "control" / "provenance"

# (relative-path predicate, reason). First match wins.
EXCLUDE = [
    (lambda r: r.startswith("tmp/csdsearches_Rigaku/"),
     "CSD cell-check results: list the deposited structures matching the cell, with formula and space group"),
    (lambda r: r.startswith("log/crysalispro_redLOG"),
     "CrysAlisPro data-reduction log: contains the automatic space-group determination and merging statistics"),
    (lambda r: r.startswith("expinfo/") and "datared" in r,
     "data-reduction settings/report: record the selected space group and reduction parameters"),
    (lambda r: r.startswith("plots_red/"),
     "data-reduction plots (Rint, lattice, profile fitting): derived results of a completed reduction"),
    (lambda r: r.startswith("tmp/"),
     "data-reduction intermediates (background images, profile tables, incident/polarisation corrections, resume files)"),
    (lambda r: r.startswith("original/") and r.endswith(".tabbin"),
     "peak-hunt and profile-fit tables: outputs of processing, not acquisition metadata"),
    (lambda r: r == "expinfo/struct_alias_db.ini",
     "structure alias database entry created by the reduction workflow"),
    (lambda r: r.startswith("log/crysalispro_userLOG") or r.startswith("log/crysalispro_ccdLOG"),
     "CrysAlisPro operator and instrument session logs: contain the reduction commands issued at the instrument"),
]
KEEP_NOTE = {
    "frames/": "raw diffraction frames (.rodhypix) and the crystal snapshots taken during collection",
    ".par": "CrysAlisPro parameter file: instrument geometry, goniometer and detector calibration used for "
            "collection (also carries the indexing cell; kept because the vendor reader may need it)",
    ".run": "run list: scan definitions actually collected",
    "expinfo/": "experiment settings written at collection time (sample formula, crystal description, "
                "data-collection strategy, coverage)",
    "kaboom/": "instrument collision models (geometry only)",
    ".ccd|.modulegeo|.CAP_shape|.acc|.ini|.lock|.sel_od|.sum|.runbup|.bup":
        "instrument and detector calibration, crystal-shape and accessory definitions, collection summaries",
    "movie/": "crystal video recorded during centring",
    "bup/": "backups of the run list and crystal-shape files made during collection",
}


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(root: Path, rels: list) -> dict:
    out = {}
    for r in rels:
        p = root / r
        st = p.stat()
        out[r] = {"sha256": sha256_of(p), "size": st.st_size,
                  "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime))}
    return out


def stage_alanine() -> None:
    dst = STG / "alanine" / "inputs" / "alanine"
    if dst.exists():
        print("staging/alanine exists - refusing to overwrite")
        sys.exit(2)
    kept, excluded = [], []
    for folder, _, files in os.walk(SRC_ALA):
        for f in files:
            p = Path(folder) / f
            rel = p.relative_to(SRC_ALA).as_posix()
            reason = next((why for pred, why in EXCLUDE if pred(rel)), None)
            if reason:
                excluded.append({"path": rel, "size": p.stat().st_size, "reason": reason})
            else:
                kept.append(rel)
    for rel in kept:
        t = dst / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SRC_ALA / rel, t)
    src_manifest = manifest(SRC_ALA, kept)
    dst_manifest = manifest(dst, kept)
    mismatch = [r for r in kept if src_manifest[r]["sha256"] != dst_manifest[r]["sha256"]]
    rec = {
        "source": str(SRC_ALA), "staged_to": str(dst),
        "staged_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "kept_files": len(kept), "kept_bytes": sum(v["size"] for v in dst_manifest.values()),
        "excluded_files": len(excluded), "excluded_bytes": sum(e["size"] for e in excluded),
        "copy_hash_mismatches": mismatch,
        "source_folder_name_note": "the source folder is named after a CCDC deposition number; the staged "
                                   "copy is named alanine so the number is not exposed to participants",
        "keep_rationale": KEEP_NOTE, "excluded": excluded,
    }
    (PROV / "alanine_staging.json").write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")
    (PROV / "alanine_source_manifest.json").write_text(json.dumps(src_manifest, ensure_ascii=False, indent=1),
                                                       encoding="utf-8")
    (STG / "alanine" / "manifest.json").write_text(json.dumps(dst_manifest, ensure_ascii=False, indent=1),
                                                   encoding="utf-8")
    print(f"alanine: kept {len(kept)} files / {rec['kept_bytes']/1e6:.1f} MB, "
          f"excluded {len(excluded)} files / {rec['excluded_bytes']/1e6:.1f} MB, mismatches {len(mismatch)}")
    by_reason: dict = {}
    for e in excluded:
        by_reason[e["reason"]] = by_reason.get(e["reason"], 0) + 1
    for k, v in by_reason.items():
        print(f"   excluded {v:3d}  {k[:95]}")


def stage_nu() -> None:
    dst = STG / "nu1000" / "inputs"
    if dst.exists():
        print("staging/nu1000 exists - refusing to overwrite")
        sys.exit(2)
    dst.mkdir(parents=True)
    pairs = [("zr-NU-1000-with-guest.hkl", "start.hkl"), ("zr-NU-1000-with-guest.ins", "start.ins"),
             ("uploads/图片.png", "synthesis.png")]
    rec = {"source": str(SRC_NU), "staged_to": str(dst),
           "staged_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "files": []}
    for s, d in pairs:
        shutil.copy2(SRC_NU / s, dst / d)
        a, b = sha256_of(SRC_NU / s), sha256_of(dst / d)
        rec["files"].append({"source_rel": s, "staged_as": d, "sha256": a, "identical": a == b,
                             "size": (dst / d).stat().st_size})
    ins = (dst / "start.ins").read_text(encoding="utf-8", errors="replace").splitlines()
    cards: dict = {}
    for line in ins:
        k = line.split()[0].upper() if line.split() else ""
        cards[k] = cards.get(k, 0) + 1

    def first(prefix: str):
        return next((l for l in ins if l.upper().startswith(prefix)), None)

    rec["start_ins_summary"] = {
        "title": ins[0] if ins else None, "cell_line": first("CELL"), "sfac_line": first("SFAC"),
        "latt_line": first("LATT"), "n_symm": cards.get("SYMM", 0), "hklf_line": first("HKLF"),
        "cards_present": sorted(k for k in cards if k),
        "has_abin_or_fab": any(k in cards for k in ("ABIN", "FAB")),
        "has_twin_basf": any(k in cards for k in ("TWIN", "BASF")),
        "has_omit_shel": any(k in cards for k in ("OMIT", "SHEL")),
        "atom_lines": [l for l in ins if l.split() and l.split()[0] in ("O", "N", "C", "H") and len(l.split()) >= 5],
        "note": "placeholder atoms at the origin only; SFAC lists no metal; the model is effectively empty, "
                "so the task starts from cell, symmetry and reflections",
    }
    hkl = dst / "start.hkl"
    with hkl.open("rb") as fh:
        n = sum(1 for _ in fh)
    rec["start_hkl_summary"] = {"lines": n, "bytes": hkl.stat().st_size}
    (PROV / "nu1000_staging.json").write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")
    (STG / "nu1000" / "manifest.json").write_text(
        json.dumps({f["staged_as"]: {"sha256": f["sha256"], "size": f["size"]} for f in rec["files"]}, indent=1),
        encoding="utf-8")
    print("nu1000:", [(f["staged_as"], f["identical"]) for f in rec["files"]], "| hkl lines", n)


if __name__ == "__main__":
    PROV.mkdir(parents=True, exist_ok=True)
    stage_alanine()
    stage_nu()
