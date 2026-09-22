"""Structural sanity checks of every delivered primary CIF (bond lengths by element pair, hydrogen geometry,
displacement parameters), as a fact base for the expert review section of the report. Read-only; writes
control/analysis/structure_checks.json and prints a compact table. Literature reference for L-alanine at 23 K
(Destro, Marsh, Bianchi, J. Phys. Chem. 1988, 92, 966): N-Calpha 1.4875, Calpha-C 1.5350, Calpha-Cbeta 1.5240,
C-O 1.2484 / 1.2647 A."""
import json
import math
import sys
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
RUNS = ["r01_H1_nu1000", "r04_C0_nu1000", "r05_H0_nu1000", "r07_P1_nu1000", "r09_P1N_nu1000", "r10_C0N_nu1000",
        "r02_H1_alanine", "r03_C0_alanine", "r06_H0_alanine", "r08_P1_alanine", "r11_C0N_alanine", "r12_P1N_alanine"]
LIT_ALA = {"N-C": 1.4875, "C-C(carboxyl)": 1.5350, "C-C(methyl)": 1.5240, "C-O": (1.2484, 1.2647)}
COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "Zr": 1.75}


def analyse(cif: Path) -> dict:
    import gemmi
    st = gemmi.read_small_structure(str(cif))
    cell = st.cell
    sg = gemmi.find_spacegroup_by_name(st.spacegroup_hm) if st.spacegroup_hm else None
    sites = st.sites
    # expand all symmetry images inside and around the cell, then search neighbours of each independent site
    ops = sg.operations() if sg else gemmi.GroupOps([gemmi.Op("x,y,z")])
    images = []
    for s in sites:
        for op in ops:
            fx, fy, fz = op.apply_to_xyz([s.fract.x, s.fract.y, s.fract.z])
            images.append((s.label, s.element.name, fx % 1.0, fy % 1.0, fz % 1.0))
    bonds = []
    for s in sites:
        el = s.element.name
        p0 = cell.orthogonalize(gemmi.Fractional(s.fract.x, s.fract.y, s.fract.z))
        for lab2, el2, fx, fy, fz in images:
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        p1 = cell.orthogonalize(gemmi.Fractional(fx + dx, fy + dy, fz + dz))
                        d = p0.dist(p1)
                        cutoff = (COV.get(el, 1.0) + COV.get(el2, 1.0)) * 1.25
                        if 0.5 < d < cutoff:
                            key = tuple(sorted([el, el2]))
                            bonds.append((s.label, lab2, "-".join(key), round(d, 4)))
    by_pair = {}
    for a, b, k, d in bonds:
        by_pair.setdefault(k, []).append(d)
    summary = {k: {"n": len(v), "min": min(v), "max": max(v), "mean": round(sum(v) / len(v), 4)} for k, v in by_pair.items()}
    n_h = sum(1 for s in sites if s.element.name == "H")
    ueq = []
    npd = []
    aniso_ratio = []
    for s in sites:
        if s.aniso.nonzero():
            a = s.aniso
            ev = sorted(gemmi.SMat33d(a.u11, a.u22, a.u33, a.u12, a.u13, a.u23).calculate_eigenvalues()) if hasattr(gemmi, "SMat33d") else None
            if ev:
                if ev[0] <= 0:
                    npd.append(s.label)
                if ev[0] > 0:
                    aniso_ratio.append((s.label, round(ev[2] / ev[0], 1)))
            ueq.append((s.label, round((a.u11 + a.u22 + a.u33) / 3, 4)))
        elif s.u_iso:
            ueq.append((s.label, round(s.u_iso, 4)))
    non_h_ueq = [u for l, u in ueq if not l.upper().startswith("H")]
    per_atom = {}
    for a, b, k, d in bonds:
        per_atom.setdefault(a, []).append((b, d))
    return {"cell": [round(cell.a, 4), round(cell.b, 4), round(cell.c, 4), round(cell.alpha, 2), round(cell.beta, 2), round(cell.gamma, 2)],
            "space_group": st.spacegroup_hm, "n_sites": len(sites), "n_H": n_h, "bond_summary": summary,
            "ueq_nonH_min": min(non_h_ueq) if non_h_ueq else None, "ueq_nonH_max": max(non_h_ueq) if non_h_ueq else None,
            "npd": npd, "max_aniso_ratio": max(aniso_ratio, key=lambda t: t[1]) if aniso_ratio else None,
            "per_atom_bonds": {a: sorted(v, key=lambda t: t[1]) for a, v in per_atom.items()}}


def main() -> int:
    out = {}
    for rid in RUNS:
        rep = json.loads((CTRL / "verifier" / rid / "report.json").read_text(encoding="utf-8"))
        arena = Path(json.loads((CTRL / "runs" / rid / "manifest.json").read_text(encoding="utf-8"))["arena"])
        cif = arena / rep["delivery"]["picked"][".cif"]
        try:
            res = analyse(cif)
        except Exception as e:  # noqa: BLE001
            res = {"error": f"{type(e).__name__}: {e}"}
        res["cif"] = str(cif)
        out[rid] = res
        bs = res.get("bond_summary", {})
        print(f"== {rid} {res.get('space_group')} sites={res.get('n_sites')} H={res.get('n_H')} Ueq(nonH) {res.get('ueq_nonH_min')}..{res.get('ueq_nonH_max')} NPD={res.get('npd')} max_aniso={res.get('max_aniso_ratio')}")
        for k in sorted(bs):
            v = bs[k]
            print(f"     {k:6s} n={v['n']:3d} {v['min']:.3f}..{v['max']:.3f} mean {v['mean']:.3f}")
    (CTRL / "analysis" / "structure_checks.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
