"""Half-dataset correlation, resolution completeness, and numerical self-tests."""
import json
import numpy as np
from framework import (ROOT, INPUT, CELL, OPS, VOLUME, canonical, reciprocal_s2,
                       load_data, structure_factor, fourier_map, orbit)
from refine import Model
from restraints import Restraints


def run():
    rows = []
    for line in (INPUT / "start.hkl").read_text().splitlines():
        h = [int(line[i:i+4]) for i in (0, 4, 8)]
        if h == [0, 0, 0]:
            continue
        rows.append(h+[float(line[12:20]), float(line[20:28])])
    a = np.array(rows)
    keep = a[:, 4] > 0
    a = a[keep]
    h = a[:, :3].astype(int)
    i, sig = a[:, 3], a[:, 4]
    can = canonical(h)
    unique, inv, counts = np.unique(can, axis=0, return_inverse=True,
                                    return_counts=True)
    rng = np.random.default_rng(20260922)
    # Randomly order observations within each group, then alternating allocation
    # guarantees both halves whenever at least two observations exist.
    order = np.lexsort((rng.random(len(inv)), inv))
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    rank = np.arange(len(inv))-np.repeat(starts, counts)
    half = np.empty(len(inv), dtype=int)
    half[order] = rank % 2
    half_means = []
    for k in (0, 1):
        use = half == k
        sw = np.bincount(inv[use], weights=1/sig[use]**2, minlength=len(unique))
        si = np.bincount(inv[use], weights=i[use]/sig[use]**2, minlength=len(unique))
        half_means.append(np.divide(si, sw, out=np.full(len(unique), np.nan),
                                    where=sw > 0))
    d = 1/np.sqrt(4*reciprocal_s2(unique))
    max_h = int(np.ceil(39.19/d.min()))
    max_l = int(np.ceil(16.61/d.min()))
    grid = np.array(np.meshgrid(np.arange(-max_h, max_h+1),
                               np.arange(-max_h, max_h+1),
                               np.arange(-max_l, max_l+1),
                               indexing="ij")).reshape(3, -1).T
    s2 = reciprocal_s2(grid)
    grid = grid[(s2 > 0) & (s2 <= 1/(4*d.min()**2)+1e-10)]
    expected = np.unique(canonical(grid), axis=0)
    expected_d = 1/np.sqrt(4*reciprocal_s2(expected))
    data = load_data()
    shells = []
    for lo, hi in [(.691, .8), (.8, .9), (.9, 1), (1, 1.1), (1.1, 1.2),
                   (1.2, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 15), (15, 100)]:
        sel = (d >= lo) & (d < hi)
        valid = sel & np.isfinite(half_means[0]) & np.isfinite(half_means[1])
        cc = (np.corrcoef(half_means[0][valid], half_means[1][valid])[0, 1]
              if valid.sum() > 3 else None)
        expected_n = int(((expected_d >= lo) & (expected_d < hi)).sum())
        shells.append(dict(d_min=lo, d_max=hi, measured=int(sel.sum()),
                           expected=expected_n,
                           completeness=float(sel.sum()/expected_n) if expected_n else None,
                           cc_half=None if cc is None else float(cc),
                           paired_halves=int(valid.sum())))
    observed_set = {tuple(v) for v in unique}
    missing_low = [v.tolist() for v, dd in zip(expected, expected_d)
                   if dd > 8 and tuple(v) not in observed_set]
    audit = dict(shells=shells, volume=float(VOLUME), missing_d_above_8=missing_low,
                 wavelength=0.68883, random_seed=20260922,
                 note="CC1/2 from random inverse-variance half-dataset merges.")
    (ROOT / "resolution_audit.json").write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))

    atoms = json.loads((ROOT / "riding_atoms.json").read_text())
    atoms = [a for a in atoms if a["element"] != "H"]
    select = np.arange(20, len(data["hkl"]), 31)
    model = Model(atoms, data["hkl"][select], data["s2"][select],
                  aniso="all", zr_fp=True)
    p = model.p0.copy()
    fc, jac = model.compute(p)
    direct = structure_factor(data["hkl"][select], model.output_atoms(p))
    sf_error = float(np.max(np.abs(fc-direct)))
    selected_parameters = np.linspace(0, len(p)-2, 40, dtype=int)
    errors = []
    for ix in selected_parameters:
        step = 1e-6
        pp, pm = p.copy(), p.copy()
        pp[ix] += step
        pm[ix] -= step
        numerical = (model.compute(pp)[0]-model.compute(pm)[0])/(2*step)
        denom = max(np.max(np.abs(numerical)), 1)
        errors.append(float(np.max(np.abs(numerical-jac[:, ix]))/denom))
    restraints = Restraints(model)
    rr, jj = restraints.calculate(p)
    gradient = jj.T@rr
    rest_errors = []
    for ix in selected_parameters:
        step = 1e-6
        pp, pm = p.copy(), p.copy()
        pp[ix] += step
        pm[ix] -= step
        rp = restraints.calculate(pp)[0]
        rm = restraints.calculate(pm)[0]
        numerical = (rp@rp-rm@rm)/(4*step)
        rest_errors.append(float(abs(numerical-gradient[ix]) /
                                 max(abs(numerical), 1)))
    # Independent orbit invariance check on computed structure factors.
    htest = data["hkl"][select][:50]
    reference = structure_factor(htest, atoms)
    sym_error = max(float(np.max(np.abs(structure_factor(htest@op, atoms)-reference)))
                    for op in OPS)
    test_h = np.array([[1, 0, 0], [2, -1, 1]])
    test_c = np.array([2., 3.])
    density = fourier_map(test_h, test_c, shape=(24, 24, 16))
    transformed = np.fft.ifftn(density).real*VOLUME
    fft_error = float(np.max(np.abs(transformed[tuple(test_h.T)]-test_c)))
    result = dict(symmetry_group_order=len(OPS), sf_implementation_error=sf_error,
                  max_relative_sf_jacobian_error=max(errors),
                  max_relative_restraint_objective_gradient_error=max(rest_errors),
                  symmetry_equivalence_error=sym_error, fourier_roundtrip_error=fft_error)
    result["passed"] = (sf_error < 1e-8 and max(errors) < 1e-4
                        and max(rest_errors) < 1e-3 and sym_error < 1e-8
                        and fft_error < 1e-9)
    (ROOT / "numerical_self_tests.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise RuntimeError("Numerical audit failed")


if __name__ == "__main__":
    run()
