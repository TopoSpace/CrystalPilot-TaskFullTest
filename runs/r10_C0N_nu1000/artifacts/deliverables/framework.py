"""Independent diffraction calculations using only Python, NumPy and SciPy.

No crystallographic software or library is used. Fractional coordinates are
column vectors; the cell matrix maps them to Cartesian Angstrom coordinates.
"""
from pathlib import Path
import ast
import hashlib
import json
import time

import numpy as np
from scipy import ndimage, optimize

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / "inputs"
CELL = np.array([[39.19, -39.19 / 2, 0],
                 [0, 39.19 * np.sqrt(3) / 2, 0], [0, 0, 16.61]])
INV = np.linalg.inv(CELL)
VOLUME = np.linalg.det(CELL)


def parse_symmetry():
    text = (INPUT / "start.ins").read_text()
    expressions = ["x,y,z"]
    for line in text.splitlines():
        if line.startswith("SYMM "):
            expressions.append(line[5:].strip())
    # This input has only linear expressions with integral coefficients.
    def linear(node):
        if isinstance(node, ast.Name) and node.id in ("x", "y", "z"):
            return np.eye(3, dtype=int)[("x", "y", "z").index(node.id)]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -linear(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return linear(node.left) + linear(node.right)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub):
            return linear(node.left) - linear(node.right)
        raise ValueError("Unexpected symmetry expression")
    ops = [np.array([linear(ast.parse(c, mode="eval").body)
                     for c in expr.split(",")]) for expr in expressions]
    # Positive LATT in this supplied file specifies inversion.
    ops += [-m for m in ops.copy()]
    ops = np.unique(np.array(ops), axis=0)
    assert len(ops) == 24
    for a in ops:
        assert np.allclose(a.T @ CELL.T @ CELL @ a, CELL.T @ CELL)
        for b in ops:
            assert any(np.array_equal(a @ b, c) for c in ops)
    return ops


OPS = parse_symmetry()


def equiv_hkl(hkl):
    return np.einsum("ni,kij->nkj", hkl, OPS)


def canonical(hkl):
    eq = equiv_hkl(hkl)
    # Unique numerical encoding sufficient for all input indices.
    m = int(np.max(np.abs(eq))) + 1
    code = (eq[:, :, 0] + m) * (2*m+1)**2
    code += (eq[:, :, 1] + m) * (2*m+1) + eq[:, :, 2] + m
    return eq[np.arange(len(eq)), np.argmax(code, axis=1)]


def reciprocal_s2(hkl):
    return np.sum((hkl @ INV)**2, axis=1) / 4


def orbit(x, tol=1e-5):
    xyz = (OPS @ np.asarray(x)) % 1
    xyz[np.isclose(xyz, 1, atol=1e-10)] = 0
    _, index = np.unique(np.round(xyz / tol).astype(int), axis=0,
                         return_index=True)
    return xyz[np.sort(index)]


# Conventional neutral-atom four-Gaussian approximations. The constants are
# explicit calculation inputs, not loaded from external programs or databases.
SCATTER = {
    "C": ([2.31, 1.02, 1.5886, 0.865], [20.8439, 10.2075, 0.5687, 51.6512], .2156),
    "O": ([3.0485, 2.2868, 1.5463, .867], [13.2771, 5.7011, .3239, 32.9089], .2508),
    "Zr": ([17.8765, 10.948, 5.41732, 3.65721],
           [1.27618, 11.916, .117622, 87.6627], 2.06929),
    "H": ([.489918, .262003, .196767, .049879],
          [20.6593, 7.74039, 49.5519, 2.20159], .001305),
}


def form_factor(element, s2):
    a, b, c = SCATTER[element]
    return np.exp(-np.outer(s2, b)) @ np.asarray(a) + c


def structure_factor(hkl, atoms):
    s2 = reciprocal_s2(hkl)
    fc = np.zeros(len(hkl))
    for at in atoms:
        if "u_cart" in at:
            xyz = OPS @ np.asarray(at["xyz"])
            mult = len(orbit(at["xyz"]))
            rot = np.einsum("ij,kjl,lm->kim", CELL, OPS, INV)
            g = np.einsum("ni,kij->nkj", hkl @ INV, rot)
            dw = np.exp(-2*np.pi**2 * np.einsum(
                "nki,ij,nkj->nk", g, np.array(at["u_cart"]), g))
            phase = (np.cos(2*np.pi*(hkl@xyz.T))*dw).sum(axis=1) * mult/24
        else:
            xyz = orbit(at["xyz"])
            phase = np.cos(2 * np.pi * (hkl @ xyz.T)).sum(axis=1)
            phase *= np.exp(-8 * np.pi**2 * at.get("u", .04) * s2)
        fc += at.get("occ", 1) * (form_factor(at["element"], s2)
                                 + at.get("fp", 0)) * phase
    return fc


def load_data():
    return dict(np.load(ROOT / "merged.npz"))


def merge():
    begin = time.time()
    rows = []
    invalid = 0
    for line in (INPUT / "start.hkl").read_text().splitlines():
        try:
            h = [int(line[i:i+4]) for i in (0, 4, 8)]
            vals = [float(line[i:i+8]) for i in (12, 20)]
        except ValueError:
            invalid += 1
            continue
        if h == [0, 0, 0]:
            continue
        rows.append(h + vals)
    raw = np.asarray(rows)
    hkl = raw[:, :3].astype(int)
    intensity, sigma = raw[:, 3], raw[:, 4]
    keep = np.isfinite(intensity) & np.isfinite(sigma) & (sigma > 0)
    hkl, intensity, sigma = hkl[keep], intensity[keep], sigma[keep]
    can = canonical(hkl)
    unique, inv, counts = np.unique(can, axis=0, return_inverse=True,
                                    return_counts=True)
    w = 1 / sigma**2
    sw = np.bincount(inv, weights=w)
    mean = np.bincount(inv, weights=w * intensity) / sw
    simple = np.bincount(inv, weights=intensity) / counts
    residual = intensity - mean[inv]
    chi2 = np.bincount(inv, weights=w * residual**2)
    internal_sigma = 1 / np.sqrt(sw)
    # Do not reduce uncertainties below the counting-statistics estimate.
    external_sigma = internal_sigma * np.sqrt(
        np.maximum(1, chi2 / np.maximum(counts - 1, 1)))
    s2 = reciprocal_s2(unique)
    d = 1 / np.sqrt(4 * s2)
    np.savez_compressed(ROOT / "merged.npz", hkl=unique, i=mean,
                        sig=external_sigma, sig_internal=internal_sigma,
                        multiplicity=counts, s2=s2, d=d)
    np.savetxt(ROOT / "merged.hkl",
               np.column_stack((unique, mean, external_sigma)),
               fmt="%4d%4d%4d%14.6f%14.6f",
               header="h k l I sigma; inverse-variance merge in supplied 6/mmm")
    shells = []
    for lo, hi in [(0, .8), (.8, 1), (1, 1.2), (1.2, 1.5),
                   (1.5, 2), (2, 3), (3, 5), (5, 100)]:
        select = (d >= lo) & (d < hi)
        obs = select[inv]
        if not select.any():
            continue
        shells.append(dict(d_min=lo, d_max=hi, unique=int(select.sum()),
                           observations=int(obs.sum()),
                           mean_i_sigma=float(np.mean(mean[select] /
                                                      external_sigma[select])),
                           rint=float(np.sum(np.abs(intensity[obs] -
                                       simple[inv[obs]])) /
                                      np.sum(np.abs(intensity[obs])))))
    extent = np.max(np.abs(hkl), axis=0)
    grid = np.array(np.meshgrid(*(np.arange(-x, x+1) for x in extent),
                               indexing="ij")).reshape(3, -1).T
    gd2 = 4 * reciprocal_s2(grid)
    grid = grid[(gd2 > 0) & (gd2 <= 1 / d.min()**2 + 1e-10)]
    expected = np.unique(canonical(grid), axis=0)
    completeness = len(unique) / len(expected)
    stats = dict(observations=len(hkl), unique=len(unique), invalid_lines=invalid,
                 hkl_extent=extent.tolist(), d_min=float(d.min()),
                 d_max=float(d.max()), multiplicity=float(counts.mean()),
                 completeness_to_measured_dmin=float(completeness),
                 possible_unique_within_measured_index_box=len(expected),
                 rint=float(np.sum(np.abs(intensity - simple[inv])) /
                            np.sum(np.abs(intensity))),
                 r_weighted_equivalent_chi2=float(chi2.sum() /
                                                  (len(hkl)-len(unique))),
                 positive=int((mean > 0).sum()),
                 observed_gt_2sig=int((mean > 2 * external_sigma).sum()),
                 shells=shells, seconds=time.time()-begin,
                 input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in INPUT.iterdir() if p.is_file()})
    (ROOT / "data_statistics.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), flush=True)


def expand_coefficients(hkl, coefficients):
    eq = equiv_hkl(hkl).reshape(-1, 3)
    c = np.repeat(coefficients, len(OPS))
    unique, index = np.unique(eq, axis=0, return_index=True)
    return unique, c[index]


def fourier_map(hkl, coefficients, shape=(192, 192, 96)):
    indices, coef = expand_coefficients(hkl, coefficients)
    array = np.zeros(shape, dtype=complex)
    ix = tuple((indices % np.array(shape)).T)
    array[ix] = coef
    return np.fft.fftn(array).real / VOLUME


def unique_peaks(density, n=80, min_distance=.6):
    maxima = ndimage.maximum_filter(density, size=3, mode="wrap")
    inds = np.argwhere((density == maxima) & (density > 0))
    vals = density[tuple(inds.T)]
    order = np.argsort(vals)[::-1]
    accepted, expanded = [], []
    for ix in order:
        pos = inds[ix] / np.array(density.shape)
        if expanded:
            diff = pos - np.vstack(expanded)
            diff -= np.round(diff)
            dist = np.linalg.norm(diff @ CELL.T, axis=1)
            if dist.min() < min_distance:
                continue
        accepted.append((float(vals[ix]), pos))
        expanded.append(orbit(pos, tol=1e-4))
        if len(accepted) >= n:
            break
    return accepted


def patterson():
    data = load_data()
    # Positive observed intensities are used without artificial sharpening.
    density = fourier_map(data["hkl"], np.maximum(data["i"], 0))
    np.save(ROOT / "patterson.npy", density)
    peaks = unique_peaks(density, n=100)
    out = [{"height": v, "xyz": p.tolist()} for v, p in peaks]
    (ROOT / "patterson_peaks.json").write_text(json.dumps(out, indent=2))
    for i, (v, p) in enumerate(peaks):
        print(f"{i:3d} {v:12.5f} {p} distance {np.linalg.norm(CELL@p):.3f}")


if __name__ == "__main__":
    import sys
    {"merge": merge, "patterson": patterson}[sys.argv[1]]()
