"""Map-peak diagnostics and local chemical-neighbour distances."""
import itertools
import json
import numpy as np
from scipy.optimize import minimize
from framework import (ROOT, CELL, OPS, orbit, load_data, structure_factor,
                       fourier_map, unique_peaks)


def representative(x):
    for pos in orbit(x):
        for shift in itertools.product((-1, 0, 1), repeat=2):
            p = pos.copy()
            p[:2] += shift
            if (p[1] >= -1e-6 and p[0] >= 2*p[1]-1e-6
                    and p[0]-.5*p[1] <= .500001 and p[2] <= .500001):
                return p
    return np.array(x)


def snap_site(x, distance=.22):
    """Project onto a nearby fixed symmetry subspace in Cartesian metric."""
    x = np.array(x)
    candidates = []
    for op in OPS:
        delta = op @ x - x
        delta -= np.round(delta)
        if np.linalg.norm(CELL @ delta) < 2*distance:
            candidates.append(op)
    # Simultaneous projection onto detected affine symmetry constraints.
    if candidates:
        a = np.vstack([op-np.eye(3) for op in candidates])
        b = np.concatenate([np.round((op-np.eye(3)) @ x) for op in candidates])
        ac = a @ np.linalg.inv(CELL)
        corr = np.linalg.lstsq(ac, b-a@x, rcond=None)[0]
        x += np.linalg.solve(CELL, corr)
    return representative(x % 1)


def distance(a, b):
    dx = np.asarray(a) - np.asarray(b)
    dx -= np.round(dx)
    return np.linalg.norm(dx @ CELL.T, axis=-1)


def main():
    data = load_data()
    atoms = json.loads((ROOT / "heavy_atoms.json").read_text())
    zr = np.vstack([orbit(a["xyz"]) for a in atoms])
    rows = json.loads((ROOT / "heavy_fourier_peaks.json").read_text())
    points = [snap_site(r["xyz"]) for r in rows]
    for j, p in enumerate(points[:45]):
        neighbours = []
        for k, q in enumerate(points[:45]):
            dd = distance(p, orbit(q))
            dd = dd[dd > .2]
            if len(dd) and dd.min() < 1.8:
                neighbours.append(f"{k}:{dd.min():.2f}")
        print(f'{j:2d} {rows[j]["height"]:7.2f} '
              f'{p[0]:.6f} {p[1]:.6f} {p[2]:.6f} '
              f'mult {len(orbit(p)):2d} Zr {distance(p, zr).min():.3f} '
              + " ".join(neighbours))
    fc = structure_factor(data["hkl"], atoms)
    fo = np.sqrt(np.maximum(data["i"], 0))
    use = data["i"] > 2*data["sig"]
    scale = np.sum(fo[use]*np.abs(fc[use])) / np.sum(fc[use]**2)
    print("Amplitude scale", scale, "R1",
          np.sum(np.abs(fo[use]-scale*np.abs(fc[use])))/np.sum(fo[use]))
    print("Top observed", sorted(zip(data["i"], data["hkl"].tolist()), reverse=True)[:20])


if __name__ == "__main__":
    main()
