"""Heavy-atom hypotheses derived from the measured Patterson vectors."""
import json
import time
import numpy as np
from scipy.optimize import least_squares
from framework import (ROOT, OPS, load_data, form_factor, reciprocal_s2,
                       structure_factor, fourier_map, unique_peaks, orbit)


def atoms_from_parameters(p, kind):
    x, q, z, u1, u2 = p[:5]
    if kind == 0:
        xyz1, xyz2 = [x, 0, z], [.5 + q, 2*q, 0]
    else:
        xyz1, xyz2 = [.5+q, 2*q, z], [x, 0, 0]
    return [dict(label="Zr1", element="Zr", xyz=xyz1, u=u1),
            dict(label="Zr2", element="Zr", xyz=xyz2, u=u2)]


def calculate(p, h, s2, f, kind):
    ats = atoms_from_parameters(p, kind)
    ans = np.zeros(len(h))
    for at, mult in zip(ats, [12, 6]):
        xyz = OPS @ np.array(at["xyz"])
        ans += (np.cos(2*np.pi * (h @ xyz.T)).sum(axis=1) * mult/24
                * f * np.exp(-8*np.pi**2 * at["u"] * s2))
    return ans


def main():
    t = time.time()
    d = load_data()
    use = (d["d"] > 1.25) & (d["d"] < 6)
    h, i, sig, s2 = (d[k][use] for k in ("hkl", "i", "sig", "s2"))
    f = form_factor("Zr", s2)
    weights = 1 / np.sqrt(sig**2 + (.10 * np.maximum(i, 0))**2 + .2**2)
    results = []
    for kind in (0, 1):
        for z in (.109, .19):
            p0 = [.454, .036, z, .05, .05]
            fc = calculate(p0, h, s2, f, kind)
            k = np.sum(weights**2*i*fc**2) / np.sum(weights**2*fc**4)
            p0 += [np.log(max(k, 1e-7))]
            def fun(p):
                fc = calculate(p, h, s2, f, kind)
                return (np.exp(p[5])*fc**2-i)*weights
            fit = least_squares(fun, p0,
                                bounds=([.40, .01, .03, .002, .002, -20],
                                        [.49, .08, .24, .4, .4, 5]),
                                max_nfev=150, ftol=1e-10, xtol=1e-10)
            atoms = atoms_from_parameters(fit.x, kind)
            r = dict(kind=kind, parameters=fit.x.tolist(), cost=fit.cost,
                     atoms=atoms, success=bool(fit.success))
            results.append(r)
            print(json.dumps(r), flush=True)
    best = min(results, key=lambda r: r["cost"])
    (ROOT / "heavy_trials.json").write_text(json.dumps(results, indent=2))
    atoms = best["atoms"]
    (ROOT / "heavy_atoms.json").write_text(json.dumps(atoms, indent=2))
    fc = structure_factor(d["hkl"], atoms)
    k = np.exp(best["parameters"][5])
    keep = d["d"] > 1.25
    # Observed Fourier with heavy-atom phases and difference Fourier.
    fo = np.sqrt(np.maximum(d["i"], 0) / k)
    coef = fo * np.sign(fc)
    density = fourier_map(d["hkl"][keep], coef[keep])
    diff = fourier_map(d["hkl"][keep], (coef-fc)[keep])
    np.save(ROOT / "heavy_fourier.npy", density)
    np.save(ROOT / "heavy_difference.npy", diff)
    peaks = unique_peaks(density, n=100, min_distance=.7)
    rows = [dict(height=v, xyz=p.tolist()) for v, p in peaks]
    (ROOT / "heavy_fourier_peaks.json").write_text(json.dumps(rows, indent=2))
    for j, (v, p) in enumerate(peaks):
        print(f"peak {j:3d} {v:9.4f} {p}", flush=True)
    print("SECONDS", time.time()-t)


if __name__ == "__main__":
    main()
