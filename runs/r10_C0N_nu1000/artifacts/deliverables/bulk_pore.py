"""Low-dimensional pore-density sensitivity test with held-out intensities."""
import argparse
import json
import numpy as np
from scipy.optimize import least_squares
from framework import ROOT, CELL, VOLUME, load_data, structure_factor
from solvent_density import make_mask, free_flags
from refine import metrics


def basis_density(atoms, hkl, kind):
    shape = (96, 96, 48)
    mask = make_mask(atoms, shape, 1.65)
    coords = np.indices(shape).reshape(3, -1).T / np.array(shape)
    xy = coords[:, :2]
    shifts = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1],
                                 indexing="ij")).reshape(2, -1).T
    def radial(centres):
        r = np.full(len(xy), 100.)
        for c in centres:
            for shift in shifts:
                delta = xy-c-shift
                cart = delta @ CELL[:2, :2].T
                r = np.minimum(r, np.linalg.norm(cart, axis=1))
        return r.reshape(shape)
    big = radial([np.zeros(2)])
    small = radial([np.array([1/3, 2/3]), np.array([2/3, 1/3])])
    partition = 1 / (1+np.exp(np.clip((big-small-11.31)/.5, -50, 50)))
    if kind == "flat":
        radial_basis = [mask]
    elif kind == "two":
        radial_basis = [mask*partition, mask*(1-partition)]
    elif kind == "layered":
        radial_basis = [mask*partition, mask*(1-partition)]
    else:
        radial_basis = [mask*partition*np.exp(-.5*((big-r)/3.5)**2)
                        for r in (0, 7, 13)]
        radial_basis += [mask*(1-partition)*np.exp(-.5*((small-r)/2.5)**2)
                         for r in (0, 4)]
    z = coords[:, 2].reshape(shape)
    if kind in ("flat", "two"):
        z_basis = [np.ones(shape)]
    else:
        z_basis = []
        for centre in (0, .25, .5):
            layer = np.zeros(shape)
            for zz in np.unique([centre, 1-centre]):
                dz = z-zz
                dz -= np.round(dz)
                layer += np.exp(-.5*(dz*16.61/1.8)**2)
            z_basis.append(layer)
    basis = [r*z for r in radial_basis for z in z_basis]
    ix = tuple((hkl % np.array(shape)).T)
    transforms = np.array([(np.fft.ifftn(b).real*VOLUME)[ix] for b in basis]).T
    means = np.array([b.mean()*VOLUME for b in basis])
    return transforms, means


def run(args):
    atoms = json.loads((ROOT / f"{args.input}_atoms.json").read_text())
    log = json.loads((ROOT / f"{args.input}_log.json").read_text())
    scale = log["all_data"]["scale"]
    data = load_data()
    fc = structure_factor(data["hkl"], atoms)
    basis, electrons = basis_density(atoms, data["hkl"], args.kind)
    flags = free_flags(data)
    use = (data["d"] > 2.5) & ~flags
    weight = 1/np.sqrt(data["sig"][use]**2 +
                        (.1*np.maximum(data["i"][use], 0))**2)
    trials = []
    rng = np.random.default_rng(423)
    for trial in range(6):
        start = np.full(basis.shape[1], .1) if trial == 0 else rng.uniform(0, .4, basis.shape[1])
        def fun(p):
            total = fc[use] + basis[use]@p
            return np.r_[(scale*total**2-data["i"][use])*weight, p/.5]
        def jac(p):
            total = fc[use] + basis[use]@p
            return np.vstack([2*scale*total[:, None]*basis[use]*weight[:, None],
                              np.eye(len(p))/.5])
        fit = least_squares(fun, start, jac=jac, bounds=(0, .8), max_nfev=200)
        trials.append((fit.cost, fit.x))
    cost, p = min(trials, key=lambda v:v[0])
    fs = basis@p
    results = dict(kind=args.kind, input=args.input, parameters=p.tolist(),
                   electrons=float(electrons@p), cost=cost)
    for label, select in [("work_low", (data["d"] > 2.5) & ~flags),
                           ("free_low", (data["d"] > 2.5) & flags),
                           ("all", data["d"] > 1.1)]:
        results[label] = metrics({k:v[select] for k,v in data.items()},
                                 (fc+fs)[select], scale)
    print(json.dumps(results, indent=2))
    (ROOT / f"{args.name}_log.json").write_text(json.dumps(results, indent=2))
    np.savez_compressed(ROOT / f"{args.name}.npz", hkl=data["hkl"], fs=fs,
                        free=flags, scale=scale, basis=basis, density_coefficients=p)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="highangle")
    parser.add_argument("--kind", choices=["flat", "two", "layered", "radial"], default="two")
    parser.add_argument("--name", default="bulk_two")
    run(parser.parse_args())
