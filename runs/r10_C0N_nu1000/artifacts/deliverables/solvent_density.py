"""Experimental pore-only density correction, with independent holdout checks.

This is an original alternating-projection algorithm, not a molecular solvent
model. The supplied reflection intensities are never overwritten.
"""
import argparse
import json
import time
import numpy as np
from scipy.spatial import cKDTree
from framework import (ROOT, CELL, VOLUME, orbit, load_data, structure_factor,
                       expand_coefficients, reciprocal_s2)
from refine import metrics


def make_mask(atoms, shape, radius):
    xyz = np.vstack([orbit(a["xyz"]) for a in atoms if a["element"] != "H"])
    shifts = np.array(np.meshgrid(*([[-1, 0, 1]]*3), indexing="ij")).reshape(3, -1).T
    expanded = (xyz[None] + shifts[:, None]).reshape(-1, 3) @ CELL.T
    tree = cKDTree(expanded)
    ix = np.indices(shape, dtype=float).reshape(3, -1).T / np.array(shape)
    dist = tree.query(ix @ CELL.T, workers=1)[0].reshape(shape)
    mask = np.clip((dist-radius)/.25+.5, 0, 1)
    return mask


def free_flags(data, fraction=.10):
    # Reproducible stratified selection in resolution order, independent of I.
    rng = np.random.default_rng(20260922)
    order = np.argsort(data["s2"])
    flags = np.zeros(len(order), dtype=bool)
    for block in np.array_split(order, max(1, len(order)//100)):
        flags[rng.choice(block, size=max(1, round(len(block)*fraction)),
                         replace=False)] = True
    return flags


def run(args):
    t = time.time()
    atoms = json.loads((ROOT / f"{args.input}_atoms.json").read_text())
    stats = json.loads((ROOT / f"{args.input}_log.json").read_text())
    data = load_data()
    scale = stats["all_data"]["scale"]
    fc = structure_factor(data["hkl"], atoms)
    flags = free_flags(data)
    use = (data["d"] >= args.dmin) & (~flags if args.holdout else True)
    h = data["hkl"][use]
    fo = np.sqrt(np.maximum(data["i"][use], 0)/scale)
    grid_shape = (128, 128, 64)
    mask = make_mask(atoms, grid_shape, args.radius)
    rho = mask * args.initial_density
    indices, _ = expand_coefficients(h, np.zeros(len(h)))
    # Build the exact mapping from unique observed coefficients to all symmetry
    # equivalents once, preserving the multiplicities of special reflections.
    from framework import equiv_hkl
    eq = equiv_hkl(h).reshape(-1, 3)
    _, ids = np.unique(eq, axis=0, return_index=True)
    mapping = np.repeat(np.arange(len(h)), 24)[ids]
    ix = tuple((indices % np.array(grid_shape)).T)
    all_ix = tuple((data["hkl"] % np.array(grid_shape)).T)
    # Isotropic Cartesian smoothing kernel, including the hexagonal metric.
    mesh = np.meshgrid(*(np.fft.fftfreq(n)*n for n in grid_shape), indexing="ij")
    grid_h = np.stack(mesh, axis=-1).reshape(-1, 3)
    kernel = np.exp(-8*np.pi**2*.55**2*reciprocal_s2(grid_h)).reshape(grid_shape)
    history = []
    for iteration in range(args.iterations):
        fs_grid = np.fft.ifftn(rho).real * VOLUME
        fs = fs_grid[all_ix]
        total = fc[use]+fs[use]
        delta = fo*np.where(total >= 0, 1, -1)-total
        coeff = np.zeros(grid_shape)
        coeff[ix] = delta[mapping]
        gradient = np.fft.fftn(coeff).real / VOLUME
        if args.smooth > 0:
            smoothed = np.fft.fftn(np.fft.ifftn(rho)*kernel).real
            gradient -= args.smooth*(rho-smoothed)
        rho = mask*np.clip(rho+args.step*gradient, 0, args.cap)
        if iteration % 20 == 0 or iteration == args.iterations-1:
            log = dict(iteration=iteration, electrons=float(rho.mean()*VOLUME),
                       maximum_density=float(rho.max()))
            for label, selection in [
                ("work", (data["d"] >= args.dmin) & ~flags),
                ("free", (data["d"] >= args.dmin) & flags),
                ("all", data["d"] >= 1.1)]:
                dd = {k:v[selection] for k,v in data.items()}
                log[label] = metrics(dd, (fc+fs)[selection], scale)
            history.append(log)
            print(json.dumps(log), flush=True)
    fs_grid = np.fft.ifftn(rho).real * VOLUME
    fs = fs_grid[all_ix]
    # Retain the Fourier transform of the density rather than artificially
    # truncating its high-angle contribution.
    np.savez_compressed(ROOT / f"{args.name}.npz", hkl=data["hkl"], fs=fs,
                        free=flags, mask=mask.astype(np.float32),
                        rho=rho.astype(np.float32), scale=scale)
    (ROOT / f"{args.name}_log.json").write_text(json.dumps(dict(
        arguments=vars(args), history=history, volume=float(VOLUME),
        pore_fraction=float(mask.mean()), seconds=time.time()-t), indent=2))
    print("SECONDS", time.time()-t, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="highangle")
    parser.add_argument("--name", default="pore_trial")
    parser.add_argument("--dmin", type=float, default=2.5)
    parser.add_argument("--radius", type=float, default=1.65)
    parser.add_argument("--iterations", type=int, default=201)
    parser.add_argument("--step", type=float, default=.8)
    parser.add_argument("--smooth", type=float, default=.05)
    parser.add_argument("--cap", type=float, default=1.2)
    parser.add_argument("--initial-density", type=float, default=.1)
    parser.add_argument("--holdout", action="store_true")
    run(parser.parse_args())
