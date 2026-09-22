"""Symmetry-constrained least squares against measured intensities."""
import argparse
import json
import time
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares
from framework import (ROOT, CELL, INV, OPS, orbit, load_data, form_factor,
                       structure_factor, fourier_map, unique_peaks)
from inspect_maps import snap_site, distance


def initial_atoms():
    heavy = json.loads((ROOT / "heavy_atoms.json").read_text())
    for at in heavy:
        at["u"] = .035
    peaks = json.loads((ROOT / "heavy_fourier_peaks.json").read_text())
    atoms = heavy.copy()
    for j, ix in enumerate([2, 3, 4, 5, 6, 7], start=1):
        atoms.append(dict(label=f"O{j}", element="O",
                          xyz=snap_site(peaks[ix]["xyz"]).tolist(), u=.06))
    for j, ix in enumerate([8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 21], start=1):
        atoms.append(dict(label=f"C{j}", element="C",
                          xyz=snap_site(peaks[ix]["xyz"]).tolist(), u=.08))
    pos = np.array(atoms[14]["xyz"])  # C7, para position of phenyl ring
    other = np.array(atoms[10]["xyz"])  # C3, opposite phenyl position
    vec = pos - other
    vec *= 1.48 / np.linalg.norm(CELL @ vec)
    atoms.append(dict(label="C12", element="C", xyz=(pos+vec).tolist(), u=.08))
    return atoms


class Model:
    def __init__(self, atoms, hkl, s2, aniso="none", zr_fp=False):
        self.atoms = atoms
        self.hkl = hkl
        self.s2 = s2
        self.sites = []
        self.fp_index = None
        p, lower, upper = [], [], []
        self.f = {e: form_factor(e, s2) for e in {a["element"] for a in atoms}}
        for at in atoms:
            x = np.array(at["xyz"])
            stabilizers = []
            for op in OPS:
                dx = op@x-x
                dx -= np.round(dx)
                if np.linalg.norm(CELL @ dx) < 1e-4:
                    stabilizers.append(op)
            const = np.vstack([(op-np.eye(3)) @ INV for op in stabilizers])
            q = INV @ null_space(const)
            # Singular values of exactly-zero matrices should leave 3 DOF.
            n = q.shape[1]
            anisotropic = aniso == "all" or (aniso == "heavy" and at["element"] == "Zr")
            if anisotropic:
                ub0 = []
                for a, b in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]:
                    u = np.zeros((3, 3))
                    u[a, b] = u[b, a] = 1 if a == b else 1/np.sqrt(2)
                    ub0.append(u)
                ub0 = np.array(ub0)
                rs = np.einsum("ij,kjl,lm->kim", CELL, stabilizers, INV)
                projected = np.array([np.mean(rs @ u @ rs.transpose(0, 2, 1),
                                              axis=0) for u in ub0])
                matrix = projected.reshape(6, 9).T
                left, sv, _ = np.linalg.svd(matrix, full_matrices=False)
                ub = left[:, sv > 1e-7].T.reshape(-1, 3, 3)
            else:
                ub = np.eye(3)[None]
            nu = len(ub)
            initial_u = np.array(at.get("u_cart", np.eye(3)*at.get("u", .05)))
            pu = np.linalg.lstsq(ub.reshape(nu, 9).T,
                                initial_u.ravel(), rcond=None)[0]
            start = len(p)
            p.extend([0.0]*n + pu.tolist())
            lower.extend([-0.7]*n + ([-.7]*nu if anisotropic else [.001]))
            upper.extend([.7]*n + [.7]*nu)
            rot_q = np.einsum("kij,jm->kim", OPS, q)
            dphase = 2*np.pi*np.einsum("ni,kij->nkj", hkl, rot_q)
            rotations = np.einsum("ij,kjl,lm->kim", CELL, OPS, INV)
            g = np.einsum("ni,kij->nkj", hkl @ INV, rotations)
            d_exponent = -2*np.pi**2*np.einsum("nki,mij,nkj->nkm", g, ub, g)
            self.sites.append(dict(x=x, q=q, start=start, n=n,
                                   stabilizer=len(stabilizers), dp=dphase,
                                   at=at, ub=ub, nu=nu, du=d_exponent,
                                   aniso=anisotropic))
        if zr_fp:
            self.fp_index = len(p)
            p += [next((a.get("fp", -5) for a in atoms if a["element"] == "Zr"), -5)]
            lower += [-18]
            upper += [5]
        p += [-7.0]
        lower += [-20]
        upper += [5]
        self.p0 = np.array(p)
        self.bounds = np.array(lower), np.array(upper)
        self.last_p = None

    def compute(self, p):
        if self.last_p is not None and np.array_equal(p, self.last_p):
            return self.fc, self.jac
        fc = np.zeros(len(self.hkl))
        jac = np.zeros((len(self.hkl), len(p)))
        for site in self.sites:
            j, n, nu = site["start"], site["n"], site["nu"]
            pos = site["x"] + site["q"] @ p[j:j+n]
            xyz = OPS @ pos
            phase = 2*np.pi*(self.hkl @ xyz.T)
            is_zr = site["at"]["element"] == "Zr"
            fp = p[self.fp_index] if is_zr and self.fp_index is not None else site["at"].get("fp", 0)
            fac = ((self.f[site["at"]["element"]] + fp) * site["at"].get("occ", 1)
                   / site["stabilizer"])
            dw = np.exp(np.clip(np.einsum("nkm,m->nk", site["du"],
                                          p[j+n:j+n+nu]), -500, 50))
            cosdw = np.cos(phase)*dw
            f = cosdw.sum(axis=1)*fac
            fc += f
            if n:
                jac[:, j:j+n] = -np.einsum("nk,nkj->nj",
                                           np.sin(phase)*dw, site["dp"])*fac[:, None]
            jac[:, j+n:j+n+nu] = np.einsum("nk,nkm->nm", cosdw, site["du"])*fac[:, None]
            if is_zr and self.fp_index is not None:
                jac[:, self.fp_index] += (cosdw.sum(axis=1) *
                                          site["at"].get("occ", 1) /
                                          site["stabilizer"])
        self.last_p, self.fc, self.jac = p.copy(), fc, jac
        return fc, jac

    def output_atoms(self, p):
        atoms = []
        for site in self.sites:
            j, n, nu = site["start"], site["n"], site["nu"]
            at = site["at"].copy()
            at["xyz"] = (site["x"] + site["q"] @ p[j:j+n]).tolist()
            umat = np.einsum("m,mij->ij", p[j+n:j+n+nu], site["ub"])
            at["u"] = float(np.trace(umat)/3)
            if site["aniso"]:
                at["u_cart"] = umat.tolist()
                at["u_eigenvalues"] = np.linalg.eigvalsh(umat).tolist()
            else:
                at.pop("u_cart", None)
            at["multiplicity"] = 24 // site["stabilizer"]
            at["position_dof"] = n
            if at["element"] == "Zr" and self.fp_index is not None:
                at["fp"] = float(p[self.fp_index])
            atoms.append(at)
        return atoms

    def adp_restraints(self, p):
        res, jac = [], []
        for site in self.sites:
            if not site["aniso"]:
                continue
            j, n, nu = site["start"], site["n"], site["nu"]
            umat = np.einsum("m,mij->ij", p[j+n:j+n+nu], site["ub"])
            vals, vecs = np.linalg.eigh(umat)
            for k, val in enumerate(vals):
                row = np.zeros(len(p))
                if val < .002:
                    r = (val-.002)/.001
                    row[j+n:j+n+nu] = np.einsum("i,mij,j->m",
                                                 vecs[:, k], site["ub"],
                                                 vecs[:, k])/.001
                elif val > .6:
                    r = (val-.6)/.02
                    row[j+n:j+n+nu] = np.einsum("i,mij,j->m",
                                                 vecs[:, k], site["ub"],
                                                 vecs[:, k])/.02
                else:
                    r = 0.0
                res.append(r)
                jac.append(row)
        return np.array(res), np.array(jac).reshape(-1, len(p))


def metrics(data, fc, scale):
    fo = np.sqrt(np.maximum(data["i"], 0))
    fcalc = np.abs(fc) * np.sqrt(scale)
    observed = data["i"] > 2*data["sig"]
    r1 = np.abs(fo-fcalc)
    w = 1/data["sig"]**2
    return dict(n=len(fo), observed=int(observed.sum()),
                scale=float(scale),
                r1_observed=float(r1[observed].sum()/fo[observed].sum()),
                r1_all=float(r1.sum()/fo.sum()),
                wr2_sigma=float(np.sqrt(np.sum(w*(data["i"]-scale*fc**2)**2) /
                                         np.sum(w*data["i"]**2))))


def make_maps(data, atoms, scale, name, dmin=1.2):
    fc = structure_factor(data["hkl"], atoms)
    keep = data["d"] >= dmin
    fo = np.sqrt(np.maximum(data["i"], 0)/scale)
    for kind, coeff in [("fo", fo*np.sign(fc)),
                        ("diff", fo*np.sign(fc)-fc),
                        ("2fofc", 2*fo*np.sign(fc)-fc)]:
        density = fourier_map(data["hkl"][keep], coeff[keep])
        np.save(ROOT / f"{name}_{kind}.npy", density)
        peaks = unique_peaks(density, n=100, min_distance=.6)
        rows = [dict(height=v, xyz=snap_site(p, .12).tolist()) for v, p in peaks]
        (ROOT / f"{name}_{kind}_peaks.json").write_text(json.dumps(rows, indent=2))
    return metrics(data, fc, scale)


def run(args):
    t = time.time()
    atoms = (json.loads((ROOT / args.input).read_text()) if args.input else
             initial_atoms())
    atoms = [a for a in atoms if a["element"] != "H"]
    if args.fix_zr_fp is not None:
        if args.zr_fp:
            raise ValueError("Cannot both fix and refine the Zr scattering offset")
        for a in atoms:
            if a["element"] == "Zr":
                a["fp"] = args.fix_zr_fp
    (ROOT / f"{args.name}_initial.json").write_text(json.dumps(atoms, indent=2))
    all_data = load_data()
    select = (all_data["d"] >= args.dmin) & (all_data["d"] <= args.dmax)
    if args.holdout:
        from solvent_density import free_flags
        select &= ~free_flags(all_data)
    data = {k: v[select] for k, v in all_data.items()}
    model = Model(atoms, data["hkl"], data["s2"], aniso=args.aniso, zr_fp=args.zr_fp)
    if args.restraints:
        from restraints import Restraints
        geometry = Restraints(model)
    else:
        geometry = None
    def extra(p):
        r, j = model.adp_restraints(p)
        if geometry is not None:
            rg, jg = geometry.calculate(p)
            return np.r_[r, rg], np.vstack([j, jg])
        return r, j
    p = model.p0.copy()
    def h_contribution(ats):
        if not args.hydrogen:
            return np.zeros(len(data["hkl"]))
        from hydrogens import add_hydrogens
        hats = [a for a in add_hydrogens(ats) if a["element"] == "H"]
        return structure_factor(data["hkl"], hats)
    fh = h_contribution(atoms)
    fc = model.compute(p)[0] + fh
    obs = data["i"] > 2*data["sig"]
    amp = np.sqrt(np.maximum(data["i"], 0))
    k = np.sum(amp[obs]*np.abs(fc[obs])) / np.sum(fc[obs]**2)
    p[-1] = np.log(k*k)
    log = []
    for cycle in range(args.cycles):
        fh = h_contribution(model.output_atoms(p))
        fc = model.compute(p)[0] + fh
        ic = np.exp(p[-1])*fc**2
        pp = (np.maximum(data["i"], 0)+2*ic)/3
        w = 1 / np.sqrt(data["sig"]**2 + (args.weight*pp)**2)
        def residual(p):
            fc, _ = model.compute(p)
            fc = fc + fh
            return np.r_[(np.exp(p[-1])*fc**2-data["i"])*w,
                         extra(p)[0]]
        def jacobian(p):
            fc, jac = model.compute(p)
            fc = fc + fh
            ans = jac * (2*np.exp(p[-1])*fc*w)[:, None]
            ans[:, -1] = np.exp(p[-1])*fc**2*w
            return np.vstack([ans, extra(p)[1]])
        fit = least_squares(residual, p, jac=jacobian, bounds=model.bounds,
                            x_scale="jac", max_nfev=args.max_nfev,
                            ftol=2e-9, xtol=2e-9, gtol=2e-6)
        p = fit.x
        ats = model.output_atoms(p)
        m = metrics(data, model.compute(p)[0] + fh, np.exp(p[-1]))
        m.update(cycle=cycle, cost=fit.cost, nfev=fit.nfev, optimality=fit.optimality)
        log.append(m)
        print(json.dumps(m), flush=True)
        (ROOT / f"{args.name}_atoms.json").write_text(json.dumps(ats, indent=2))
    scale = np.exp(p[-1])
    atoms = model.output_atoms(p)
    if args.hydrogen:
        from hydrogens import add_hydrogens
        atoms = add_hydrogens(atoms)
    (ROOT / f"{args.name}_atoms.json").write_text(json.dumps(atoms, indent=2))
    m_all = make_maps(all_data, atoms, scale, args.name)
    print("ALL DATA", json.dumps(m_all), flush=True)
    holdout_metrics = None
    if args.holdout:
        from solvent_density import free_flags
        test = ((all_data["d"] >= args.dmin) & (all_data["d"] <= args.dmax)
                & free_flags(all_data))
        holdout_metrics = metrics({k:v[test] for k,v in all_data.items()},
                                  structure_factor(all_data["hkl"][test], atoms),
                                  scale)
        print("HOLDOUT", json.dumps(holdout_metrics), flush=True)
    for at in atoms:
        print(at["label"], at["xyz"], "U", at["u"], "mult", len(orbit(at["xyz"])))
    (ROOT / f"{args.name}_log.json").write_text(json.dumps(dict(
        cycles=log, all_data=m_all, dmin=args.dmin, weight=args.weight,
        dmax=args.dmax, aniso=args.aniso, zr_fp=args.zr_fp,
        restraints=args.restraints, n_parameters=len(p),
        hydrogen=args.hydrogen,
        fixed_zr_fp=args.fix_zr_fp,
        holdout=holdout_metrics,
        seconds=time.time()-t), indent=2))
    np.savez_compressed(ROOT / f"{args.name}_normal.npz", jac=fit.jac,
                        residual=fit.fun, parameters=p)
    print("SECONDS", time.time()-t, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--name", default="isotropic")
    parser.add_argument("--dmin", type=float, default=1.1)
    parser.add_argument("--dmax", type=float, default=100)
    parser.add_argument("--weight", type=float, default=.08)
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--max-nfev", type=int, default=120)
    parser.add_argument("--aniso", choices=["none", "heavy", "all"], default="none")
    parser.add_argument("--zr-fp", action="store_true")
    parser.add_argument("--restraints", action="store_true")
    parser.add_argument("--hydrogen", action="store_true")
    parser.add_argument("--holdout", action="store_true")
    parser.add_argument("--fix-zr-fp", type=float)
    run(parser.parse_args())
