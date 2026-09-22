"""Explicit chemical-geometry and displacement restraints for the found linker.

Only measured framework sites are used. No reference structure coordinates
are present. Distances are in Angstrom and displacement tensors in Angstrom^2.
"""
import itertools
import numpy as np
from framework import CELL, INV, OPS


class Restraints:
    def __init__(self, model):
        self.model = model
        self.bonds = []
        self.angles = []
        self.planes = []
        identity = np.eye(3)
        self.origin_points = [(i, identity, np.zeros(3)) for i in range(len(model.sites))]
        neighbours = {}
        atoms = model.atoms
        for i, a in enumerate(atoms):
            if a["element"] not in ("C", "O"):
                continue
            xi = np.array(a["xyz"])
            near = []
            for j, b in enumerate(atoms):
                if b["element"] not in ("C", "O") or a["element"] == b["element"] == "O":
                    continue
                seen = set()
                for op in OPS:
                    xj0 = op @ np.array(b["xyz"])
                    shift = np.round(xi - xj0)
                    xj = xj0 + shift
                    key = tuple(np.round(xj, 7))
                    if key in seen:
                        continue
                    seen.add(key)
                    d = np.linalg.norm(CELL @ (xi-xj))
                    if not .7 < d < 1.85:
                        continue
                    labels = {a["label"], b["label"]}
                    if "O" in (a["element"], b["element"]):
                        ideal, sig = 1.26, .025
                    elif labels in ({"C3", "C6"}, {"C7", "C12"}):
                        ideal, sig = 1.49, .03
                    else:
                        ideal, sig = 1.40, .03
                    point = (j, op, shift)
                    self.bonds.append((self.origin_points[i], point, ideal, sig))
                    near.append((point, ideal))
            if a["element"] == "C":
                neighbours[i] = near
        for i, near in neighbours.items():
            for (a, la), (b, lb) in itertools.combinations(near, 2):
                self.angles.append((a, b, np.sqrt(la*la+lb*lb+la*lb), .06))
        labels = {a["label"]: i for i, a in enumerate(atoms)}
        phenyl = [self.origin_points[labels[s]] for s in
                  ["C2", "C3", "C4", "C5", "C7", "C10"]]
        self.planes.append((phenyl, .035))
        # Expand the pyrene fragment around its measured molecular centre.
        centre = np.array([.5, .25, .5])
        pyrene = []
        seen = set()
        for label in ("C1", "C6", "C8", "C9", "C11"):
            j = labels[label]
            for op in OPS:
                x0 = op @ np.array(atoms[j]["xyz"])
                shift = np.round(centre-x0)
                x = x0+shift
                key = tuple(np.round(x, 7))
                if key not in seen and np.linalg.norm(CELL@(x-centre)) < 4.8:
                    seen.add(key)
                    pyrene.append((j, op, shift))
        if len(pyrene) != 16:
            raise ValueError(f"Pyrene plane has {len(pyrene)} atoms, expected 16")
        self.planes.append((pyrene, .035))
        carbox = [self.origin_points[labels["C12"]]]
        carbox += [p for p, _ in neighbours[labels["C12"]]]
        self.planes.append((carbox, .04))

    def point(self, item, p):
        i, op, shift = item
        site = self.model.sites[i]
        j, n = site["start"], site["n"]
        pos = site["x"] + site["q"] @ p[j:j+n]
        x = CELL @ (op @ pos + shift)
        jac = np.zeros((3, len(p)))
        jac[:, j:j+n] = CELL @ op @ site["q"]
        return x, jac

    def tensor(self, item, p):
        i, op, shift = item
        s = self.model.sites[i]
        j, n, nu = s["start"], s["n"], s["nu"]
        rot = CELL @ op @ INV
        basis = rot @ s["ub"] @ rot.T
        u = np.einsum("m,mij->ij", p[j+n:j+n+nu], basis)
        return u, (j+n, nu, basis)

    def calculate(self, p):
        residual, jacobian = [], []
        for a, b, target, sigma in self.bonds + self.angles:
            xa, ja = self.point(a, p)
            xb, jb = self.point(b, p)
            dx = xa-xb
            d = np.linalg.norm(dx)
            residual.append((d-target)/sigma)
            jacobian.append((dx/d) @ (ja-jb)/sigma)
        for group, sigma in self.planes:
            xj = [self.point(g, p) for g in group]
            xyz = np.array([x for x, _ in xj])
            centred = xyz-xyz.mean(axis=0)
            _, _, vh = np.linalg.svd(centred)
            normal = vh[-1]
            # Variable-projection derivative: sufficient for the least-squares
            # gradient; normal re-estimated at every residual evaluation.
            mean_j = np.mean([j for _, j in xj], axis=0)
            for x, j in xj:
                residual.append(float((x-xyz.mean(axis=0)) @ normal)/sigma)
                jacobian.append(normal @ (j-mean_j)/sigma)
        for a, b, _, _ in self.bonds:
            ua, (ia, na, ba) = self.tensor(a, p)
            ub, (ib, nb, bb) = self.tensor(b, p)
            xa, _ = self.point(a, self.model.p0)
            xb, _ = self.point(b, self.model.p0)
            direction = (xa-xb)/np.linalg.norm(xa-xb)
            row = np.zeros(len(p))
            row[ia:ia+na] += np.einsum("i,mij,j->m", direction, ba, direction)/.008
            row[ib:ib+nb] -= np.einsum("i,mij,j->m", direction, bb, direction)/.008
            residual.append(float(direction @ (ua-ub) @ direction)/.008)
            jacobian.append(row)
            # Weak neighbouring-atom tensor similarity stabilizes poorly
            # determined transverse components without making them isotropic.
            for i, j in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]:
                row = np.zeros(len(p))
                row[ia:ia+na] += ba[:, i, j]/.08
                row[ib:ib+nb] -= bb[:, i, j]/.08
                residual.append((ua[i, j]-ub[i, j])/.08)
                jacobian.append(row)
        return np.array(residual), np.array(jacobian)
