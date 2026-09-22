"""Chemically inferred aromatic riding H atoms; node protons remain unknown."""
import numpy as np
from framework import CELL, INV, orbit


def add_hydrogens(atoms):
    result = [a.copy() for a in atoms if a["element"] != "H"]
    for atom in atoms:
        if atom["element"] != "C":
            continue
        xyz = np.array(atom["xyz"])
        vectors = []
        for other in atoms:
            if other["element"] not in ("C", "O"):
                continue
            diff = orbit(other["xyz"])-xyz
            diff -= np.round(diff)
            cart = diff @ CELL.T
            dist = np.linalg.norm(cart, axis=1)
            for v, d in zip(cart, dist):
                if .5 < d < 1.8:
                    vectors.append(v/d)
        if len(vectors) != 2:
            continue
        direction = -np.sum(vectors, axis=0)
        direction /= np.linalg.norm(direction)
        pos = xyz + INV @ (direction*.95)
        result.append(dict(label="H"+atom["label"][1:], element="H",
                           xyz=pos.tolist(), u=1.2*atom["u"],
                           parent=atom["label"], calculated=True, occ=1.0))
    return result
