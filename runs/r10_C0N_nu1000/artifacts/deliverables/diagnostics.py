"""Numerical model diagnostics independent of the optimizer."""
import argparse
import json
from collections import Counter
import numpy as np
from framework import ROOT, CELL, orbit, load_data, structure_factor
from inspect_maps import distance
from refine import metrics


def run(name):
    atoms = json.loads((ROOT / f"{name}_atoms.json").read_text())
    stats = json.loads((ROOT / f"{name}_log.json").read_text())
    k = stats["all_data"]["scale"]
    data = load_data()
    fc = structure_factor(data["hkl"], atoms)
    out = []
    for lo, hi in [(.69, 1), (1, 1.1), (1.1, 1.2), (1.2, 1.5),
                   (1.5, 2), (2, 3), (3, 5), (5, 100)]:
        sel = (data["d"] >= lo) & (data["d"] < hi)
        m = metrics({key: val[sel] for key, val in data.items()}, fc[sel], k)
        m.update(dmin=lo, dmax=hi)
        out.append(m)
        print(lo, hi, json.dumps(m))
    print("LARGEST INTENSITY RESIDUALS")
    ids = np.argsort(np.abs(data["i"]-k*fc**2))[::-1][:40]
    for ix in ids:
        print(data["hkl"][ix], f'd={data["d"][ix]:.2f}',
              f'I={data["i"][ix]:.2f} Ic={k*fc[ix]**2:.2f}',
              f'sig={data["sig"][ix]:.2f}')
    print("GEOMETRY")
    geo = []
    for at in atoms:
        p = np.array(at["xyz"])
        neighbours = []
        for other in atoms:
            dd = distance(p, orbit(other["xyz"]))
            limit = 2.65 if "Zr" in (at["element"], other["element"]) else 1.85
            for dist in dd[(dd > .2) & (dd < limit)]:
                neighbours.append([other["label"], round(float(dist), 4)])
        print(at["label"], sorted(neighbours, key=lambda x:x[1]),
              "Ueig", np.round(np.linalg.eigvalsh(at.get("u_cart",
                                                  np.eye(3)*at["u"])), 4))
        geo.append(dict(label=at["label"], neighbours=neighbours))
    formula = Counter()
    for at in atoms:
        formula[at["element"]] += len(orbit(at["xyz"]))*at.get("occ", 1)
    print("CELL CONTENT", dict(formula))
    (ROOT / f"{name}_diagnostics.json").write_text(json.dumps(
        dict(shells=out, geometry=geo, cell_content=dict(formula)), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    run(parser.parse_args().name)
