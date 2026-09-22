"""Independent small CIF reader for round-trip validation of our own exports."""
import json
import re
import shlex
import numpy as np
from framework import ROOT, CELL, INV, OPS, orbit, load_data, structure_factor
from refine import metrics


def tokenize_cif(text):
    lines = iter(text.splitlines())
    for line in lines:
        if line.startswith(";"):
            paragraph = [line[1:]]
            for part in lines:
                if part.startswith(";"):
                    break
                paragraph.append(part)
            else:
                raise ValueError("Unterminated CIF text field")
            yield "\n".join(paragraph)
        else:
            yield from shlex.split(line, comments=True, posix=True)


def parse_cif(path):
    tokens = list(tokenize_cif(path.read_text()))
    fields, loops = {}, []
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t.startswith("data_"):
            i += 1
        elif t == "loop_":
            i += 1
            keys = []
            while i < len(tokens) and tokens[i].startswith("_"):
                keys.append(tokens[i])
                i += 1
            values = []
            while i < len(tokens) and not (tokens[i].startswith("_") or
                                           tokens[i] == "loop_" or
                                           tokens[i].startswith("data_")):
                values.append(tokens[i])
                i += 1
            if not keys or len(values) % len(keys):
                raise ValueError("Invalid CIF loop row width")
            loops.append([dict(zip(keys, values[j:j+len(keys)]))
                          for j in range(0, len(values), len(keys))])
        elif t.startswith("_"):
            if t in fields:
                raise ValueError("Duplicate CIF field: "+t)
            fields[t] = tokens[i+1]
            i += 2
        else:
            raise ValueError("Unexpected token: "+t[:50])
    return fields, loops


def value(text):
    return float(re.sub(r"\(\d+\)$", "", text))


def read_structure(path):
    fields, loops = parse_cif(path)
    loop = next(x for x in loops if "_atom_site_label" in x[0])
    tensors = next(x for x in loops if "_atom_site_aniso_label" in x[0])
    dispersions = next(x for x in loops if "_atom_type_symbol" in x[0])
    fp = {x["_atom_type_symbol"]:value(x["_atom_type_scat_dispersion_real"])
          for x in dispersions}
    d = np.diag(np.linalg.norm(INV, axis=1))
    u = {}
    for row in tensors:
        mat = np.zeros((3, 3))
        for (i, j), tag in zip([(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)],
                              ["11", "22", "33", "23", "13", "12"]):
            mat[i, j] = mat[j, i] = value(row["_atom_site_aniso_U_"+tag])
        u[row["_atom_site_aniso_label"]] = CELL @ d @ mat @ d @ CELL.T
    atoms = []
    for row in loop:
        label = row["_atom_site_label"]
        element = row["_atom_site_type_symbol"]
        atom = dict(label=label, element=element,
                    xyz=[value(row["_atom_site_fract_"+axis]) for axis in "xyz"],
                    u=value(row["_atom_site_U_iso_or_equiv"]),
                    occ=value(row["_atom_site_occupancy"]), fp=fp[element])
        if label in u:
            atom["u_cart"] = u[label]
        if len(orbit(atom["xyz"])) != int(row["_atom_site_symmetry_multiplicity"]):
            raise ValueError("Exported special-position multiplicity mismatch: "+label)
        atoms.append(atom)
    return fields, atoms


def run():
    data = load_data()
    results = []
    for name, filename in [("final_all_data", "nu1000_framework.cif"),
                           ("highangle_check", "nu1000_highangle_check.cif"),
                           ("no_dispersion", "nu1000_no_dispersion_check.cif")]:
        fields, atoms = read_structure(ROOT / filename)
        reference = json.loads((ROOT / f"{name}_atoms.json").read_text())
        log = json.loads((ROOT / f"{name}_log.json").read_text())
        fc = structure_factor(data["hkl"], atoms)
        expected = structure_factor(data["hkl"], reference)
        max_error = float(np.max(np.abs(fc-expected)))
        select = (data["d"] >= log["dmin"]) & (data["d"] <= log["dmax"])
        m = metrics({k:v[select] for k,v in data.items()}, fc[select],
                    log["all_data"]["scale"])
        r_error = abs(m["r1_observed"]-value(fields["_refine_ls_R_factor_gt"]))
        tensor_error = max(float(np.max(np.abs(np.array(a["u_cart"])-b["u_cart"])))
                           for a, b in zip(reference, atoms) if "u_cart" in a)
        row = dict(file=filename, atom_sites=len(atoms),
                   max_Fcalc_roundtrip_error=max_error, R1_roundtrip_error=r_error,
                   max_cartesian_tensor_roundtrip_error=tensor_error)
        row["passed"] = max_error < 1e-3 and r_error < 2e-6 and tensor_error < 1e-7
        results.append(row)
    fields, loops = parse_cif(ROOT / "nu1000_calculated.fcf")
    rows = next(x for x in loops if "_refln_index_h" in x[0])
    if len(rows) != len(data["hkl"]):
        raise ValueError("Incorrect FCF reflection count")
    model = json.loads((ROOT / "final_all_data_atoms.json").read_text())
    scale = json.loads((ROOT / "final_all_data_log.json").read_text())["all_data"]["scale"]
    ic = scale*structure_factor(data["hkl"], model)**2
    i_error, sig_error, ic_error = 0., 0., 0.
    for index, row in enumerate(rows):
        assert [int(row["_refln_index_"+s]) for s in "hkl"] == data["hkl"][index].tolist()
        i_error = max(i_error, abs(value(row["_refln_F_squared_meas"])-data["i"][index]))
        sig_error = max(sig_error, abs(value(row["_refln_F_squared_sigma"])-data["sig"][index]))
        ic_error = max(ic_error, abs(value(row["_refln_F_squared_calc"])-ic[index]))
    result = dict(cif_roundtrips=results, fcf_reflections=len(rows),
                  fcf_I_error=i_error, fcf_sigma_error=sig_error,
                  fcf_Icalc_error=ic_error,
                  passed=bool(all(x["passed"] for x in results)
                              and max(i_error, sig_error, ic_error) < 1e-4),
                  note="Self-written format and numerical checks, not checkCIF or a standard validation service.")
    (ROOT / "export_validation.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise RuntimeError("Export validation failed")


if __name__ == "__main__":
    run()
