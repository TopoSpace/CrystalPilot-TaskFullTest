"""Check the deliverable CIF, RES, HKL and FAB without relying on visual inspection."""
import argparse
import io
import json
from pathlib import Path
import hashlib
import gemmi
import numpy as np
from cctbx import xray, miller
from cctbx.array_family import flex
from iotbx import reflection_file_reader

parser = argparse.ArgumentParser()
parser.add_argument("directory")
args = parser.parse_args()
directory = Path(args.directory).resolve()
base = directory / "nu1000_framework"
block = gemmi.cif.read_file(str(base.with_suffix(".cif"))).sole_block()
assert block.find_value("_cell_measurement_temperature") == "?"
assert block.find_value("_diffrn_ambient_temperature") == "?"
assert block.find_value("_space_group_IT_number") == "191"
text = base.with_suffix(".res").read_text()
xs = xray.structure.from_shelx(file=io.StringIO(
    "\n".join(l for l in text.splitlines() if not l.startswith("ABIN"))))
xs.scattering_type_registry(table="it1992")
counts = {}
for sc in xs.scatterers():
    counts[sc.scattering_type] = counts.get(sc.scattering_type, 0) + sc.multiplicity()*sc.occupancy
assert counts == {"Zr": 18.0, "O": 96.0, "C": 264.0, "H": 132.0}
disp = {l.split()[1].upper(): float(l.split()[2])
        for l in text.splitlines() if l.startswith("DISP ")}
for sc in xs.scatterers():
    sc.fp = disp.get(sc.scattering_type.upper(), 0)
    sc.fdp = 0
    sc.flags.set_use_fp_fdp(True)
obs, fc_listed = reflection_file_reader.any_reflection_file(
    str(base.with_suffix(".fcf"))).as_miller_arrays()
fc_listed = fc_listed.customized_copy(crystal_symmetry=xs).map_to_asu().sort()
fc_atoms = fc_listed.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
indices, values = [], []
for line in base.with_suffix(".fab").read_text().splitlines():
    t = line.split()
    h = tuple(map(int, t[:3]))
    if h == (0, 0, 0):
        break
    indices.append(h)
    values.append(complex(float(t[3]), float(t[4])))
fm = miller.array(miller.set(xs, flex.miller_index(indices), False),
                  flex.complex_double(values)).map_to_asu().sort().common_set(fc_listed)
assert fm.size() == 4395 == obs.size()
assert list(fm.indices()) == list(fc_listed.indices())
atom_error = flex.mean(flex.abs(flex.abs(fc_atoms.data())-flex.abs(fc_listed.data())))
total_error = flex.mean(flex.abs(flex.abs(fc_atoms.data()+fm.data())-flex.abs(fc_listed.data())))
# LIST 6 includes the ABIN contribution in Fc for this SHELXL build.
assert min(atom_error, total_error) < 0.02
parsed_cif = xray.structure.from_cif(file_path=str(base.with_suffix(".cif")))
assert len(parsed_cif) == 1
doc_model = next(iter(parsed_cif.values()))
assert len(doc_model.scatterers()) == len(xs.scatterers())
coordinate_rounding = []
for a, b in zip(doc_model.scatterers(), xs.scatterers()):
    assert a.label.upper() == b.label.upper()
    delta = np.array(a.site)-b.site
    coordinate_rounding.append(float(np.max(np.abs(delta))))
    # CIF coordinates are rounded to the precision of their reported s.u.
    assert coordinate_rounding[-1] <= 5.1e-5
tags = (
    "_refine_ls_R_factor_gt", "_refine_ls_R_factor_all",
    "_refine_ls_wR_factor_ref", "_refine_ls_goodness_of_fit_ref",
    "_refine_ls_number_reflns", "_refine_ls_number_parameters",
    "_refine_ls_number_restraints", "_refine_ls_shift/su_max",
    "_refine_diff_density_max", "_refine_diff_density_min",
    "_diffrn_reflns_av_R_equivalents",
)
result = {
    "cif_and_res_agree_within_output_rounding": True,
    "max_CIF_RES_fractional_rounding": max(coordinate_rounding),
    "temperature_unknown": True,
    "space_group": "P 6/m m m (191)", "formula_per_cell": counts,
    "fab_reflections": fm.size(), "mean_Fc_error_atom_only": atom_error,
    "mean_Fc_error_with_ABIN": total_error,
    "refinement": {t: block.find_value(t) for t in tags},
    "hkl_sha256": hashlib.sha256(base.with_suffix(".hkl").read_bytes()).hexdigest(),
}
(directory / "verification.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
