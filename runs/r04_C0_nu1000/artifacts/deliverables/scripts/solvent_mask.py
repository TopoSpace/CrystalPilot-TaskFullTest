"""Use the installed smtbx solvent-mask engine and export SHELXL ABIN factors.

LIST 6 Fo^2 has SHELXL's anomalous contribution removed. The real f-prime
values must still be included; this is verified against its tabulated Fc.
The source HKL is never modified.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from cctbx import xray
from cctbx.array_family import flex
from iotbx import reflection_file_reader, ccp4_map
from smtbx import masks

parser = argparse.ArgumentParser()
parser.add_argument("model")
parser.add_argument("output")
parser.add_argument("--probe", type=float, default=1.3)
parser.add_argument("--shrink", type=float, default=1.3)
parser.add_argument("--grid", type=float, default=0.25)
parser.add_argument("--complete", action="store_true")
parser.add_argument("--cycles", type=int, default=30)
args = parser.parse_args()
base = Path(args.model).resolve()
output = Path(args.output).resolve()
text = base.with_suffix(".res").read_text()
xs = xray.structure.from_shelx(filename=str(base.with_suffix(".res")))
xs.scattering_type_registry(table="it1992")
disp = {
    line.split()[1].upper(): float(line.split()[2])
    for line in text.splitlines() if line.startswith("DISP ")
}
for sc in xs.scatterers():
    sc.fp = disp.get(sc.scattering_type.upper(), 0)
    sc.fdp = 0
    sc.flags.set_use_fp_fdp(True)
arrays = reflection_file_reader.any_reflection_file(
    str(base.with_suffix(".fcf"))
).as_miller_arrays()
obs, fc_listed = arrays
obs = obs.customized_copy(crystal_symmetry=xs)
fc = obs.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
errors = flex.abs(flex.abs(fc.data()) - flex.abs(fc_listed.data()))
print("FC cross-check: mean/max absolute discrepancy:",
      flex.mean(errors), flex.max(errors), flush=True)
if flex.mean(errors) > 0.1:
    raise ValueError("FC calculation does not match SHELXL; do not mask")
obs = obs.map_to_asu().sort()
mask = masks.mask(xs, obs, use_set_completion=args.complete)
mask.compute(
    solvent_radius=args.probe, shrink_truncation_radius=args.shrink,
    grid_step=args.grid, resolution_factor=None, use_space_group_symmetry=True,
)
print("Grid:", mask.crystal_gridding.n_real(), flush=True)
print("Voids:", mask.n_voids(), flush=True)
fmask = mask.structure_factors(max_cycles=args.cycles)
fmask = fmask.common_set(obs)
assert fmask.size() == obs.size()
print("Mask volume:", mask.solvent_accessible_volume)
print("Engine total electron count:", mask.f_000_s)
with output.with_suffix(".fab").open("w") as stream:
    for hkl, z in zip(fmask.indices(), fmask.data()):
        stream.write("%4d %4d %4d %16.8f %16.8f\n" %
                     (*hkl, z.real, z.imag))
    stream.write("0 0 0 0 0\n")
    stream.write("smtbx mask of unassigned pore density; no solvent identities assigned.\n")
    stream.write("Source: " + str(base) + "\n")
fmask.as_mtz_dataset(column_root_label="Fmask").mtz_object().write(
    str(output.with_suffix(".mtz"))
)
ccp4_map.write_ccp4_map(
    file_name=str(output.with_suffix(".ccp4")),
    unit_cell=xs.unit_cell(), space_group=xs.space_group(),
    map_data=mask.mask.data.as_double(),
    labels=flex.std_string(["smtbx solvent-accessible region mask"]),
)
summary = dict(
    source=str(base), engine="smtbx.masks.mask",
    use_set_completion=args.complete, solvent_radius=args.probe,
    shrink_truncation_radius=args.shrink, grid_step=args.grid,
    grid=list(mask.crystal_gridding.n_real()), max_cycles=args.cycles,
    fcalc_check_mean=flex.mean(errors), fcalc_check_max=flex.max(errors),
    void_volume_A3=mask.solvent_accessible_volume,
    void_fraction=mask.solvent_accessible_volume / xs.unit_cell().volume(),
    electrons_per_cell=mask.f_000_s,
    n_voids=mask.n_voids(),
    n_reflections=fmask.size(), scale_factor=mask.scale_factor,
    caveat="Unassigned density; electron count is model-dependent, not a solvent formula.",
)
output.with_suffix(".mask.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
