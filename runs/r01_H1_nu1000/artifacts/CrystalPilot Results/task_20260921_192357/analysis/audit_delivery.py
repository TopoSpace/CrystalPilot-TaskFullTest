"""Read-only checks of this delivery; no structure generation or refinement."""

import json
from pathlib import Path

import gemmi
from cctbx import crystal, miller
from cctbx.array_family import flex
from iotbx import reflection_file_reader


root = Path(__file__).resolve().parents[1]
block = gemmi.cif.read_file(str(root / "final.cif")).sole_block()
cell = [
    gemmi.cif.as_number(block.find_value("_cell_" + name))
    for name in (
        "length_a", "length_b", "length_c",
        "angle_alpha", "angle_beta", "angle_gamma",
    )
]
symmetry = crystal.symmetry(unit_cell=cell, space_group_symbol="P 6/m m m")
raw = reflection_file_reader.any_reflection_file(
    str(root / "final.hkl") + "=hklf4"
).as_miller_arrays(
    crystal_symmetry=symmetry, merge_equivalents=False
)[0]
raw_asu = raw.map_to_asu()
targets = [
    (-4, 34, 3), (-1, 26, 7),
    (0, 1, 0), (-1, 2, 0), (0, 2, 0), (0, 0, 1),
]
mapped_targets = miller.set(
    symmetry, flex.miller_index(targets), anomalous_flag=False
).map_to_asu().indices()

fcf = gemmi.cif.read_file(str(root / "final.fcf")).sole_block()
fcf_indices = flex.miller_index(
    [
        tuple(map(int, row))
        for row in fcf.find(
            ["_refln_index_h", "_refln_index_k", "_refln_index_l"]
        )
    ]
)
fcf_asu = set(
    miller.set(symmetry, fcf_indices, anomalous_flag=False)
    .map_to_asu()
    .indices()
)
counts = {tuple(h): 0 for h in mapped_targets}
for hkl in raw_asu.indices():
    if hkl in counts:
        counts[hkl] += 1

missing = []
for original, canonical in zip(targets, mapped_targets):
    missing.append(
        {
            "checkcif_hkl": original,
            "canonical_hkl": canonical,
            "d_A": symmetry.unit_cell().d(original),
            "raw_symmetry_equivalent_observations": counts[canonical],
            "in_fcf": canonical in fcf_asu,
        }
    )

atom_rows = [
    {
        "label": row[0],
        "element": row[1],
        "occupancy": gemmi.cif.as_number(row[2]),
        "Ueq_A2": gemmi.cif.as_number(row[3]),
    }
    for row in block.find(
        [
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_occupancy",
            "_atom_site_U_iso_or_equiv",
        ]
    )
]
selected = raw.resolution_filter(d_min=1.0)
result = {
    "scope": "Read-only analysis of final.cif/fcf/hkl; no model changed",
    "cell": cell,
    "raw_observations": raw.size(),
    "observations_d_ge_1_A": selected.size(),
    "negative_intensities_d_ge_1_A": int((selected.data() < 0).count(True)),
    "fcf_rows": len(fcf_indices),
    "missing_reflection_audit": missing,
    "atom_table": atom_rows,
}
print(json.dumps(result, ensure_ascii=False, indent=2))
