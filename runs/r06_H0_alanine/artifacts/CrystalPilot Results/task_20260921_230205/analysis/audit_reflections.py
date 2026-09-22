"""Read-only diagnostics of this task's DIALS reflection tables."""

import json
from pathlib import Path

import numpy as np
from dials.array_family import flex
from dxtbx.format.FormatROD import FormatROD
from dxtbx.model.experiment_list import ExperimentListFactory


ROOT = Path(__file__).resolve().parents[3]
FRAMES = ROOT / ".crystalpilot" / "frames"
TARGETS = {
    (0, 2, 2),
    (0, 2, 3),
    (0, 2, 4),
    (0, 3, 1),
    (0, 3, 3),
    (1, 0, 2),
    (1, 3, 1),
}
FIELDS = [
    "id",
    "panel",
    "miller_index",
    "xyzcal.px",
    "xyzobs.px.value",
    "intensity.sum.value",
    "intensity.sum.variance",
    "intensity.prf.value",
    "intensity.prf.variance",
    "intensity.scale.value",
    "intensity.scale.variance",
    "inverse_scale_factor",
    "partiality",
    "profile.correlation",
    "flags",
]

for stem in ["integrated", "scaled"]:
    path = FRAMES / (stem + ".refl")
    table = flex.reflection_table.from_file(str(path))
    print("\nTABLE", stem, "ROWS", len(table))
    print("FLAG ENUM", table.flags.names)
    for row in table.rows():
        hkl = tuple(abs(i) for i in row["miller_index"])
        if hkl in TARGETS:
            record = {key: row[key] for key in FIELDS if key in row}
            record["flag_names"] = [
                name
                for name, value in table.flags.names.items()
                if int(value) and (int(row["flags"]) & int(value)) == int(value)
            ]
            print(json.dumps(record, default=str))
    expt_path = FRAMES / (stem + ".expt")
    if expt_path.exists():
        experiments = ExperimentListFactory.from_json_file(
            str(expt_path), check_format=False
        )
        for index, experiment in enumerate(experiments):
            print("EXPERIMENT", index, experiment.imageset.get_template())
            print("CELL", experiment.crystal.get_unit_cell())
            print("SG", experiment.crystal.get_space_group().info())

# Read original pixels at predicted spots; no reflection is edited or excluded.
experiments = ExperimentListFactory.from_json_file(
    str(FRAMES / "integrated.expt"), check_format=False
)
table = flex.reflection_table.from_file(str(FRAMES / "integrated.refl"))
print("\nRAW PIXEL PATCHES (11 x 11 pixels, central image)")
for row in table.rows():
    hkl = tuple(abs(i) for i in row["miller_index"])
    if hkl not in TARGETS or row["partiality"] < 0.9:
        continue
    experiment = experiments[row["id"]]
    image_set = experiment.imageset
    x, y, z = row["xyzcal.px"]
    frame_index = int(z)
    if not 0 <= frame_index < len(image_set):
        continue
    # Use the supplied ROD decoder without initializing fresh instrument models.
    reader = FormatROD.__new__(FormatROD)
    reader._image_file = image_set.get_path(frame_index)
    reader._txt_header = FormatROD._read_ascii_header(reader._image_file)
    image = FormatROD.get_raw_data(reader).as_numpy_array()
    offset_x, offset_y = experiment.detector[row["panel"]].get_raw_image_offset()
    x += offset_x
    y += offset_y
    ix, iy = int(x), int(y)
    patch = image[iy - 5 : iy + 6, ix - 5 : ix + 6]
    print(
        json.dumps(
            {
                "hkl": row["miller_index"],
                "original_image": image_set.get_path(frame_index),
                "panel": row["panel"],
                "predicted_xy": [x, y],
                "min": int(patch.min()),
                "max": int(patch.max()),
                "mean": float(patch.mean()),
                "median": float(np.median(patch)),
                "zero_fraction": float(np.mean(patch == 0)),
                "negative_fraction": float(np.mean(patch < 0)),
            }
        )
    )
