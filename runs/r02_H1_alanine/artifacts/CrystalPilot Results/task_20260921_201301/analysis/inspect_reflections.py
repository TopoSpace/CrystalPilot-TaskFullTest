"""Read-only audit of DIALS reflection flags and suspect low-angle data."""
import json
from pathlib import Path

from dials.array_family import flex


ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / ".crystalpilot" / "frames"
TARGETS = {(0, 2, 2), (0, 2, 3), (1, 3, 1), (1, 0, 2), (0, 3, 3)}
FLAG_NAMES = [
    "indexed", "strong", "integrated_sum", "integrated_prf",
    "bad_for_scaling", "excluded_for_scaling", "outlier_in_scaling",
    "scaled", "dont_integrate",
]


def inspect(name):
    table = flex.reflection_table.from_file(str(WORK / name))
    output = []
    for i, hkl in enumerate(table["miller_index"]):
        if tuple(abs(v) for v in hkl) not in TARGETS:
            continue
        row = {"row": i, "hkl": list(hkl)}
        for key in (
            "id", "panel", "xyzcal.px", "xyzobs.px.value", "partiality",
            "intensity.sum.value", "intensity.sum.variance",
            "intensity.prf.value", "intensity.prf.variance",
            "intensity.scale.value", "intensity.scale.variance",
            "inverse_scale_factor",
        ):
            if key in table:
                value = table[key][i]
                row[key] = list(value) if isinstance(value, tuple) else value
        row["flags"] = [
            name for name in FLAG_NAMES
            if hasattr(table.flags, name)
            and table.get_flags(getattr(table.flags, name))[i]
        ]
        output.append(row)
    return {"file": name, "total_rows": len(table), "selected": output}


result = [inspect("integrated.refl"), inspect("scaled.refl")]
destination = Path(__file__).with_name("suspect_reflections.json")
destination.write_text(json.dumps(result, indent=2), encoding="utf-8")
for section in result:
    print(section["file"], "rows", section["total_rows"])
    for row in section["selected"]:
        print(json.dumps(row))
