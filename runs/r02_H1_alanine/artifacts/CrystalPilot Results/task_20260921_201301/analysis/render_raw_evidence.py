"""Render original detector counts, without editing data or any model."""
import json
from pathlib import Path

import dxtbx
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[3]
FRAMES = ROOT / "inputs" / "alanine" / "frames"
CASES = [
    ("pgw240033_Mo_3_46.rodhypix", 0, 341, 381, "(0,2,-2): near-zero observation"),
    ("pgw240033_Mo_3_77.rodhypix", 1, 190, 382, "(0,-2,2): strong equivalent"),
    ("pgw240033_Mo_3_71.rodhypix", 0, 292, 441, "(1,3,-1): near-zero observation"),
    ("pgw240033_Mo_3_72.rodhypix", 0, 292, 322, "(-1,3,-1): strong equivalent"),
]
fig, axes = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
records = []
for ax, (filename, panel, x, y, title) in zip(axes.flat, CASES):
    image = dxtbx.load(str(FRAMES / filename))
    data = image.get_raw_data()[panel].as_numpy_array()
    x0, x1 = max(0, x - 75), min(data.shape[1], x + 76)
    y0, y1 = max(0, y - 75), min(data.shape[0], y + 76)
    ax.imshow(
        np.log1p(np.maximum(data[y0:y1, x0:x1], 0)),
        origin="lower", cmap="gray_r", vmin=0, vmax=6,
        extent=(x0, x1, y0, y1),
    )
    ax.add_patch(Rectangle((x - 5, y - 5), 10, 10, fill=False, edgecolor="red"))
    ax.set_title(title + "\n" + filename + f", panel {panel}", fontsize=9)
    ax.set_xlabel("Detector x (pixel)")
    ax.set_ylabel("Detector y (pixel)")
    box = data[max(0, y - 5):y + 6, max(0, x - 5):x + 6]
    records.append({
        "original_frame": str(FRAMES / filename), "panel": panel,
        "predicted_xy_rounded": [x, y], "box_max_counts": float(box.max()),
        "box_sum_counts": float(box.sum()), "box_zero_fraction": float((box == 0).mean()),
        "note": "Single-frame 11x11 raw box, not background-corrected integrated intensity.",
    })
fig.savefig(Path(__file__).with_name("raw_detector_evidence.png"), dpi=160)
Path(__file__).with_name("raw_detector_evidence.json").write_text(
    json.dumps(records, indent=2), encoding="utf-8"
)
print(json.dumps(records, indent=2))
