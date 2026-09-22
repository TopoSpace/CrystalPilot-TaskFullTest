"""Numerical audit of a DIALS experiment/reflection pair (no plots)."""
import json
from pathlib import Path
import sys
import numpy as np
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex

stem = Path(sys.argv[1])
expts = ExperimentListFactory.from_json_file(str(stem.with_suffix('.expt')), check_format=False)
refl = flex.reflection_table.from_file(str(stem.with_suffix('.refl')))
report = {'experiments': [], 'n_reflections': len(refl), 'columns': list(refl.keys())}
for i, e in enumerate(expts):
    r = refl.select(refl['id'] == i)
    item = {'id': i, 'n_reflections': len(r), 'unit_cell': e.crystal.get_unit_cell().parameters(), 'space_group': str(e.crystal.get_space_group().info()), 'images': e.scan.get_image_range()}
    try:
        item['unit_cell_esd'] = e.crystal.get_cell_parameter_sd()
        item['volume_esd'] = e.crystal.get_cell_volume_sd()
    except (AttributeError, RuntimeError):
        pass
    if 'xyzcal.mm' in r and 'xyzobs.mm.value' in r:
        delta = r['xyzcal.mm'].as_numpy_array() - r['xyzobs.mm.value'].as_numpy_array()
        inlier = r.get_flags(r.flags.used_in_refinement).as_numpy_array()
        if np.any(inlier):
            item['rmsd_inliers_mm_mm_rad'] = np.sqrt(np.mean(delta[inlier]**2, axis=0)).tolist()
            item['refinement_inliers'] = int(inlier.sum())
    if 'partiality' in r:
        item['partiality_quantiles'] = np.quantile(r['partiality'].as_numpy_array(), [0, .1, .5, .9, 1]).tolist()
    report['experiments'].append(item)
text = json.dumps(report, indent=2)
stem.with_suffix('.audit.json').write_text(text, encoding='utf-8')
print(text)
