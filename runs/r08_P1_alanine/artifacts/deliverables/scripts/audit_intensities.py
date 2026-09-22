"""Recompute merging statistics from exported signed HKLF4 intensities.
Using HKLF4 avoids the MTZ reader's mapping to a non-anomalous reciprocal ASU.
"""
import io
import json
from pathlib import Path
import sys
from iotbx import merging_statistics
from cctbx.array_family import flex
from cctbx import crystal, miller
from dxtbx.model.experiment_list import ExperimentListFactory

work = Path(sys.argv[1]).resolve()
expt = ExperimentListFactory.from_json_file(str(work / 'scaled.expt'), check_format=False)[0]
cs = crystal.symmetry(unit_cell=expt.crystal.get_unit_cell(), space_group=expt.crystal.get_space_group())
indices, data, sigmas = [], [], []
for line in (work / 'alanine.hkl').read_text().splitlines():
    if len(line) < 28:
        continue
    hkl = tuple(int(line[i:i+4]) for i in (0,4,8))
    if hkl == (0,0,0):
        break
    indices.append(hkl)
    data.append(float(line[12:20]))
    sigmas.append(float(line[20:28]))
i_obs = miller.array(miller.set(cs, flex.miller_index(indices), anomalous_flag=True), data=flex.double(data), sigmas=flex.double(sigmas)).set_observation_type_xray_intensity().set_info(miller.array_info(source_type='shelx_hklf', labels=['I','SIGI']))
print('Unmerged signed HKLF4 records:', i_obs.size())
report = {'source': 'alanine.hkl', 'label': i_obs.info().label_string(), 'cell': i_obs.unit_cell().parameters(), 'space_group': str(i_obs.space_group_info()), 'statistics': {}}
out = io.StringIO()
for cutoff in [None, 0.73, 0.75, 0.77, 0.80, 0.84]:
    name = 'full' if cutoff is None else str(cutoff)
    stats = merging_statistics.dataset_statistics(i_obs, d_min=cutoff, anomalous=False, n_bins=10, sigma_filtering=None, use_internal_variance=False)
    report['statistics'][name] = stats.as_dict()
    out.write('\n### Resolution cutoff ' + name + ' angstrom\n')
    stats.show(out=out)
    stats.overall.show_summary()
for cutoff in [None, 0.80]:
    name = 'full' if cutoff is None else str(cutoff)
    a = i_obs.resolution_filter(d_min=cutoff).as_anomalous_array().eliminate_sys_absent().merge_equivalents(use_internal_variance=False).array()
    minus = a.anomalous_differences()
    report['statistics']['anomalous_' + name] = {'unique': a.size(), 'completeness': a.completeness(), 'friedel_pairs': minus.size(), 'mean_deltaI_over_sigma_abs': flex.mean(flex.abs(minus.data()) / minus.sigmas())}
(work / 'intensity_audit.json').write_text(json.dumps(report, indent=2, default=str), encoding='utf-8')
(work / 'intensity_audit.txt').write_text(out.getvalue(), encoding='utf-8')
print('Anomalous summaries:', {k: v for k,v in report['statistics'].items() if k.startswith('anomalous')})
