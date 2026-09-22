"""Select the six production scans (616 frames), leaving pre-experiment scans out."""
import json
from dxtbx.model.experiment_list import ExperimentListFactory, ExperimentList
from run_command import ROOT, run

work = ROOT / 'work' / 'dials'
experiments = ExperimentListFactory.from_json_file(str(work / 'imported.expt'), check_format=False)
production = ExperimentList([e for e in experiments if not __import__('pathlib').Path(e.imageset.get_path(0)).name.startswith('pre_')])
production.as_file(str(work / 'production.expt'))
metadata = []
for n, e in enumerate(production):
    entry = {'id': n, 'template': e.imageset.get_template(), 'images': e.scan.get_image_range(), 'oscillation': e.scan.get_oscillation(), 'beam': e.beam.to_dict(), 'goniometer': e.goniometer.to_dict(), 'detector': e.detector.to_dict()}
    metadata.append(entry)
    print(n, entry['images'], entry['oscillation'], 'wavelength', e.beam.get_wavelength(), 'panels', len(e.detector))
(work / 'imported_geometry.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
assert sum(len(e.imageset) for e in production) == 616
raise SystemExit(run('work/dials', '04_find_spots', ['dials.find_spots', 'production.expt', 'nproc=4']))
