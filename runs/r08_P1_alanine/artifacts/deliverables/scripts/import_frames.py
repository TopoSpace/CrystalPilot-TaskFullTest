"""Inventory and import only this task's raw frames; inputs remain read-only."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from run_command import ROOT, run

frames = sorted((ROOT / 'inputs' / 'alanine' / 'frames').glob('*.rodhypix'))
manifest = [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in frames]
(ROOT / 'raw_frames_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print('Raw frame counts per sweep:', dict(Counter(p.stem.split('_')[-2] for p in frames)))
phil = ROOT / 'work' / 'dials' / 'raw_import.phil'
phil.write_text('input {\n' + ''.join('  directory = "' + p.as_posix() + '"\n' for p in frames) + '}\n', encoding='utf-8')
sys.exit(run('work/dials', '03_import_phil', ['dials.import', 'raw_import.phil']))
