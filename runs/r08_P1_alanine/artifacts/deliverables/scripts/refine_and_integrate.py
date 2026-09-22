"""Refine independent scan orientations with a common static cell, then integrate."""
from copy import deepcopy
from dxtbx.model.experiment_list import ExperimentListFactory
from run_command import ROOT, run

work = ROOT / 'work' / 'dials'
experiments = ExperimentListFactory.from_json_file(str(work / 'refined_static.expt'), check_format=False)
for e in experiments:
    e.crystal = deepcopy(e.crystal)
experiments.as_file(str(work / 'scan_input.expt'))
commands = [
    ('09_refine_scan', ['dials.refine', 'scan_input.expt', 'refined_static.refl', 'scan_varying=True', 'crystal.fix=cell', 'detector.fix=all', 'beam.fix=all', 'interval_width_degrees=36', 'output.experiments=refined.expt', 'output.reflections=refined.refl', 'output.log=refined_scan.log']),
    ('10_integrate', ['dials.integrate', 'refined.expt', 'refined.refl', 'nproc=4', 'prediction.d_min=0.65']),
    ('11_symmetry', ['dials.symmetry', 'integrated.expt', 'integrated.refl', 'systematic_absences.small_molecule=True', 'output.html=None']),
]
for label, command in commands:
    code = run('work/dials', label, command)
    if code:
        raise SystemExit(code)
