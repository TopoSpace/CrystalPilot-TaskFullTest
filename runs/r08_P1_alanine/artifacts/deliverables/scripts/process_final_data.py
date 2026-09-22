"""Final raw-data reduction using the jointly refined common cell and panel geometry."""
from run_command import run

commands = [
    ('13_integrate_joint', ['dials.integrate', '../dials/refined_joint.expt', '../dials/refined_joint.refl', 'nproc=4', 'prediction.d_min=0.65']),
    ('14_symmetry', ['dials.symmetry', 'integrated.expt', 'integrated.refl', 'output.html=None']),
    ('15_scale', ['dials.scale', 'symmetrized.expt', 'symmetrized.refl', '../dials/scale.phil']),
    ('16_export_shelx', ['dials.export', 'scaled.expt', 'scaled.refl', 'format=shelx', 'intensity=scale', 'shelx.composition=C3H7NO2', 'shelx.hklout=alanine.hkl', 'shelx.ins=alanine.ins']),
    ('17_shelxt', [r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxt.exe', 'alanine']),
]
for label, command in commands:
    code = run('work/final_data', label, command)
    if code:
        raise SystemExit(code)
