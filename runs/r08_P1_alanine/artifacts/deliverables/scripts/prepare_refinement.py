"""Prepare the chemically assigned non-hydrogen model from SHELXT coordinates."""
from pathlib import Path
import shutil
import numpy as np
from dxtbx.model.experiment_list import ExperimentListFactory
from run_command import ROOT, run

crystal = ExperimentListFactory.from_json_file(str(ROOT / 'work/dials/refined_joint_esd.expt'), check_format=False)[0].crystal
cell = crystal.get_unit_cell().parameters()
esds = crystal.get_cell_parameter_sd()
if len(esds) != 6 or not all(x > 0 for x in esds[:3]):
    raise RuntimeError('Crystallographic cell uncertainties not recovered')
print('Joint-refinement cell:', cell, 'esds:', esds, 'volume:', crystal.get_unit_cell().volume(), 'volume esd:', crystal.get_cell_volume_sd())
old = ExperimentListFactory.from_json_file(str(ROOT / 'work/dials/refined_joint.expt'), check_format=False)[0].crystal
assert np.allclose(cell, old.get_unit_cell().parameters(), atol=1e-10)
source = ROOT / 'work/final_data/alanine_a.res'
mapping = {'O001': ('O1',4), 'O002': ('O2',4), 'O003': ('N1',3), 'C004': ('C1',1), 'C005': ('C2',1), 'C006': ('C3',1)}
atoms = []
for line in source.read_text().splitlines():
    fields = line.split()
    if fields and fields[0] in mapping:
        name, sfac = mapping[fields[0]]
        xyz = tuple(float(x) for x in fields[2:5])
        atoms.append((name, sfac, xyz, float(fields[6])))
print('Distances below 1.7 A for proposed non-H model:')
for i,a in enumerate(atoms):
    for b in atoms[:i]:
        d = crystal.get_unit_cell().distance(a[2],b[2])
        if d < 1.7:
            print(a[0], b[0], round(d,4))
header = 'TITL alanine raw-data solution, non-H diagnostic refinement\n'
header += 'CELL 0.71073 ' + ' '.join(f'{x:.6f}' for x in cell[:3]) + ' 90 90 90\n'
header += 'ZERR 4 ' + ' '.join(f'{x:.8f}' for x in esds[:3]) + ' 0 0 0\n'
header += '''LATT -1
SYMM 1/2-X,-Y,1/2+Z
SYMM -X,1/2+Y,1/2-Z
SYMM 1/2+X,1/2-Y,-Z
SFAC C H N O
UNIT 12 28 4 8
TEMP -123.15
SHEL 999 0.77
L.S. 20
ANIS
WGHT 0.05
FVAR 1.0
BOND $H
FMAP 2
PLAN 20
ACTA
LIST 6
'''
body = ''.join(f'{name:4s} {sfac} ' + ' '.join(f'{x:.6f}' for x in xyz) + f' 11.00000 {u:.6f}\n' for name,sfac,xyz,u in atoms)
work = ROOT / 'work/refine/01b_nonH'
work.mkdir(parents=True, exist_ok=True)
(work / 'alanine.ins').write_text(header + body + 'HKLF 4\nEND\n', encoding='ascii')
shutil.copy2(ROOT / 'work/final_data/alanine.hkl', work / 'alanine.hkl')
raise SystemExit(run('work/refine/01b_nonH','19_shelxl_nonH',[r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe','alanine']))
