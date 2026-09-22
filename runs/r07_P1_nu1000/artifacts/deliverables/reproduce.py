"""Reproduce the supplied final framework refinement without changing inputs.
Run with the provided scientific Python. Optional: --shelxl PATH --out NEW_DIR.
"""
from pathlib import Path
import argparse,subprocess,shutil,json,hashlib,sys,re
import gemmi,numpy as np
from annotate_cif import annotate
BASE=Path(__file__).resolve().parent
p=argparse.ArgumentParser(); p.add_argument('--shelxl'); p.add_argument('--out',default=str(BASE/'reproduced'))
a=p.parse_args()
if a.shelxl: exe=a.shelxl
elif (BASE.parent/'environment.json').exists(): exe=json.loads((BASE.parent/'environment.json').read_text())['software']['shelxl']
else: exe='shelxl'
out=Path(a.out).resolve()
if out.exists(): raise RuntimeError('Choose a new output directory; refusing to overwrite '+str(out))
out.mkdir(parents=True)
name='nu1000_framework'
for ext in ['ins','hkl','fab']: shutil.copyfile(BASE/(name+'.'+ext),out/(name+'.'+ext))
command=[exe,name,'-t4']
with (out/'console.log').open('w') as f: result=subprocess.run(command,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=300)
if not (out/(name+'.cif')).exists(): raise RuntimeError('SHELXL did not produce a CIF; see console.log')
annotate(out/(name+'.cif'))
ref=gemmi.cif.read_file(str(BASE/(name+'.cif'))).sole_block(); new=gemmi.cif.read_file(str(out/(name+'.cif'))).sole_block()
keys=['_refine_ls_R_factor_gt','_refine_ls_R_factor_all','_refine_ls_wR_factor_ref','_refine_ls_goodness_of_fit_ref','_refine_ls_number_reflns','_refine_ls_shift/su_max','_refine_diff_density_max','_refine_diff_density_min']
rs=gemmi.make_small_structure_from_block(ref); ns=gemmi.make_small_structure_from_block(new)
rd={s.label:s for s in rs.sites}; changes=[]
for s in ns.sites:
    r=rd[s.label]; delta=np.array([s.fract.x-r.fract.x,s.fract.y-r.fract.y,s.fract.z-r.fract.z]); delta-=np.rint(delta)
    changes.append(float(np.linalg.norm(np.array(rs.cell.orth.mat)@delta)))
record=dict(command=command,exit_code=result.returncode,reference={k:ref.find_value(k) for k in keys},reproduced={k:new.find_value(k) for k in keys},max_site_shift_angstrom=max(changes),hkl_sha256=hashlib.sha256((out/(name+'.hkl')).read_bytes()).hexdigest(),unknown_temperature=new.find_value('_diffrn_ambient_temperature'))
def aniso(block):
    labels=list(block.find_values('_atom_site_aniso_label'))
    cols=[list(block.find_values('_atom_site_aniso_U_'+ij)) for ij in ['11','22','33','23','13','12']]
    return {label:np.array([float(re.sub(r'\([^)]*\)','',c[i])) for c in cols]) for i,label in enumerate(labels)}
ra=aniso(ref); na=aniso(new)
record['same_atom_labels']=set(rd)=={s.label for s in ns.sites}
record['max_adp_component_difference']=max(float(np.max(np.abs(ra[k]-na[k]))) for k in ra)
record['hkl_matches_reference']=record['hkl_sha256']==hashlib.sha256((BASE/(name+'.hkl')).read_bytes()).hexdigest()
record['passed']=(result.returncode==0 and record['same_atom_labels'] and record['hkl_matches_reference']
    and abs(float(ref.find_value(keys[0]))-float(new.find_value(keys[0])))<.0002
    and abs(float(ref.find_value(keys[2]))-float(new.find_value(keys[2])))<.0005
    and max(changes)<.002 and record['max_adp_component_difference']<.001
    and new.find_value('_diffrn_ambient_temperature')=='?' and new.find_value('_cell_measurement_temperature')=='?')
(out/'reproduction_check.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2))
if not record['passed']: sys.exit(1)
