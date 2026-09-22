"""Local scientific/file validation (not a claim of official online checkCIF)."""
import hashlib
import itertools
import json
import re
import shutil
import subprocess
import time
from collections import Counter
import gemmi
import numpy as np
from run_command import ROOT

out=ROOT/'deliverables'
work=ROOT/'work/validation_list4'
work.mkdir(parents=True,exist_ok=True)
report={}
num=gemmi.cif.as_number

def parse_cif(path):
    doc=gemmi.cif.read_file(str(path))
    block=doc.sole_block()
    structure=gemmi.make_small_structure_from_block(block)
    assert len(structure.sites)>0
    return block,structure

b,structure=parse_cif(out/'alanine.cif')
original,_=parse_cif(ROOT/'work/refine/10_publication/alanine.cif')
for tag in ['_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z','_atom_site_U_iso_or_equiv','_atom_site_aniso_U_11','_atom_site_aniso_U_22','_atom_site_aniso_U_33']:
    assert list(b.find_values(tag))==list(original.find_values(tag)),tag
for tag in ['_shelx_res_file','_shelx_hkl_file']:
    assert b.find_value(tag)==original.find_value(tag),tag
report['cif_parses_gemmi']=True
report['model_and_embedded_data_unchanged_by_metadata_enrichment']=True
formula=Counter(s.element.name for s in structure.sites)
assert dict(formula)=={'O':2,'N':1,'H':7,'C':3}
report['asymmetric_unit_formula']=dict(formula)
labels=list(b.find_values('_atom_site_label'))
coords=np.array([[num(b.find_values('_atom_site_fract_'+v)[i]) for v in 'xyz'] for i in range(len(labels))])
ucell=np.array([num(b.find_value('_cell_length_'+v)) for v in 'abc'])
positions=dict(zip(labels,coords))
report['bond_lengths_A']={a+'-'+c:float(np.linalg.norm((positions[a]-positions[c])*ucell)) for a,c in [('O1','C1'),('O2','C1'),('N1','C2'),('C1','C2'),('C2','C3'),('N1','H1A'),('N1','H1B'),('N1','H1C')]}
u={}
for i,label in enumerate(b.find_values('_atom_site_aniso_label')):
    u11,u22,u33,u23,u13,u12=[num(b.find_values('_atom_site_aniso_U_'+v)[i]) for v in ['11','22','33','23','13','12']]
    u[label]=np.array([[u11,u12,u13],[u12,u22,u23],[u13,u23,u33]])
report['nonH_adp_eigenvalues_A2']={k:np.linalg.eigvalsh(v).tolist() for k,v in u.items()}
assert all(min(e)>0 for e in report['nonH_adp_eigenvalues_A2'].values())
uiso=list(b.find_values('_atom_site_U_iso_or_equiv'))
for i,label in enumerate(labels):
    if label not in u:
        u[label]=np.eye(3)*num(uiso[i])
assert all(np.linalg.eigvalsh(v).min()>0 for v in u.values())
report['all_adps_positive']=True
# Independently recompute structure-factor amplitudes from CIF model, IT92
# neutral-atom factors, supplied anomalous corrections, and all symmetry mates.
fdoc=gemmi.cif.read_file(str(out/'alanine.fcf')).sole_block()
hkl=np.array([[int(fdoc.find_values('_refln_index_'+v)[i]) for v in 'hkl'] for i in range(len(fdoc.find_values('_refln_index_h')))])
fo2=np.array([num(x) for x in fdoc.find_values('_refln_F_squared_meas')])
sig=np.array([num(x) for x in fdoc.find_values('_refln_F_squared_sigma')])
assert int(fdoc.find_value('_shelx_refln_list_code'))==4
fc2=np.array([num(x) for x in fdoc.find_values('_refln_F_squared_calc')])
fc=np.sqrt(fc2)
q=hkl/ucell
stol2=(q*q).sum(axis=1)/4
factors={}
for i,e in enumerate(b.find_values('_atom_type_symbol')):
    e=gemmi.cif.as_string(e)
    fp=num(b.find_values('_atom_type_scat_dispersion_real')[i]); fdp=num(b.find_values('_atom_type_scat_dispersion_imag')[i])
    factors[e]=np.array([gemmi.Element(e).it92.calculate_sf(x) for x in stol2])+fp+1j*fdp
ops=[gemmi.Op(gemmi.cif.as_string(x)) for x in b.find_values('_space_group_symop_operation_xyz')]
assert len(ops)==4
calculated=np.zeros(len(hkl),dtype=complex)
for i,label in enumerate(labels):
    element=gemmi.cif.as_string(b.find_values('_atom_site_type_symbol')[i])
    occupancy=num(b.find_values('_atom_site_occupancy')[i])
    for op in ops:
        r=np.array(op.rot,dtype=float)/op.DEN
        t=np.array(op.tran,dtype=float)/op.DEN
        pos=r@coords[i]+t
        umat=r@u[label]@r.T
        damp=np.exp(-2*np.pi*np.pi*np.einsum('ij,jk,ik->i',q,umat,q))
        calculated+=occupancy*factors[element]*damp*np.exp(2j*np.pi*(hkl@pos))
fc_independent=abs(calculated)
rel=float(abs(fc-fc_independent).sum()/fc.sum())
report['independent_Fcalc_vs_FCF_amplitude_R']=rel
report['independent_Fcalc_note']='CIF-rounded coordinates/ADPs, IT92 scattering factors, and supplied anomalous terms; a small nonzero difference from SHELXL full precision is expected.'
assert rel<0.005,rel
fobs=np.sqrt(np.maximum(fo2,0)); gt=fo2>2*sig
r1all=float(abs(fobs-fc).sum()/fobs.sum()); r1gt=float(abs(fobs[gt]-fc[gt]).sum()/fobs[gt].sum())
pvalues=(np.maximum(fo2,0)+2*fc2)/3
weights=1/(sig**2+(0.0375*pvalues)**2+0.0505*pvalues)
wr2=float(np.sqrt(np.sum(weights*(fo2-fc2)**2)/np.sum(weights*fo2**2)))
gof=float(np.sqrt(np.sum(weights*(fo2-fc2)**2)/(len(hkl)-68)))
report['FCF']={'list_code':4,'reflections':len(hkl),'observed':int(gt.sum()),'R1_all_recomputed':r1all,'R1_gt_recomputed':r1gt,'wR2_recomputed':wr2,'GoF_recomputed':gof,'refinement_reflections':int(b.find_value('_refine_ls_number_reflns')),'note':'LIST 4 retains Friedel-distinct observations and is compared directly with the refinement statistics.'}
assert len(hkl)==int(b.find_value('_refine_ls_number_reflns')),report['FCF']
assert abs(r1all-num(b.find_value('_refine_ls_R_factor_all')))<0.0002,report['FCF']
assert abs(r1gt-num(b.find_value('_refine_ls_R_factor_gt')))<0.0002,report['FCF']
assert abs(wr2-num(b.find_value('_refine_ls_wR_factor_ref')))<0.001,report['FCF']
assert abs(gof-num(b.find_value('_refine_ls_goodness_of_fit_ref')))<0.01,report['FCF']
# Resolution sensitivity: compare full-precision RES coordinates, not rounded CIF.
def res_coordinates(path):
    result={}
    for line in path.read_text().split('HKLF')[0].splitlines():
        fields=line.split()
        if len(fields)>=7 and fields[0] in labels:
            result[fields[0]]=np.array([float(x) for x in fields[2:5]])
    return result
refxyz=res_coordinates(ROOT/'work/refine/07_final/alanine.res')
report['refinement_controls']={}
for stage in ['02b_riding_H','03_masked_riding','04_NH_omit','05_NH_free','06_weighted','07_final','08_cutoff_080','09_full_measured']:
    path=ROOT/'work/refine'/stage/'alanine.res'
    text=path.read_text()
    metrics=[line[4:] for line in text.splitlines() if line.startswith(('REM wR2','REM R1','REM Highest'))]
    xyz=res_coordinates(path)
    shift={label:float(np.linalg.norm((xyz[label]-refxyz[label])*ucell)) for label in xyz}
    report['refinement_controls'][stage]={'summary':metrics,'max_nonH_shift_A_vs_final':max(shift[x] for x in shift if not x.startswith('H')),'max_H_shift_A_vs_final':max((shift[x] for x in shift if x.startswith('H')),default=None)}
# Read-only-input check against import-time raw-frame hashes.
manifest=json.loads((ROOT/'raw_frames_manifest.json').read_text())
failed=[]
for entry in manifest:
    path=ROOT/entry['path']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256'] or path.stat().st_size!=entry['bytes']:
        failed.append(entry['path'])
assert not failed,failed
report['raw_input_hashes']={'checked':len(manifest),'changed':failed}
# A clean-room final refinement rerun, isolated from the delivered files.
rerun=work/'shelxl_rerun'
rerun.mkdir(exist_ok=False)
for ext in ['ins','hkl']:
    shutil.copy2(out/f'alanine.{ext}',rerun/f'alanine.{ext}')
exe=r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe'
t0=time.time()
p=subprocess.run([exe,'alanine'],cwd=rerun,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
(rerun/'console.log').write_bytes(p.stdout)
br,_=parse_cif(rerun/'alanine.cif')
metric_tags=['_refine_ls_R_factor_gt','_refine_ls_R_factor_all','_refine_ls_wR_factor_ref','_refine_ls_goodness_of_fit_ref','_refine_ls_number_parameters','_refine_ls_number_reflns']
report['clean_room_refinement_rerun']={'argv':[exe,'alanine'],'cwd':str(rerun),'exit_code':p.returncode,'seconds':time.time()-t0,'metrics':{k:br.find_value(k) for k in metric_tags}}
assert p.returncode==0
assert all(br.find_value(k)==b.find_value(k) for k in metric_tags)
# Try the supplied local PLATON validator. No online service is contacted here.
platon=work/'platon'
platon.mkdir(exist_ok=False)
for ext in ['cif','fcf']:
    shutil.copy2(out/f'alanine.{ext}',platon/f'alanine.{ext}')
argv=[r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\platon.exe','-u','alanine.cif']
t0=time.time()
try:
    p=subprocess.run(argv,cwd=platon,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=45)
    (platon/'console.log').write_bytes(p.stdout)
    pstatus={'exit_code':p.returncode,'output':p.stdout.decode(errors='replace')[-3000:]}
except subprocess.TimeoutExpired as exc:
    (platon/'console.log').write_bytes(exc.stdout or b'')
    pstatus={'timeout_seconds':45,'status':'No completed local validation report; not a pass'}
files=[x.name for x in platon.iterdir() if x.name not in ['alanine.cif','alanine.fcf','console.log']]
report['platon']={'argv':argv,'cwd':str(platon),'seconds':time.time()-t0,'generated_files':files,**pstatus}
if (platon/'alanine.chk').exists():
    check=(platon/'alanine.chk').read_text(errors='replace')
    report['platon']['written_report_has_summary']='ALERT_Level and ALERT_Type Summary' in check
    report['platon']['alerts']=re.findall(r'^\s*\d{3}_ALERT_\d_[ABCG].*$',check,re.M)
    report['platon']['alert_counts']=dict(Counter(re.findall(r'^\s*\d{3}_ALERT_\d_([ABCG])',check,re.M)))
    report['platon']['status']='Written local report inspected; process may require timeout termination. Outstanding alerts are not a validation pass.'
(work/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(out/'validation').mkdir(exist_ok=True)
shutil.copy2(work/'validation.json',out/'validation/validation.json')
print(json.dumps(report,indent=2))
