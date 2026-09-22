"""Self-written numerical, chemical, symmetry and CIF-readback checks.
This is NOT an independent standard checkCIF/PLATON validation.
"""
import json,re,shlex,itertools
from collections import Counter
import numpy as np
from raw_frames import OUT
from solve import SIGNS,TRANS
from refine import StructureFactors,fit_stage,unpack,rvalues

def parse_cif(path):
 lines=path.read_text(encoding='utf8').splitlines();tokens=[];i=0
 while i<len(lines):
  if lines[i].startswith(';'):
   block=[lines[i][1:]];i+=1
   while i<len(lines) and not lines[i].startswith(';'):block.append(lines[i]);i+=1
   if i==len(lines):raise ValueError('Unterminated CIF text')
   tokens.append('\n'.join(block))
  else:tokens.extend(shlex.split(lines[i],comments=True,posix=True))
  i+=1
 tags={};loops=[];i=0
 while i<len(tokens):
  t=tokens[i]
  if t.lower().startswith('data_'):i+=1;continue
  if t=='loop_':
   i+=1;names=[]
   while i<len(tokens) and tokens[i].startswith('_'):names.append(tokens[i]);i+=1
   vals=[]
   while i<len(tokens) and not (tokens[i].startswith('_') or tokens[i] in ['loop_','stop_'] or tokens[i].lower().startswith('data_')):vals.append(tokens[i]);i+=1
   if not names or len(vals)%len(names):raise ValueError('Incomplete CIF loop')
   loops.append((names,[vals[j:j+len(names)] for j in range(0,len(vals),len(names))]))
  elif t.startswith('_'):
   if t in tags:raise ValueError('Duplicate CIF tag '+t)
   tags[t]=tokens[i+1];i+=2
  else:raise ValueError('Unexpected CIF token '+t)
 return tags,loops

def number(t):return float(re.sub(r'\(\d+\)$','',t))
def find_loop(loops,key):return next((n,r) for n,r in loops if key in n)
def angle(u,v):return float(np.rad2deg(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))

def main():
 m=json.loads((OUT/'refined_model.json').read_text());cell=np.array(m['cell']);p=np.array(m['parameters']);x=np.array(m['coordinates']);u=np.array(m['U_cartesian'])
 r=np.loadtxt(OUT/'refinement_reflections.csv',delimiter=',',skiprows=1);h=r[:,:3];I=r[:,3];sig=r[:,4]
 sf=StructureFactors(h,cell,True);fc=sf.calculate(x,u);calc=m['scale']*abs(fc)**2
 checks={};checks['stored_calculated_intensities_match']=bool(np.allclose(calc,r[:,5],rtol=2e-9,atol=2e-6))
 checks['formula_and_H_count']=Counter(m['types'])==Counter({'C':3,'N':1,'O':2,'H':7})
 checks['all_ADP_positive_definite']=bool(np.min(np.linalg.eigvalsh(u))>0)
 checks['no_extreme_nonH_ADP']=bool(np.max(np.linalg.eigvalsh(u[:6]))<.1)
 checks['all_values_finite']=bool(np.isfinite(x).all() and np.isfinite(u).all() and np.isfinite(r).all())
 closure=[]
 for a,b in itertools.product(range(4),repeat=2):
  s=SIGNS[a]*SIGNS[b];t=(SIGNS[a]*TRANS[b]+TRANS[a])%1
  closure.append(any(np.array_equal(s,c) and np.allclose((t-d+.5)%1-.5,0) for c,d in zip(SIGNS,TRANS)))
 checks['four_symmetry_operations_close']=all(closure)
 absent=np.array([np.eye(3,dtype=int)[i]*n for i in range(3) for n in [1,3,5,7]])
 f_abs=StructureFactors(absent,cell,True).calculate(x,u)
 checks['axial_odd_reflections_cancel']=bool(np.max(abs(f_abs))<1e-8)
 checks['Friedel_intensities_equal']=bool(np.allclose(abs(fc)**2,abs(StructureFactors(-h,cell,True).calculate(x,u))**2))
 # CIF readback is intentionally independent of the writer, but remains task-local.
 tags,loops=parse_cif(OUT/'alanine.cif');names,rows=find_loop(loops,'_atom_site_label')
 cx=np.array([[number(row[names.index('_atom_site_fract_'+a)]) for a in 'xyz'] for row in rows])
 checks['CIF_has_13_atoms']=len(rows)==13
 ctol=np.array([[.5001*10**(-len(re.sub(r'\(\d+\)$','',row[names.index('_atom_site_fract_'+a)]).split('.')[-1])) for a in 'xyz'] for row in rows])
 checks['CIF_coordinates_roundtrip']=bool(np.all(abs(cx-x)<=ctol))
 checks['CIF_R1_matches']=abs(number(tags['_refine_ls_R_factor_gt'])-m['statistics']['R1_gt2sigma'])<5.1e-7
 _,fl=parse_cif(OUT/'alanine.fcf');fn,fr=find_loop(fl,'_refln_index_h');fr=np.array(fr,float)
 checks['FCF_reflection_count']=len(fr)==len(r)
 checks['FCF_intensities_roundtrip']=bool(np.allclose(fr[:,3],I/m['scale'],atol=1e-7,rtol=1e-8))
 checks['FCF_calculated_intensities_roundtrip']=bool(np.allclose(fr[:,5],abs(fc)**2,atol=1e-6,rtol=1e-8))
 # Bond angles and bonded-atom rigid-bond displacement differences.
 cart=x*cell;bonds=[(0,1),(0,2),(0,3),(2,4),(2,5)];angles=[];rigid=[]
 for a,b,c in [(1,0,2),(1,0,3),(2,0,3),(0,2,4),(0,2,5),(4,2,5)]:
  angles.append(dict(atoms=[m['labels'][i] for i in [a,b,c]],degrees=angle(cart[a]-cart[b],cart[c]-cart[b])))
 for a,b in bonds:
  v=cart[b]-cart[a];v/=np.linalg.norm(v)
  rigid.append(dict(atoms=[m['labels'][a],m['labels'][b]],delta_U_parallel_A2=float(v@(u[a]-u[b])@v)))
 contacts=[];hbonds=[]
 for op,shift in itertools.product(range(4),itertools.product([-1,0,1],repeat=3)):
  pos=(x*SIGNS[op]+TRANS[op]+np.array(shift))*cell
  for a,b in itertools.product(range(6),repeat=2):
   if op==0 and shift==(0,0,0):continue
   dist=float(np.linalg.norm(cart[a]-pos[b]))
   if dist<3.5:contacts.append(dict(atom1=m['labels'][a],atom2=m['labels'][b],operation=op+1,cell_shift=list(shift),distance=dist))
  for ih,io in itertools.product([7,8,9],[4,5]):
   dha=angle(cart[1]-cart[ih],pos[io]-cart[ih]);ha=float(np.linalg.norm(pos[io]-cart[ih]));da=float(np.linalg.norm(pos[io]-cart[1]))
   if ha<2.7 and da<3.6 and dha>110:hbonds.append(dict(donor='N1',H=m['labels'][ih],acceptor=m['labels'][io],operation=op+1,cell_shift=list(shift),D_H=float(np.linalg.norm(cart[1]-cart[ih])),H_A=ha,D_A=da,D_H_A=dha))
 contacts.sort(key=lambda z:z['distance']);hbonds.sort(key=lambda z:z['H_A'])
 checks['no_inter_molecular_nonH_contact_below_2A']=not contacts or contacts[0]['distance']>2
 # An additional full-data step quantifies actual parameter shift, not merely termination status.
 fit,extra,_,_=fit_stage(p,'aniso',h,I,sig,cell,True,max_nfev=35)
 esd=np.sqrt(np.diag(np.load(OUT/'refinement_covariance.npy')))
 shifts=abs(fit.x-p)/np.maximum(esd,1e-10)
 checks['extra_cycle_max_shift_below_0p01_su']=bool(shifts.max()<.01)
 # Conditional five-fold predictive check. Starting orientation/hand comes from
 # the solved structure; every non-H xyz/U and other parameter is refit to training I.
 rng=np.random.default_rng(220927);fold=np.empty(len(h),int);fold[rng.permutation(len(h))]=np.arange(len(h))%5
 pred=np.zeros(len(h));cv=[]
 for k in range(5):
  train=fold!=k;test=~train
  f,st,_,_=fit_stage(p,'aniso',h[train],I[train],sig[train],cell,True,max_nfev=80)
  xx,uu,sc,_=unpack(f.x,'aniso',cell,True);ft=StructureFactors(h[test],cell,True).calculate(xx,uu)
  pred[test]=sc*abs(ft)**2
  cv.append(dict(fold=k,n_test=int(test.sum()),test_statistics=rvalues(I[test],sig[test],ft,sc,0),train_statistics=st))
 weights=r[:,8];sel=I>2*sig;fo=np.sqrt(np.maximum(I,0));fcv=np.sqrt(pred)
 cross=dict(method='Five-fold conditional parameter-refinement prediction; solved orientation/hand shared, not a blind solution test',seed=220927,folds=cv,R1_gt2sigma=float(abs(fo[sel]-fcv[sel]).sum()/fo[sel].sum()),R1_all=float(abs(fo-fcv).sum()/fo.sum()),wR2=float(np.sqrt(np.sum(weights*(I-pred)**2)/np.sum(weights*I*I))))
 np.savetxt(OUT/'crossvalidation_reflections.csv',np.c_[h,fold,I,sig,pred],delimiter=',',header='h,k,l,heldout_fold,Fo2,sigma,heldout_predicted_Fc2',comments='',fmt='%.10g')
 merged=np.load(OUT/'merged.npy');ab=merged[merged[:,7]>0]
 st=m['statistics'];diff=json.loads((OUT/'difference_statistics.json').read_text())
 quality=dict(R1_observed_below_0p05=st['R1_gt2sigma']<.05,wR2_below_0p12=st['wR2']<.12,GoF_below_2=st['goodness_of_fit']<2,difference_extrema_below_0p5=max(abs(diff['min_e_A3']),abs(diff['max_e_A3']))<.5)
 result=dict(quality_targets=quality,quality_targets_met=all(quality.values()),checks=checks,all_basic_checks_pass=all(checks.values()),independent_standard_validation=False,ADP_eigenvalues=np.linalg.eigvalsh(u[:6]).tolist(),Ueq=(np.trace(u[:6],axis1=1,axis2=2)/3).tolist(),angles=angles,rigid_bond_differences=rigid,intermolecular_contacts=contacts,hydrogen_bonds=hbonds,max_additional_parameter_shift_over_su=float(shifts.max()),cross_validation=cross,axial_absence_max_abs_I_over_sigma=float(np.max(abs(ab[:,3]/ab[:,4]))),axial_absence_n=len(ab),warnings=['Absolute configuration undetermined; Friedel pairs merged','No independent standard crystallographic validation','No empirical/physical absorption correction beyond run/module relative scale','Hydrogen positions ride; H displacement multipliers imposed','Geometry and intensity errors include unmodelled systematic components','High-resolution completeness is incomplete at 0.70 A'])
 (OUT/'validation.json').write_text(json.dumps(result,indent=2));print('BASIC CHECKS',json.dumps(checks),flush=True);print('Cross-validated R1',cross['R1_gt2sigma'],'extra-cycle max shift/su',shifts.max(),flush=True)
 if not all(checks.values()):raise RuntimeError('Basic validation failed; inspect validation.json')
 if not all(quality.values()):raise RuntimeError('Numerical quality targets not met; preserve outputs and inspect validation.json')
if __name__=='__main__':main()
