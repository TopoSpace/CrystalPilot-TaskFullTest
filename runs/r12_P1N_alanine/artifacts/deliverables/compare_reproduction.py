"""Compare a fresh raw-data solve with the delivered model up to symmetry/hand."""
import json,itertools
import numpy as np
from raw_frames import OUT
from solve import SIGNS,TRANS
from refine import StructureFactors

def main():
 newdir=OUT/'reproduction_verified'
 a=json.loads((OUT/'refined_model.json').read_text());b=json.loads((newdir/'refined_model.json').read_text())
 cell=np.array(a['cell']);xa=np.array(a['coordinates'])[:6];xb=np.array(b['coordinates'])[:6]
 best=1e9;details=None
 for hand,op,shift,swap in itertools.product([-1,1],range(4),itertools.product([0,.5],repeat=3),[False,True]):
  xp=hand*(xb*SIGNS[op]+TRANS[op])+np.array(shift)
  if swap:xp=xp[[0,1,2,3,5,4]]
  dv=((xp-xa+.5)%1-.5)*cell;cost=float(np.sqrt(np.mean(np.sum(dv*dv,axis=1))))
  if cost<best:best=cost;details=dict(hand=hand,symmetry_operation=op+1,origin_shift=list(shift),oxygen_exchange=swap,per_atom_displacements_A=np.linalg.norm(dv,axis=1).tolist())
 ma=np.load(OUT/'merged.npy');mb=np.load(newdir/'merged.npy');same=bool(ma.shape==mb.shape and np.allclose(ma,mb,rtol=1e-10,atol=1e-10,equal_nan=True))
 h=ma[ma[:,7]==0,:3];sf=StructureFactors(h,cell,True)
 fa=abs(sf.calculate(np.array(a['coordinates']),np.array(a['U_cartesian'])))
 fb=abs(sf.calculate(np.array(b['coordinates']),np.array(b['U_cartesian'])))
 amp_r=float(abs(fa-fb).sum()/fa.sum())
 result=dict(merged_data_reproduced=same,nonH_RMS_A_after_symmetry_and_hand=best,matching_transform=details,calculated_amplitude_relative_discrepancy=amp_r,main_statistics=a['statistics'],fresh_statistics=b['statistics'],fresh_validation=json.loads((newdir/'validation.json').read_text())['all_basic_checks_pass'])
 result['passed']=bool(same and best<.002 and amp_r<.002 and result['fresh_validation'])
 (OUT/'reproducibility_comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
 if not result['passed']:raise RuntimeError('Fresh reproduction does not yet agree; inspect results')
if __name__=='__main__':main()
