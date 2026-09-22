"""Conditional geometry covariance and held-out peak prediction test."""
import json
import numpy as np
from scipy.optimize import least_squares
from raw_frames import OUT
from geometry import pixel_rays,ub_runs
from geometry_probe import rotate_z

SCALE=np.r_[[.1]*3,[.01]*18,[1,5,5,.2,.2,.2]*2]

def residual(p,obs):
 q=rotate_z(pixel_rays(p,obs[:,4],obs[:,5],obs[:,7])-[1,0,0],obs[:,6])
 pred=np.einsum('nij,nj->ni',ub_runs(p)[obs[:,0].astype(int)-1],obs[:,1:4])
 return (q-pred).ravel()

def main():
 p=np.array(json.loads((OUT/'geometry.json').read_text())['parameters'])
 obs=np.loadtxt(OUT/'calibration_centroids.csv',delimiter=',',skiprows=1)
 err=np.linalg.norm(residual(p,obs).reshape(-1,3),axis=1)
 use=err<.006;obs=obs[use]
 r=residual(p,obs);j=np.zeros((len(r),len(p)))
 for i,s in enumerate(SCALE):
  pp=p.copy();pm=p.copy();pp[i]+=s*1e-4;pm[i]-=s*1e-4
  j[:,i]=(residual(pp,obs)-residual(pm,obs))/2e-4
 _,sv,vt=np.linalg.svd(j,full_matrices=False)
 rank=int(np.sum(sv>sv[0]*1e-8));covq=np.linalg.pinv(j.T@j,rcond=1e-16)*(r@r/(len(r)-rank));cov=covq*SCALE[:,None]*SCALE[None,:]
 np.save(OUT/'geometry_covariance.npy',cov)
 rng=np.random.default_rng(220926);test=rng.random(len(obs))<.2
 f=least_squares(lambda q:residual(p+q*SCALE,obs[~test])/.001,np.zeros(len(p)),loss='soft_l1',f_scale=.6,max_nfev=100)
 pt=p+f.x*SCALE;held=residual(pt,obs[test]).reshape(-1,3);train=residual(pt,obs[~test]).reshape(-1,3)
 # Peak-prediction discrepancy in reciprocal units is not an angular pixel RMS.
 result=dict(n_parameters=len(p),jacobian_rank=rank,scaled_singular_values=sv.tolist(),cell_standard_uncertainty=np.sqrt(np.diag(cov)[:3]).tolist(),conditional_uncertainty=True,uncertainty_note='Linearized covariance with residual variance; ignores external wavelength, detector-geometry systematics and correlated centroid errors',n_calibration=len(obs),n_holdout=int(test.sum()),holdout_reciprocal_rms=float(np.sqrt(np.mean(np.sum(held*held,axis=1)))),train_reciprocal_rms=float(np.sqrt(np.mean(np.sum(train*train,axis=1)))),split_fit_cell=pt[:3].tolist(),full_cell=p[:3].tolist(),split_cell_shift=(pt[:3]-p[:3]).tolist(),holdout_seed=220926)
 (OUT/'geometry_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
