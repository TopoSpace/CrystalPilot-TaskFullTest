from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from geometry import *

OUT=Path(__file__).resolve().parent
p=np.load(OUT/"peaks.npy")
initial=np.array([47,406.52,380.45,.03374,49.92498,.29559,-.36049])
hf=calibrated_q(p,initial)@np.linalg.inv(UB0).T
h=np.rint(hf).astype(int)
e=np.linalg.norm(hf-h,axis=1)
print("initial indexing",np.quantile(e,[0,.5,.9,.95,.99]),np.mean(e<.12))
rows=[]
for key,indices in pd.DataFrame(np.c_[p[:,0],h]).groupby([0,1,2,3]).groups.items():
    ii=np.array(list(indices))
    ii=ii[e[ii]<.14]
    if len(ii)==0:
        continue
    pp=p[ii]
    if pp[:,4].sum()<150 or np.ptp(pp[:,6])>3:
        continue
    avg=np.average(pp,axis=0,weights=pp[:,4])
    avg[4]=pp[:,4].sum()
    rows.append(np.r_[avg,key[1:]])
g=np.array(rows)
print("centroids",len(g))
gh=g[:,10:13]
weights=np.minimum(np.sqrt(g[:,4]/300),3)


def residual(z):
    q=calibrated_q(g,z[9:])
    return ((q-gh@z[:9].reshape(3,3).T)*weights[:,None]).ravel()


z0=np.r_[UB0.ravel(),initial]
lo=np.r_[np.full(9,-.3),40,390,365,-2,47,-2,-.360491]
hi=np.r_[np.full(9,.3),55,425,395,2,53,2,-.360489]
fit=least_squares(residual,z0,bounds=(lo,hi),loss="soft_l1",f_scale=.0004,
                  x_scale="jac",max_nfev=150,verbose=1,ftol=1e-11,xtol=1e-11,gtol=1e-11)
ub=fit.x[:9].reshape(3,3)
pr=fit.x[9:]
err=calibrated_q(g,pr)-gh@ub.T
direct=WAVELENGTH*np.linalg.inv(ub).T
cell=np.linalg.norm(direct,axis=0)
angles=np.degrees(np.arccos(np.clip((direct.T@direct)/np.outer(cell,cell),-1,1)))
print("UB",ub,"params",pr,"cell",cell,"angles",angles,"rms",np.sqrt(np.mean(err**2)),sep="\n")
result=dict(ub=ub.tolist(),params=pr.tolist(),cell=cell.tolist(),angles=angles.tolist(),
            residual_q_rms=float(np.sqrt(np.mean(err**2))),centroids=len(g),
            parameter_order=["distance","center_x","center_y","theta_zero","kappa_alpha","kappa_zero","omega_zero"])
(OUT/"geometry.json").write_text(json.dumps(result,indent=2))
np.savetxt(OUT/"indexed_centroids.csv",np.c_[g,err],delimiter=",",
           header="run,frame,x,y,sum,maximum,omega,theta,kappa,phi,h,k,l,dqx,dqy,dqz",comments="")
