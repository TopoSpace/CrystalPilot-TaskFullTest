"""Independent diffraction geometry calibration against observed spot centroids."""
import json,time
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from raw_frames import OUT
from geometry_probe import rays,rz,rotate_z
LAMBDA=.71073
CELL0=np.array([5.7936114,5.9671562,12.2910746])
UB0=np.array([[-.0036164,.1093447,-.0225878],[.0317839,-.0440864,-.051439],[-.1181851,-.0152375,-.0131508]])

def load_data():
 p=np.load(OUT/'peaks.npy');meta=json.loads((OUT/'frame_metadata.json').read_text())
 ix=p[:,2].astype(int)
 angles=np.array([[(m['omega0']+m['omega1'])/2,m['theta'],m['kappa'],m['phi']] for m in meta])
 return p,meta,angles[ix]

def initial_ub_runs(meta):
 u,ss,vt=np.linalg.svd(UB0@np.diag(CELL0/LAMBDA));u=u@vt
 axis=np.array([-np.sin(np.deg2rad(49.92498)),0,np.cos(np.deg2rad(49.92498))])
 mats=[]
 for run in range(1,7):
  m=next(m for m in meta if int(m['file'].split('_')[-2])==run)
  rk=Rotation.from_rotvec(-np.deg2rad(m['kappa'])*axis).as_matrix()
  mats.append(rk@rz(-m['phi'])@u)
 return np.array(mats)

def initial_params(meta):
 rot=Rotation.from_matrix(initial_ub_runs(meta)).as_rotvec().ravel()
 return np.r_[CELL0,rot,[47.901,241.275,380.45,0,0,0],[47.324,569.88,380.45,0,0,0]]

def detector_basis(par):
 mod=par[21:].reshape(2,6)
 normals=[];exs=[];eys=[]
 for j in range(2):
  a=np.deg2rad([-19,19.016][j]);r=Rotation.from_rotvec(np.deg2rad(mod[j,3:6])).as_matrix()
  normals.append(r@np.array([np.cos(a),-np.sin(a),0]));exs.append(r@np.array([-np.sin(a),-np.cos(a),0]));eys.append(r@np.array([0,0,-1]))
 return mod,np.array(normals),np.array(exs),np.array(eys)

def pixel_rays(par,x,y,theta):
 mod,n,ex,ey=detector_basis(par);j=(x>=400).astype(int)
 v=mod[j,0,None]*n[j]+((x-mod[j,1])*.1)[:,None]*ex[j]+((y-mod[j,2])*.1)[:,None]*ey[j]
 v=rotate_z(v,-theta)
 return v/np.linalg.norm(v,axis=1)[:,None]

def ub_runs(par):
 return Rotation.from_rotvec(par[3:21].reshape(6,3)).as_matrix()@np.diag(LAMBDA/par[:3])

def indices(par,p,ang):
 q=pixel_rays(par,p[:,3],p[:,4],ang[:,1])-[1,0,0]
 q=rotate_z(q,ang[:,0])
 ub=ub_runs(par)
 return np.einsum('nij,nj->ni',np.linalg.inv(ub)[p[:,0].astype(int)-1],q)

def group_peaks(p,ang,h):
 groups={}
 for i in range(len(p)):
  key=(int(p[i,0]),*map(int,h[i]));frame=int(p[i,1]);g=groups.setdefault(key,{})
  if frame not in g or p[i,5]>p[g[frame],5]:g[frame]=i
 rows=[]
 for key,g in groups.items():
  ids=np.array(list(g.values()));wt=p[ids,5];s=wt.sum()
  if s<180:continue
  om=ang[ids,0]
  if np.ptp(om)>4:continue
  mean=lambda z:float(np.dot(wt,z)/s)
  rows.append([*key,mean(p[ids,3]),mean(p[ids,4]),mean(om),ang[ids[0],1],s,len(ids),np.sqrt(mean((om-mean(om))**2))])
 return np.array(rows)

def calibrate():
 p,meta,ang=load_data();par=initial_params(meta);hf=indices(par,p,ang);h=np.rint(hf)
 keep=(np.linalg.norm(h-hf,axis=1)<.18)&(p[:,5]>40)
 obs=group_peaks(p[keep],ang[keep],h[keep]);print('calibration groups',len(obs),flush=True)
 np.savetxt(OUT/'calibration_centroids.csv',obs,delimiter=',',header='run,h,k,l,x,y,omega,theta,counts,nframes,omega_width',comments='')
 def residual(par,obs):
  q=pixel_rays(par,obs[:,4],obs[:,5],obs[:,7])-[1,0,0]
  q=rotate_z(q,obs[:,6]);ub=ub_runs(par)
  pred=np.einsum('nij,nj->ni',ub[obs[:,0].astype(int)-1],obs[:,1:4])
  return (q-pred).ravel()/.001
 xscale=np.r_[[.1]*3,[.01]*18,[1,5,5,.2,.2,.2]*2]
 lb=par-np.r_[[.3]*3,[.15]*18,[4,30,30,3,3,3]*2]
 ub=par+np.r_[[.3]*3,[.15]*18,[4,30,30,3,3,3]*2]
 res=least_squares(residual,par,args=(obs,),loss='soft_l1',f_scale=1,x_scale=xscale,bounds=(lb,ub),max_nfev=120)
 err=np.linalg.norm(residual(res.x,obs).reshape(-1,3),axis=1)*.001
 sel=err<.006
 res=least_squares(residual,res.x,args=(obs[sel],),loss='soft_l1',f_scale=.6,x_scale=xscale,bounds=(lb,ub),max_nfev=100)
 err=np.linalg.norm(residual(res.x,obs).reshape(-1,3),axis=1)*.001
 hf=indices(res.x,p,ang);derr=np.linalg.norm(hf-np.rint(hf),axis=1)
 print('cell',res.x[:3],'detector',res.x[21:].reshape(2,6),flush=True)
 print('centroid reciprocal error quantiles',np.percentile(err,[0,50,90,95,99,100]),flush=True)
 print('indexed fraction all/strong',np.mean(derr<.15),np.mean(derr[p[:,5]>100]<.15),flush=True)
 result=dict(parameters=res.x.tolist(),cell=res.x[:3].tolist(),wavelength=LAMBDA,ub_runs=ub_runs(res.x).tolist(),n_calibration=len(obs),n_retained=int(sel.sum()),reciprocal_error_quantiles=np.percentile(err,[0,50,90,95,99,100]).tolist(),indexed_fraction_all=float(np.mean(derr<.15)),indexed_fraction_strong=float(np.mean(derr[p[:,5]>100]<.15)),optimizer_message=res.message)
 (OUT/'geometry.json').write_text(json.dumps(result,indent=2))
 np.save(OUT/'indexed_peaks.npy',np.c_[p,ang,hf,derr])
 np.save(OUT/'calibration_residuals.npy',np.c_[obs,err])
if __name__=='__main__':calibrate()
