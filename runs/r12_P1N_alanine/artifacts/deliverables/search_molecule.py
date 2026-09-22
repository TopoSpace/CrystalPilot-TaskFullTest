"""Direct-space global search from a GENERIC ideal alanine connectivity.
This is a search hypothesis built from ideal covalent distances and tetrahedral/
trigonal angles, NOT a supplied or looked-up crystal/molecular structure.
All non-H coordinates are released in subsequent reciprocal-space refinement.
The template hand is arbitrary; diffraction here is not used to assign D/L.
"""
import json,time
import numpy as np
from scipy.optimize import differential_evolution,least_squares
from scipy.spatial.transform import Rotation
from solve import SIGNS,TRANS,form_factor
from raw_frames import OUT
TYPES=['C','N','C','C','O','O']
LABELS=['CA','N1','C1','C2','O1','O2']

def molecule(torsion):
 torsion=np.atleast_1d(torsion);b=len(torsion)
 ang=np.deg2rad(109.47);u=np.array([np.sin(ang),0,np.cos(ang)])
 v=np.array([np.cos(ang),0,-np.sin(ang)]);w=np.cross(u,v)
 xyz=np.zeros((b,6,3));xyz[:,1]=[0,0,1.48];xyz[:,2]=1.52*u
 xyz[:,3]=1.53*np.array([np.sin(ang)*np.cos(2*np.pi/3),np.sin(ang)*np.sin(2*np.pi/3),np.cos(ang)])
 ov=np.cos(torsion)[:,None]*v+np.sin(torsion)[:,None]*w
 xyz[:,4]=xyz[:,2]+1.25*(.5*u+np.sqrt(.75)*ov)
 xyz[:,5]=xyz[:,2]+1.25*(.5*u-np.sqrt(.75)*ov)
 return xyz

def coordinates(params,cell):
 pp=np.atleast_2d(params)
 rot=Rotation.from_euler('xyz',pp[:,3:6]).as_matrix()
 xyz=np.einsum('bij,bnj->bni',rot,molecule(pp[:,6]))/cell
 return xyz+pp[:,None,:3]

def main():
 t=time.time();data=np.load(OUT/'merged.npy');cell=np.array(json.loads((OUT/'geometry.json').read_text())['cell'])
 use=(data[:,7]==0)&(data[:,6]>=.85)&(data[:,3]>3*data[:,4])
 d=data[use];h=d[:,:3];fo=np.sqrt(d[:,3]);s2=np.sum((h/cell)**2,axis=1)/4
 ff=np.stack([form_factor(el,s2)*np.exp(-8*np.pi*np.pi*.02*s2) for el in TYPES],axis=1)
 hs=h[:,None,:]*SIGNS[None,:,:];tr=np.exp(2j*np.pi*(h@TRANS.T))
 def calc(coords):
  phase=np.einsum('hsj,baj->bhsa',hs,coords)
  f=(np.exp(2j*np.pi*phase)*tr[None,:,:,None]*ff[None,:,None,:]).sum(axis=(2,3))
  return abs(f)
 def objective(params):
  params=np.asarray(params);scalar=params.ndim==1
  xyz=coordinates(params if scalar else params.T,cell)
  fc=calc(xyz);sc=(fc@fo)/np.sum(fc*fc,axis=1)
  loss=np.sum((sc[:,None]*fc-fo)**2,axis=1)/np.sum(fo*fo)
  return loss[0] if scalar else loss
 bounds=[(0,.5)]*3+[(0,2*np.pi),(-np.pi/2,np.pi/2),(0,2*np.pi),(0,2*np.pi)]
 history=[];best=1e9
 for trial in range(6):
  nit=[0]
  def callback(x,convergence):
   nit[0]+=1
   if nit[0]%100==0:print('search',trial,nit[0],'loss',objective(x),'sec',round(time.time()-t,1),flush=True)
   return objective(x)<.012
  res=differential_evolution(objective,bounds,popsize=24,maxiter=1200,tol=.002,mutation=(.5,1),recombination=.8,rng=np.random.default_rng(220926+trial),vectorized=True,updating='deferred',polish=True,callback=callback)
  loss=objective(res.x);history.append(dict(trial=trial,loss=float(loss),nfev=res.nfev,message=res.message,parameters=res.x.tolist(),coordinates=(coordinates(res.x,cell)[0]%1).tolist()))
  print('search result',trial,loss,flush=True)
  if loss<best:
   best=loss;xyz=coordinates(res.x,cell)[0]
   result=dict(method='direct-space differential evolution of ideal connectivity; no reference coordinates',trial=trial,loss=float(loss),coordinates=(xyz%1).tolist(),types=TYPES,labels=LABELS,cell=cell.tolist(),parameters=res.x.tolist(),n_reflections=len(d),arbitrary_hand=True)
   (OUT/'molecular_solution.json').write_text(json.dumps(result,indent=2))
  if trial>=1 and best<.025:break
 (OUT/'molecular_search_trials.json').write_text(json.dumps(history,indent=2))
if __name__=='__main__':main()
