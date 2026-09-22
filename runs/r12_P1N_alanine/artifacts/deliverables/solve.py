"""Ab initio charge flipping in P212121, using only measured intensities.
No molecular coordinates or reference structure are supplied to the solver.
"""
import json,time
import numpy as np
from scipy import ndimage
from raw_frames import OUT
SIGNS=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])
TRANS=np.array([[0,0,0],[.5,.5,0],[0,.5,.5],[.5,0,.5]])
COEFF={
 'C':([2.31,1.02,1.5886,.865],[20.8439,10.2075,.5687,51.6512],.2156),
 'N':([12.2126,3.1322,2.0125,1.1663],[.0057,9.8933,28.9975,.5826],-11.529),
 'O':([3.0485,2.2868,1.5463,.867],[13.2771,5.7011,.3239,32.9089],.2508),
 'H':([.493002,.322912,.140191,.04081],[10.5109,26.1257,3.14236,57.7997],.003038)
}
def form_factor(el,s2):
 a,b,c=COEFF[el];return np.sum(np.array(a)[:,None]*np.exp(-np.array(b)[:,None]*s2[None,:]),axis=0)+c

def symmetry_positions(x):return (x[:,None,:]*SIGNS[None,:,:]+TRANS[None,:,:])%1

def peak_atoms(rho,cell,n=6):
 shape=np.array(rho.shape);mx=ndimage.maximum_filter(rho,size=3,mode='wrap');idx=np.argwhere((rho==mx)&(rho>0));heights=rho[tuple(idx.T)];order=np.argsort(heights)[::-1]
 atoms=[];hs=[]
 for k in order:
  pos=idx[k]/shape
  if atoms:
   sy=symmetry_positions(np.array(atoms)).reshape(-1,3);delta=sy-pos;delta-=np.rint(delta)
   if np.min(np.linalg.norm(delta*cell,axis=1))<1.05:continue
  # Parabolic sub-grid interpolation about the local maximum.
  dd=[]
  for ax in range(3):
   pp=idx[k].copy();mm=idx[k].copy();pp[ax]=(pp[ax]+1)%shape[ax];mm[ax]=(mm[ax]-1)%shape[ax]
   v0=rho[tuple(idx[k])];vp=rho[tuple(pp)];vm=rho[tuple(mm)];den=vm-2*v0+vp
   dd.append(np.clip(.5*(vm-vp)/den if den<0 else 0,-.5,.5))
  pos=(idx[k]+dd)/shape;atoms.append(pos);hs.append(heights[k])
  if len(atoms)==n:break
 return np.array(atoms),np.array(hs)

def main():
 data=np.load(OUT/'merged.npy');cell=np.array(json.loads((OUT/'geometry.json').read_text())['cell'])
 shape=np.array([32,32,64]);nvox=np.prod(shape)
 use=(data[:,7]==0)&(data[:,6]>=.75)&(data[:,3]>0)&(data[:,3]>1.5*data[:,4])
 h=data[use,:3].astype(int);I=data[use,3]
 s2=np.sum((h/cell)**2,axis=1)/4
 ff2=4*(3*form_factor('C',s2)**2+form_factor('N',s2)**2+2*form_factor('O',s2)**2)
 E=np.sqrt(I/ff2)*np.exp(1.0*s2)
 E/=np.sqrt(np.mean(E**2));ampvals=E*np.sqrt(24)
 amp=np.zeros(shape);known=np.zeros(shape,bool)
 for sign in np.array(list(__import__('itertools').product([-1,1],repeat=3))):
  hh=h*sign;idx=tuple((hh%shape).T);amp[idx]=ampvals;known[idx]=True
 grids=np.indices(shape).reshape(3,-1).T
 symidx=[]
 for s,t in zip(SIGNS,TRANS):symidx.append(np.ravel_multi_index(tuple(((grids*s+t*shape).astype(int)%shape).T),shape))
 symidx=np.array(symidx)
 freq=np.array(np.meshgrid(*[np.fft.fftfreq(n)*n for n in shape],indexing='ij'))
 sphere=np.sqrt(np.sum((freq/cell[:,None,None,None])**2,axis=0))<=1/.75
 allowed=sphere.copy();allowed[0,0,0]=False
 rng=np.random.default_rng(20260922);best=1e9;start=time.time();log=[]
 for trial in range(25):
  rho=rng.normal(size=shape);rho=rho.ravel()[symidx].mean(axis=0).reshape(shape)
  F=np.fft.fftn(rho);F=amp*np.exp(1j*np.angle(F));rho=np.fft.ifftn(F).real
  delta=rho.std()*[1.0,1.1,1.2,.9,1.3][trial%5]
  trialbest=1e9
  for it in range(800):
   modified=np.where(rho<delta,-rho,rho)
   G=np.fft.fftn(modified)
   absG=np.abs(G);r=np.sum(np.abs(absG[known]-amp[known]))/np.sum(amp[known])
   # Projection on measured moduli; unmeasured interior coefficients float.
   F=G*allowed;F[known]=amp[known]*np.exp(1j*np.angle(G[known]));F[0,0,0]=0
   rho=np.fft.ifftn(F).real
   if r<trialbest:trialbest=r
   if r<best:
    best=float(r);atoms,heights=peak_atoms(rho,cell)
    np.save(OUT/'solution_density.npy',rho)
    result=dict(trial=trial,iteration=it,charge_flip_R=best,coordinates=atoms.tolist(),peak_heights=heights.tolist(),cell=cell.tolist(),space_group='P 21 21 21',n_observed_used=len(h),seed=20260922)
    (OUT/'solution.json').write_text(json.dumps(result,indent=2))
   if it%200==0:print('CF',trial,it,'R',round(r,4),'best',round(best,4),'sec',round(time.time()-start,1),flush=True)
   if best<.22 and it>300:break
  log.append(dict(trial=trial,best=float(trialbest)))
  if best<.22:break
 (OUT/'solution_trials.json').write_text(json.dumps(log,indent=2))
 print('BEST',best,flush=True)
if __name__=='__main__':main()
