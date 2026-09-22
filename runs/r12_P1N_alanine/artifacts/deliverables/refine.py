"""Independent full-matrix least-squares refinement on F^2 with riding H.
Neutral-atom four-Gaussian form factors are explicit in solve.py.
No externally supplied atom coordinates or refinement engine is used.
"""
import json,itertools,time
import numpy as np
from scipy.optimize import least_squares
from solve import SIGNS,TRANS,form_factor
from raw_frames import OUT
TYPES=['C','N','C','C','O','O'];LABELS=['CA','N1','C1','C2','O1','O2']
H_LABELS=['HA','HN1','HN2','HN3','HM1','HM2','HM3']

def unit(x):return x/np.linalg.norm(x)

def make_hydrogens(x,cell,torsions):
 # Molecule represented with whole-molecule coordinates, not wrapped atomwise.
 c=x*cell;ca,n,co,cm,o1,o2=c
 ha=ca+.98*unit(-sum(unit(v-ca) for v in [n,co,cm]))
 hyd=[ha]
 for center,ang,length in [(n,torsions[0],.91),(cm,torsions[1],.98)]:
  axis=unit(ca-center);ref=co-ca;v=unit(ref-np.dot(ref,axis)*axis);w=np.cross(axis,v)
  for k in range(3):
   t=ang+2*np.pi*k/3
   direction=-axis/3+np.sqrt(8/9)*(v*np.cos(t)+w*np.sin(t))
   hyd.append(center+length*direction)
 return np.array(hyd)/cell

def unpack(p,mode,cell,hydrogen=True):
 x=p[:18].reshape(6,3)
 if mode=='iso':
  u=np.array([np.eye(3)*v for v in p[18:24]]);logscale=p[24];torsions=p[25:27] if hydrogen else [0,0]
 else:
  v=p[18:54].reshape(6,6);l=np.zeros((6,3,3))
  l[:,0,0]=v[:,0];l[:,1,1]=v[:,1];l[:,2,2]=v[:,2];l[:,1,0]=v[:,3];l[:,2,0]=v[:,4];l[:,2,1]=v[:,5]
  u=l@l.transpose(0,2,1);logscale=p[54];torsions=p[55:57] if hydrogen else [0,0]
 if hydrogen:
  hx=make_hydrogens(x,cell,torsions);ueq=np.trace(u,axis1=1,axis2=2)/3
  hu=np.array([ueq[0]*1.2]+[ueq[1]*1.5]*3+[ueq[3]*1.5]*3)
  x=np.r_[x,hx];u=np.r_[u,np.array([np.eye(3)*v for v in hu])]
 return x,u,float(np.exp(logscale)),torsions

class StructureFactors:
 def __init__(self,h,cell,hydrogen=True):
  self.h=np.asarray(h);self.cell=np.asarray(cell);self.hs=self.h[:,None,:]*SIGNS[None,:,:]
  self.q=self.hs/self.cell;self.trans=np.exp(2j*np.pi*(self.h@TRANS.T));s2=np.sum((self.h/self.cell)**2,axis=1)/4
  self.ff=np.stack([form_factor(el,s2) for el in TYPES+(['H']*7 if hydrogen else [])],axis=1)
 def calculate(self,x,u):
  phase=np.einsum('hsj,aj->hsa',self.hs,x)
  dw=np.exp(-2*np.pi**2*np.einsum('hsi,aij,hsj->hsa',self.q,u,self.q))
  return np.sum(np.exp(2j*np.pi*phase)*dw*self.trans[:,:,None]*self.ff[:,None,:],axis=(1,2))

def rvalues(I,sig,fc,k,npar,weights=None):
 calc=k*abs(fc)**2;obsamp=np.sqrt(np.maximum(I,0));calcamp=np.sqrt(calc)
 sel=I>2*sig
 if weights is None:weights=1/(sig**2+(.04*np.maximum(I,0))**2)
 return dict(R1_gt2sigma=float(np.sum(abs(obsamp[sel]-calcamp[sel]))/np.sum(obsamp[sel])),R1_all=float(np.sum(abs(obsamp-calcamp))/np.sum(obsamp)),wR2=float(np.sqrt(np.sum(weights*(I-calc)**2)/np.sum(weights*I*I))),goodness_of_fit=float(np.sqrt(np.sum(weights*(I-calc)**2)/(len(I)-npar))),n_reflections=len(I),n_gt2sigma=int(sel.sum()),n_parameters=npar,scale=k)

def fit_stage(p,mode,h,I,sig,cell,hydrogen=True,max_nfev=180):
 sf=StructureFactors(h,cell,hydrogen)
 weights=1/(sig**2+(.03*np.maximum(I,0))**2)
 def fun(p):
  x,u,k,_=unpack(p,mode,cell,hydrogen)
  return (k*abs(sf.calculate(x,u))**2-I)*np.sqrt(weights)
 lb=np.full(len(p),-np.inf);ub=np.full(len(p),np.inf)
 if mode=='iso':lb[18:24]=.001;ub[18:24]=.15;scale=np.r_[[.01]*18,[.015]*6,1,[1]*(len(p)-25)]
 else:
  for a in range(6):lb[18+a*6:21+a*6]=.015;ub[18+a*6:21+a*6]=.5
  scale=np.r_[[.01]*18,[.05]*36,1,[1]*(len(p)-55)]
 res=least_squares(fun,p,bounds=(lb,ub),x_scale=scale,loss='linear',max_nfev=max_nfev,ftol=1e-10,xtol=1e-10,gtol=1e-8)
 x,u,k,tors=unpack(res.x,mode,cell,hydrogen);fc=sf.calculate(x,u);stats=rvalues(I,sig,fc,k,len(res.x),weights)
 stats['optimizer_message']=res.message;stats['optimality']=res.optimality
 print(mode,'H',hydrogen,json.dumps(stats),flush=True)
 return res,stats,fc,weights

def main():
 cell=np.array(json.loads((OUT/'geometry.json').read_text())['cell']);sol=json.loads((OUT/'molecular_solution.json').read_text())
 data=np.load(OUT/'merged.npy');use=(data[:,7]==0)&np.isfinite(data[:,3])&np.isfinite(data[:,4])&(data[:,4]>0)
 data=data[use];h=data[:,:3];I=data[:,3];sig=data[:,4]
 x=np.array(sol['coordinates']);x=x[0]+(x-x[0]+.5)%1-.5
 uiso=np.full(6,.02);sf=StructureFactors(h,cell,False);trials=[];best=np.inf
 # C and N differ by only one electron. Explicitly compare N/methyl-carbon
 # assignments; the molecular hand is not known from these Friedel-merged data.
 for permutation in [[0,1,2,3,4,5],[0,3,2,1,4,5]]:
  xt=x[permutation];fc=sf.calculate(xt,np.array([np.eye(3)*v for v in uiso]));k=np.sum(np.maximum(I,0))/np.sum(abs(fc)**2)
  p=np.r_[xt.ravel(),uiso,np.log(k)]
  rr,st,ff,ww=fit_stage(p,'iso',h,I,sig,cell,False)
  chi=float(rr.fun@rr.fun);trials.append(dict(input_atom_permutation=permutation,chi_squared=chi,statistics=st))
  if chi<best:best=chi;res0,stats0,fc0,w=rr,st,ff,ww;chosen_permutation=permutation
 (OUT/'assignment_trials.json').write_text(json.dumps(dict(chosen_input_atom_permutation=chosen_permutation,trials=trials),indent=2))
 p=np.r_[res0.x,[0,0]];best=1e99
 sfh=StructureFactors(h,cell,True)
 for t1,t2 in itertools.product(np.arange(12)*2*np.pi/36,repeat=2):
  pp=p.copy();pp[-2:]=[t1,t2];xx,uu,kk,_=unpack(pp,'iso',cell,True)
  chi=np.sum(w*(I-kk*abs(sfh.calculate(xx,uu))**2)**2)
  if chi<best:best=chi;pbest=pp
 res1,stats1,fc1,w=fit_stage(pbest,'iso',h,I,sig,cell,True)
 xx,uu,kk,tors=unpack(res1.x,'iso',cell,True)
 chol=np.zeros((6,6));chol[:,:3]=np.sqrt(np.diagonal(uu[:6],axis1=1,axis2=2))
 p=np.r_[res1.x[:18],chol.ravel(),np.log(kk),tors]
 res,stats,fc,weights=fit_stage(p,'aniso',h,I,sig,cell,True)
 xx,uu,kk,tors=unpack(res.x,'aniso',cell,True)
 dof=len(I)-len(res.x);chi2=float(np.sum(res.fun**2));cov=np.linalg.pinv(res.jac.T@res.jac,rcond=1e-12)*(chi2/dof)
 esd=np.sqrt(np.diag(cov)[:18]).reshape(6,3)
 model=dict(input_atom_permutation=chosen_permutation,cell=cell.tolist(),space_group='P 21 21 21',hall='P 2ac 2ab',Z=4,formula='C3 H7 N O2',labels=LABELS+H_LABELS,types=TYPES+['H']*7,coordinates=xx.tolist(),U_cartesian=uu.tolist(),coordinate_su_nonH=esd.tolist(),parameters=res.x.tolist(),scale=kk,torsions=np.asarray(tors).tolist(),statistics=stats,stages=[stats0,stats1,stats],weight_model='1/(sigma(I)^2 + (0.03*max(I,0))^2)',absolute_configuration='undetermined; arbitrary template hand; Friedel equivalents merged; no anomalous refinement',hydrogen_model='riding ideal geometry, two refined torsions; CA-H=0.98, N-H=0.91, methyl-H=0.98 A',constraints='H riding and Uiso multipliers only; all non-H xyz and positive-definite anisotropic U free',temperature_K=150.0)
 (OUT/'refined_model.json').write_text(json.dumps(model,indent=2))
 np.save(OUT/'refinement_covariance.npy',cov)
 np.savetxt(OUT/'refinement_reflections.csv',np.c_[h,I,sig,kk*abs(fc)**2,fc.real,fc.imag,weights,res.fun],delimiter=',',header='h,k,l,Fo2,sigmaFo2,kFc2,Fc_real,Fc_imag,weight,weighted_residual',comments='',fmt='%.10g')
 (OUT/'refinement_stages.json').write_text(json.dumps([stats0,stats1,stats],indent=2))
 # Bond distances and standard uncertainties from the full coordinate covariance.
 bonds=[]
 for i,j in [(0,1),(0,2),(0,3),(2,4),(2,5)]:
  dv=(xx[j]-xx[i])*cell;distance=np.linalg.norm(dv);grad=np.zeros(len(res.x));grad[j*3:j*3+3]=dv/distance*cell;grad[i*3:i*3+3]=-dv/distance*cell
  su=float(np.sqrt(grad@cov@grad));bonds.append(dict(atom1=LABELS[i],atom2=LABELS[j],distance=float(distance),su=su))
 (OUT/'bond_geometry.json').write_text(json.dumps(bonds,indent=2));print('bonds',bonds,flush=True)
 # Difference Fourier, with unmeasured coefficients zero, no invented F000.
 shape=np.array([48,48,96]);grid=np.zeros(shape,complex)
 allh=np.unique(np.concatenate([h*s for s in itertools.product([-1,1],repeat=3)]),axis=0).astype(int)
 vals={tuple(map(int,h0)):i for i,h0 in enumerate(h)}
 calc=StructureFactors(allh,cell,True).calculate(xx,uu)
 obs=np.array([np.sqrt(max(I[vals[tuple(abs(v))]],0)/kk) for v in allh])
 coeff=(obs-abs(calc))*calc/np.maximum(abs(calc),1e-10)
 grid[tuple((allh%shape).T)]=coeff
 rho=np.fft.fftn(grid).real/np.prod(cell)
 np.save(OUT/'difference_density.npy',rho)
 imin=np.unravel_index(np.argmin(rho),rho.shape);imax=np.unravel_index(np.argmax(rho),rho.shape)
 diff=dict(min_e_A3=float(rho.min()),max_e_A3=float(rho.max()),rms_e_A3=float(rho.std()),min_fractional=(np.array(imin)/shape).tolist(),max_fractional=(np.array(imax)/shape).tolist(),grid=shape.tolist(),missing_coefficients='zero',F000='zero in difference map')
 (OUT/'difference_statistics.json').write_text(json.dumps(diff,indent=2));print('difference',diff,flush=True)
if __name__=='__main__':main()
