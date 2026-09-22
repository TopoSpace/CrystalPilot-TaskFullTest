"""Symmetry-constrained independent atom least squares and Fourier maps, written in run."""
from crystal import *
from scipy.optimize import least_squares
from scipy.linalg import null_space
import sys
AI=np.linalg.inv(A)


def site_setup(x,tolerance=.12):
    x=np.array(x,float)
    xx=np.einsum('sij,j->si',R,x)
    dif=xx-x; trans=np.rint(dif); dif-=trans
    stabilizer=R[np.linalg.norm(dif@A.T,axis=1)<tolerance]
    # Symmetrize near-special sites using local lattice translations.
    xx=np.einsum('sij,j->si',stabilizer,x)
    xx-=np.rint(xx-x)
    x=xx.mean(axis=0)%1
    constraints=np.vstack([(r-np.eye(3))@AI for r in stabilizer])
    basis=AI@null_space(constraints,rcond=1e-7)
    orbit=np.einsum('sij,j->si',R,x)%1
    _,idx=np.unique(np.round(orbit,7),axis=0,return_index=True)
    return x,basis,R[idx],stabilizer


def atom(label,el,x,U=.05,occ=1,tolerance=.12):
    x,basis,reps,stab=site_setup(x,tolerance)
    return {'label':label,'element':el,'xyz':x.tolist(),'Uiso':U,'occ':occ,
            'basis':basis.tolist(),'reps':reps.tolist(),'stabilizer':stab.tolist()}

class Model:
    def __init__(self,atoms,h):
        self.atoms=atoms;self.h=h;self.q=np.einsum('ni,ij,nj->n',h,GI,h);self.n=len(h)
        self.factors={e:formfactor(e,self.q/4) for e in set(a['element'] for a in atoms)}
        self.spec=[];self.p0=[];self.lo=[];self.hi=[]
        for a in atoms:
            basis=np.array(a['basis']);reps=np.array(a['reps']);x=np.array(a['xyz'])
            hr=np.einsum('ni,sij->nsj',h,reps)
            nd=basis.shape[1];offset=len(self.p0)
            self.spec.append((a,hr,2*np.pi*(hr@basis),2*np.pi*(hr@x),offset,nd))
            self.p0.extend([0]*nd+[a['Uiso']]);self.lo.extend([-.5]*nd+[.005]);self.hi.extend([.5]*nd+[.5])
        self.p0=np.array(self.p0+[-3.5]);self.lo=np.array(self.lo+[-8]);self.hi=np.array(self.hi+[1])

    def calc(self,p,jac=False):
        f=np.zeros(self.n);J=np.zeros((self.n,len(p))) if jac else None
        for a,hr,hb,phase,off,nd in self.spec:
            ph=phase+hb@p[off:off+nd]
            co=np.cos(ph).sum(axis=1)
            dw=np.exp(-2*np.pi**2*self.q*p[off+nd])
            fac=self.factors[a['element']]*dw*a['occ']
            term=fac*co;f+=term
            if jac:
                if nd:J[:,off:off+nd]=fac[:,None]*np.einsum('ns,nsj->nj',-np.sin(ph),hb)
                J[:,off+nd]=term*(-2*np.pi**2*self.q)
        scale=np.exp(p[-1]);f*=scale
        if jac:
            J*=scale;J[:,-1]=f
            return f,J
        return f

    def update(self,p):
        out=[]
        for a,hr,hb,phase,off,nd in self.spec:
            a=dict(a);a['xyz']=((np.array(a['xyz'])+np.array(a['basis'])@p[off:off+nd])%1).tolist();a['Uiso']=float(p[off+nd]);out.append(a)
        return out

def metrics(I,sig,fc,mask=None):
    if mask is None:mask=np.ones(len(I),bool)
    I=I[mask];sig=sig[mask];fc=fc[mask];fo=np.sqrt(np.maximum(I,0));obs=I>2*sig
    w=1/(sig**2+(.05*np.maximum(I,0))**2+.01)
    return {'n':len(I),'n_obs':int(obs.sum()),'R1_all':float(np.sum(abs(fo-abs(fc)))/np.sum(fo)),
            'R1_obs':float(np.sum(abs(fo[obs]-abs(fc[obs])))/np.sum(fo[obs])),
            'wR2':float(np.sqrt(np.sum(w*(I-fc**2)**2)/np.sum(w*I**2)))}

def refine(atoms,scale=.025,dmin=1.,dmax=8.,tag='refined',nfev=100,aweight=.08):
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q']
    mask=(q<=1/dmin**2)&(q>=1/dmax**2)
    hm=h[mask];im=I[mask];sm=sig[mask]
    model=Model(atoms,hm);p=model.p0.copy();p[-1]=np.log(scale)
    weights=1/np.sqrt(sm**2+(aweight*np.maximum(im,0))**2+.01)
    def fun(p):
        f=model.calc(p)
        return (f*f-im)*weights
    def jac(p):
        f,j=model.calc(p,True)
        return 2*f[:,None]*j*weights[:,None]
    print('Refining',tag,'atoms',len(atoms),'parameters',len(p),'reflections',len(hm),flush=True)
    opt=least_squares(fun,p,jac=jac,bounds=(model.lo,model.hi),x_scale='jac',max_nfev=nfev,ftol=1e-8)
    new=model.update(opt.x);scale=float(np.exp(opt.x[-1]));mall=Model(new,h);pa=mall.p0;pa[-1]=np.log(scale);fc=mall.calc(pa)
    report={'atoms':new,'scale':scale,'n_params':len(p),'dmin':dmin,'dmax':dmax,'metrics_fit':metrics(I,sig,fc,mask),'metrics_all':metrics(I,sig,fc),'cost':float(opt.cost),'success':bool(opt.success),'message':opt.message,'nfev':opt.nfev}
    (OUT/(tag+'.json')).write_text(json.dumps(report,indent=2));np.savez_compressed(OUT/(tag+'_fc.npz'),h=h,fc=fc)
    print(report['metrics_fit'],'scale',scale,'eval',opt.nfev,flush=True)
    for a in new:print(a['label'],a['element'],np.round(a['xyz'],6),'U',round(a['Uiso'],4),'mult',len(a['reps']),flush=True)
    return report,fc

def maps(report,fc,tag,dmin=1.,npeak=160):
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q'];m=q<1/dmin**2
    fo=np.sqrt(np.maximum(I,0));phase=np.sign(fc)
    for typ,coef in [('fo',fo*phase),('diff',fo*phase-fc),('2fo',(2*fo*phase-fc))]:
        rho=map_fft(h[m],coef[m]/report['scale'],shape=(240,240,120))
        save_peaks(peaks(rho,npeak,.7),tag+'_'+typ+'_peaks.csv')
        np.save(OUT/(tag+'_'+typ+'_map.npy'),rho.astype(np.float32))
    return rho

def representatives_near(x,target=[.5,.12,.24]):
    eq=expand_xyz(x)
    dd=eq-target;dd-=np.round(dd)
    j=np.argmin(np.linalg.norm(dd@A.T,axis=1))
    return np.array(target)+dd[j]

def inspect_peaks(file,atoms,n=50):
    pk=np.loadtxt(OUT/file,delimiter=',',skiprows=1)
    allxyz=[];labels=[]
    for a in atoms:
        xx=expand_xyz(a['xyz']);allxyz.extend(xx);labels.extend([a['label']]*len(xx))
    allxyz=np.array(allxyz)
    for i,row in enumerate(pk[:n]):
        x=representatives_near(row[:3]);ds=mindist(x,allxyz);order=np.argsort(ds)[:3]
        print(i+1,np.round(x,6),'height',round(row[3],2),'near',[(labels[j],round(ds[j],3)) for j in order],flush=True)

if __name__=='__main__':
    z=json.loads((OUT/'zr_solution.json').read_text())
    atoms=[atom('Zr'+str(i+1),'Zr',x,z['B']/(8*np.pi**2),tolerance=.02) for i,x in enumerate(z['sites'])]
    report,fc=refine(atoms,z['amplitude_scale'],tag='zr_refined',dmin=1.,dmax=7.,nfev=100)
    maps(report,fc,'zr_refined')
    inspect_peaks('zr_refined_diff_peaks.csv',report['atoms'],60)
