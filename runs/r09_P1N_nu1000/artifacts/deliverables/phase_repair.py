"""Experimental centric low-resolution sign optimization using pore support.
Ising objective is the band-limited solvent density leaking into the framework region.
Only working reflections enter sign selection; test reflections remain unobserved.
"""
from pore_density import *
from anisotropic import AnisoModel
import argparse


def run(tag='phase_repair',input_tag='framework_h_highangle',cutoff=2.0,starts=20):
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q']
    r=json.loads((OUT/(input_tag+'.json')).read_text())
    model=AnisoModel(r['atoms'],h,True,hydrogens=r.get('hydrogens_included',False));p=model.p0;p[-1]=np.log(r['scale']);p[-3]=r['Zr_fp'];p[-2]=r['Zr_fpp_squared']
    ic,frame,zr,scale,fp,t=model.calc(p,details=True)
    work=np.random.default_rng(712711).random(len(h))>.07
    low=q<=1/cutoff**2;use=low&work;ids=np.where(use)[0];hh=h[ids]
    # Real part of centric structure factor after the modeled imaginary contribution.
    fo=np.sqrt(np.maximum(I/scale**2-t*zr**2,0));ff=fo[ids];fc=frame[ids]
    pore=PoreDensity(r['atoms'],h,q,shape=(192,192,96),cutoff=cutoff)
    ahat=ifftn(1-pore.mask).real
    eqs=[np.unique(x@R,axis=0) for x in hh];mult=np.array([len(x) for x in eqs])
    n=len(hh);K=np.empty((n,n));k0=np.empty(n)
    for j,eq in enumerate(eqs):
        diff=hh[:,None,:]-eq[None,:,:];vals=ahat[tuple((diff%pore.shape).transpose(2,0,1))]
        K[:,j]=mult*np.sum(vals,axis=1)
        k0[j]=ahat[tuple((eq%pore.shape).T)].sum()
    K=(K+K.T)/2;K-=np.outer(k0,k0)/ahat[0,0,0]
    diag=np.diag(K);rng=np.random.default_rng(913419);s0=np.sign(fc);s0[s0==0]=1
    best=None;hist=[]
    for start in range(starts):
        s=s0.copy()
        if start:
            ambiguous=np.where(abs(fc)<1.5*ff)[0];flip=rng.choice(ambiguous,size=max(1,len(ambiguous)//4),replace=False);s[flip]*=-1
        c=s*ff-fc;g=K@c;E=float(c@g)
        dd=-4*s*ff*g+4*ff**2*diag;T0=max(1,float(np.percentile(abs(dd),60))*.2)
        for sweep in range(100):
            T=T0*(1e-5**(sweep/99))
            for j in rng.permutation(n):
                dc=-2*s[j]*ff[j];de=2*dc*g[j]+dc*dc*diag[j]
                if de<0 or rng.random()<np.exp(-min(de/T,700)):
                    c[j]+=dc;s[j]*=-1;g+=dc*K[:,j];E+=de
        for sweep in range(20):
            count=0
            for j in rng.permutation(n):
                dc=-2*s[j]*ff[j];de=2*dc*g[j]+dc*dc*diag[j]
                if de<-1e-7:c[j]+=dc;s[j]*=-1;g+=dc*K[:,j];E+=de;count+=1
            if count==0:break
        DC=float(-k0@c/ahat[0,0,0]);record={'start':start,'energy':float(c@K@c),'phase_flips':int(np.sum(s!=s0)),'DC_estimate':DC};hist.append(record)
        print('SIGNS',record,flush=True)
        if best is None or record['energy']<best['energy']:
            best={**record,'s':s.copy(),'c':c.copy()}
    target=np.zeros(len(h));target[ids]=best['c']
    grid=np.zeros(pore.shape,complex);grid[pore.grid_index]=target[pore.source];grid[0,0,0]=best['DC_estimate']
    pore.rho=(np.maximum(fftn(grid).real/V,0)*pore.mask).astype(np.float32)
    baseline=metrics_local(I,sig,ic,low,work)
    history=[]
    for it in range(150):
        fs=pore.coefficients();delta=np.zeros(len(h));delta[ids]=target[ids]-fs[ids]
        grid=np.zeros(pore.shape,complex);grid[pore.grid_index]=delta[pore.source]
        dr=fftn(grid).real/V
        pore.rho=(np.maximum(pore.rho+.8*dr,0)*pore.mask).astype(np.float32)
        if it%10==0 or it==149:
            fs=pore.coefficients();icorr=scale**2*((frame+fs)**2+t*zr**2)
            record={'iteration':it,**metrics_local(I,sig,icorr,low,work),'electrons':float(pore.rho.mean()*V),'maxrho':float(pore.rho.max())};history.append(record);print('RECONSTRUCT',record,flush=True)
    fs=pore.coefficients();np.savez_compressed(OUT/(tag+'.npz'),h=h,rho=pore.rho,mask=pore.mask,fs=fs,work=work,low=low,signs=best['s'],ids=ids,target=target)
    report={'input':input_tag,'cutoff':cutoff,'baseline':baseline,'sign_trials':hist,'reconstruction':history}
    (OUT/(tag+'.json')).write_text(json.dumps(report,indent=2))

def metrics_local(I,sig,ic,low,work):
    fo=np.sqrt(np.maximum(I,0));fc=np.sqrt(np.maximum(ic,0));out={}
    for lab,mask in [('work',low&work),('test',low&~work)]:out['R_'+lab]=float(np.sum(abs(fo[mask]-fc[mask]))/np.sum(fo[mask]))
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--tag',default='phase_repair');ap.add_argument('--input',default='framework_h_highangle');ap.add_argument('--cutoff',type=float,default=2.);ap.add_argument('--starts',type=int,default=20);a=ap.parse_args()
    run(a.tag,a.input,a.cutoff,a.starts)
