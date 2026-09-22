"""Low-dimensional NON-ATOMIC pore density: symmetry-related Gaussian components.
These are continuous electron-density terms, not assigned solvent/guest atoms.
All component centers are selected away from the solved framework. Withheld
reflections are never used during component or framework optimization.
"""
from restraints import *
from pore_density import PoreDensity,ifftn
import argparse

class GaussianDensity:
    def __init__(self,components,h,mask_ft):
        self.components=components;self.h=h;self.mask_ft=mask_ft;self.q=np.einsum('ni,ij,nj->n',h,GI,h)
        self.spec=[];p=[];lo=[];hi=[]
        for a in components:
            basis=np.array(a['basis']);hr=np.einsum('ni,sij->nsj',h,a['reps']);nd=basis.shape[1];off=len(p)
            self.spec.append((a,2*np.pi*(hr@basis),2*np.pi*(hr@np.array(a['xyz'])),off,nd))
            p.extend([0]*nd+[a.get('variance',.45),np.log(a.get('electrons',15.))]);lo.extend([-.65]*nd+[.08,np.log(.05)]);hi.extend([.65]*nd+[4.,np.log(250.)])
        p.extend([.15,4.]);lo.extend([0,0]);hi.extend([.8,60])
        self.p0=np.array(p);self.lo=np.array(lo);self.hi=np.array(hi)
    def calc(self,p,jac=False):
        f=np.zeros(len(self.h));J=np.zeros((len(f),len(p))) if jac else None
        for a,hb,ph0,off,nd in self.spec:
            ph=ph0+hb@p[off:off+nd];co=np.cos(ph).sum(axis=1)
            fac=np.exp(p[off+nd+1]-2*np.pi**2*self.q*p[off+nd]);term=fac*co;f+=term
            if jac:
                J[:,off:off+nd]=fac[:,None]*np.einsum('ns,nsj->nj',-np.sin(ph),hb)
                J[:,off+nd]=-2*np.pi**2*self.q*term;J[:,off+nd+1]=term
        term=self.mask_ft*np.exp(-p[-1]*self.q/4);f+=p[-2]*term
        if jac:J[:,-2]=term;J[:,-1]=-self.q/4*p[-2]*term;return f,J
        return f
    def update(self,p):
        cs=[]
        for a,hb,ph0,off,nd in self.spec:
            a=dict(a);a['xyz']=((np.array(a['xyz'])+np.array(a['basis'])@p[off:off+nd])%1).tolist();a['variance']=float(p[off+nd]);a['electrons']=float(np.exp(p[off+nd+1]));cs.append(a)
        return cs

def candidates(inp,filename,nmax=14):
    if nmax==0:return []
    pk=np.loadtxt(OUT/filename,delimiter=',',skiprows=1);framework=np.vstack([expand_xyz(a['xyz']) for a in inp['atoms']]);cs=[];expanded=[]
    for row in pk:
        x=representatives_near(row[:3]);d=float(min(mindist(x,framework)))
        if d<3.1:continue
        if expanded and min(mindist(x,np.array(expanded)))<1.5:continue
        a=atom('Q'+str(len(cs)+1),'C',x,tolerance=.15);a['element']='density';a['source_peak_height']=float(row[3]);a['framework_distance_A']=d;a['variance']=.45;a['electrons']=10.
        cs.append(a);expanded.extend(expand_xyz(a['xyz']))
        if len(cs)==nmax:break
    return cs

def run(input_tag='framework_h_highangle',tag='framework_density',ncomp=14,epochs=5,radius_scale=1.):
    inp=json.loads((OUT/(input_tag+'.json')).read_text());seed=json.loads((OUT/'framework_iso1.json').read_text())['atoms']
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q'];work=np.random.default_rng(712711).random(len(h))>.07
    fitrange=q<=1.;fit=fitrange&work
    cs=candidates(inp,input_tag+'_diff_peaks.csv',ncomp)
    print('DENSITY CANDIDATES',[(a['label'],np.round(a['xyz'],5).tolist(),len(a['reps']),a['source_peak_height']) for a in cs],flush=True)
    pore=PoreDensity(inp['atoms'],h,q,radius_scale=radius_scale)
    maskft=(ifftn(pore.mask)*V)[pore.h_index].real
    bg=GaussianDensity(cs,h[fit],maskft[fit]);bgall=GaussianDensity(cs,h,maskft)
    m=AnisoModel(inp['atoms'],h[fit],True,hydrogens=True);mall=AnisoModel(inp['atoms'],h,True,hydrogens=True);nf=len(m.p0)
    p=np.r_[m.p0,bg.p0];p[nf-1]=np.log(inp['scale']);p[nf-3]=inp['Zr_fp'];p[nf-2]=inp['Zr_fpp_squared']
    lo=np.r_[m.lo,bg.lo];hi=np.r_[m.hi,bg.hi]
    geo=Geometry(m,seed,strength=2.);simu=geo.adp_matrix();im=I[fit];sm=sig[fit]
    def calc(p,jac=False,full=False):
        mm=mall if full else m;bb=bgall if full else bg
        if jac:sol,js=bb.calc(p[nf:],True)
        else:sol=bb.calc(p[nf:])
        mm.solvent=sol
        if not jac:return mm.calc(p[:nf])
        ic,J=mm.calc(p[:nf],True);_,real,zr,scale,fp,t=mm.calc(p[:nf],details=True)
        return ic,np.column_stack([J,2*scale**2*real[:,None]*js])
    # Initially keep framework fixed while continuous density obtains its scale.
    for stage in range(epochs+1):
        current=calc(p);P=(np.maximum(im,0)+2*current)/3;w=1/np.sqrt(sm**2+(.08*P)**2+.01)
        def fun(p):return np.r_[(calc(p)-im)*w,geo.fun(p[:nf]),simu@p[:nf],10*m.adp_restraints(p[:nf])]
        def jac(p):
            ic,J=calc(p,True);rest=np.vstack([geo.jac(p[:nf]),simu,10*m.adp_restraints(p[:nf],True)])
            return np.vstack([J*w[:,None],np.pad(rest,((0,0),(0,len(p)-nf)))])
        if stage==0:
            def f0(pb):return (calc(np.r_[p[:nf],pb])-im)*w
            def j0(pb):return calc(np.r_[p[:nf],pb],True)[1][:,nf:]*w[:,None]
            opt=least_squares(f0,p[nf:],jac=j0,bounds=(lo[nf:],hi[nf:]),max_nfev=150,x_scale='jac',ftol=1e-7);p[nf:]=opt.x
        else:
            opt=least_squares(fun,p,jac=jac,bounds=(lo,hi),max_nfev=100,x_scale='jac',ftol=1e-7);p=opt.x
        ic=calc(p,full=True);rec={'stage':stage,'work':metrics(I,sig,np.sqrt(ic),fitrange&work),'test':metrics(I,sig,np.sqrt(ic),fitrange&~work),'cost':float(opt.cost),'nfev':opt.nfev}
        print('GAUSSIAN',rec,flush=True)
        if stage==0:history=[rec]
        else:history.append(rec)
        atoms=m.update(p[:nf]);components=bg.update(p[nf:]);mall.solvent=bgall.calc(p[nf:]);ic,real,zr,scale,fp,t=mall.calc(p[:nf],details=True)
        report={'atoms':atoms,'riding_hydrogens':m.update_hydrogens(p[:nf]),'scale':float(scale),'Zr_fp':float(fp),'Zr_fpp_squared':float(t),'dispersion_refined':True,'hydrogens_included':True,'geometry_strength':2.,'n_params':len(p),'dmin':1.,'dmax':50.,
            'metrics_fit':metrics(I,sig,np.sqrt(ic),fit),'metrics_work':rec['work'],'metrics_test':rec['test'],'metrics_all_range':metrics(I,sig,np.sqrt(ic),fitrange),
            'mask_radius_scale':radius_scale,'components':components,'uniform_density_e_A3':float(p[-2]),'uniform_B':float(p[-1]),'pore_volume_fraction':pore.volfrac,'history':history,'input_model':input_tag,
            'estimated_pore_electrons':float(sum(a['electrons']*len(a['reps']) for a in components)+p[-2]*pore.mask.mean()*V),'component_description':'Unassigned continuous Gaussian electron densities; NOT solvent atomic species.'}
        (OUT/(tag+'.json')).write_text(json.dumps(report,indent=2));np.savez_compressed(OUT/(tag+'_fc.npz'),h=h,ic=ic,real=real,zr=zr,scale=scale,fp=fp,t=t,solvent=mall.solvent,mask_ft=maskft,work=work,fit=fit)
    np.savez_compressed(OUT/(tag+'_covariance.npz'),cov=np.linalg.pinv(opt.jac.T@opt.jac,rcond=1e-10),parameters=p,jac=opt.jac,residual=opt.fun)
    fo=np.sqrt(np.maximum(I,0))/scale;amp=np.sqrt(real**2+t*zr**2);delta=(fo/np.maximum(amp,1e-8)-1)*real
    rho=map_fft(h[fitrange],delta[fitrange],shape=(240,240,120));np.save(OUT/(tag+'_diff_map.npy'),rho.astype(np.float32));save_peaks(peaks(rho,120,.7),tag+'_diff_peaks.csv')
    print('FINAL',report['metrics_all_range'],'pore_e',report['estimated_pore_electrons'],flush=True)
    for a in atoms:print(a['label'],np.round(a['xyz'],6),'Ueig',np.round(a['U_eigenvalues'],4),flush=True)
    for a in components:print(a['label'],np.round(a['xyz'],6),'electrons',a['electrons'],'variance',a['variance'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='framework_h_highangle');ap.add_argument('--tag',default='framework_density');ap.add_argument('--components',type=int,default=14);ap.add_argument('--epochs',type=int,default=5);ap.add_argument('--radius-scale',type=float,default=1.);a=ap.parse_args();run(a.input,a.tag,a.components,a.epochs,a.radius_scale)
