"""Independent implementation of site-symmetry constrained Cartesian ADPs and F^2 LS.
Optional effective Zr dispersion terms are empirical, NOT tabulated edge values.
"""
from refine import *

E6=np.array([[[1,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,1,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,1]],
             [[0,1/np.sqrt(2),0],[1/np.sqrt(2),0,0],[0,0,0]],[[0,0,1/np.sqrt(2)],[0,0,0],[1/np.sqrt(2),0,0]],[[0,0,0],[0,0,1/np.sqrt(2)],[0,1/np.sqrt(2),0]]])

def tensor_basis(stabilizer):
    qs=np.array([A@r@AI for r in stabilizer])
    projected=np.array([np.mean([q@e@q.T for q in qs],axis=0) for e in E6])
    _,s,vh=np.linalg.svd(projected.reshape(6,9),full_matrices=False)
    return vh[s>1e-7].reshape(-1,3,3)

def aromatic_hydrogens(atoms):
    hydrogens=[]
    for i,a in enumerate(atoms):
        if a['element']!='C':continue
        dirs=[]
        for b in atoms:
            if b['element']!='C':continue
            ex=expand_xyz(b['xyz']);d=ex-np.array(a['xyz']);d-=np.round(d);dc=d@A.T;ds=np.linalg.norm(dc,axis=1)
            for j in np.where((ds>.3)&(ds<1.85))[0]:dirs.append(dc[j]/ds[j])
        if len(dirs)!=2:continue
        vec=-np.sum(dirs,axis=0);vec*=.95/np.linalg.norm(vec)
        x=(np.array(a['xyz'])+AI@vec)%1
        h=atom('H'+a['label'][1:],'H',x,1.2*a['Uiso'],tolerance=.015)
        h['parent']=a['label'];h['parent_index']=i;hydrogens.append(h)
    return hydrogens

class AnisoModel:
    def __init__(self,atoms,h,dispersion=True,solvent=None,hydrogens=False):
        self.atoms=atoms;self.h=h;self.q=np.einsum('ni,ij,nj->n',h,GI,h);self.n=len(h);self.dispersion=dispersion
        self.solvent=np.zeros(len(h)) if solvent is None else np.array(solvent)
        self.hydrogens=hydrogens if isinstance(hydrogens,list) else (aromatic_hydrogens(atoms) if hydrogens else [])
        self.factors={e:formfactor(e,self.q/4) for e in set(a['element'] for a in atoms)}
        self.spec=[];p=[];lo=[];hi=[]
        for a in atoms:
            basis=np.array(a['basis']);reps=np.array(a['reps']);x=np.array(a['xyz'])
            hr=np.einsum('ni,sij->nsj',h,reps);kc=hr@AI
            ub=tensor_basis(a['stabilizer']);U=np.array(a.get('Ucart',np.eye(3)*a['Uiso']))
            nd=basis.shape[1];nu=len(ub);off=len(p)
            hU=np.einsum('nsi,tij,nsj->nst',kc,ub,kc)*(-2*np.pi**2)
            self.spec.append({'a':a,'hb':2*np.pi*(hr@basis),'phase':2*np.pi*(hr@x),'hU':hU,'ub':ub,'off':off,'nd':nd,'nu':nu})
            p.extend([0]*nd+np.einsum('tij,ij->t',ub,U).tolist());lo.extend([-.4]*nd+[-.8]*nu);hi.extend([.4]*nd+[.8]*nu)
        self.hspec=[]
        for ah in self.hydrogens:
            parent=self.spec[ah['parent_index']];basis=np.array(parent['a']['basis'])
            hr=np.einsum('ni,sij->nsj',h,ah['reps'])
            self.hspec.append({'parent':parent,'hb':2*np.pi*(hr@basis),'phase':2*np.pi*(hr@np.array(ah['xyz'])),'ff':formfactor('H',self.q/4),'utrace':np.trace(parent['ub'],axis1=1,axis2=2)*.4})
        self.natomic=len(p)
        if dispersion:p.extend([-4,9]);lo.extend([-15,0]);hi.extend([3,150])
        p.append(-3.5);lo.append(-8);hi.append(1)
        self.p0=np.array(p);self.lo=np.array(lo);self.hi=np.array(hi)

    def calc(self,p,jac=False,details=False):
        fp=p[-3] if self.dispersion else 0.;t=p[-2] if self.dispersion else 0.
        real=self.solvent.copy();zr=np.zeros(self.n)
        jr=np.zeros((self.n,len(p))) if jac else None;jz=np.zeros((self.n,len(p))) if jac else None
        for s in self.spec:
            a=s['a'];off=s['off'];nd=s['nd'];nu=s['nu'];el=a['element']
            ph=s['phase']+s['hb']@p[off:off+nd]
            dw=np.exp(np.clip(s['hU']@p[off+nd:off+nd+nu],-60,40))
            cos=np.cos(ph);co=(cos*dw).sum(axis=1)*a['occ']
            ff=self.factors[el]+(fp if el=='Zr' else 0)
            real+=ff*co
            if el=='Zr':zr+=co
            if jac:
                dp=np.einsum('ns,nsj->nj',-np.sin(ph)*dw,s['hb'])*a['occ']
                du=np.einsum('ns,nst->nt',cos*dw,s['hU'])*a['occ']
                deriv=np.column_stack([dp,du]);jr[:,off:off+nd+nu]=ff[:,None]*deriv
                if el=='Zr':jz[:,off:off+nd+nu]=deriv
        for hs in self.hspec:
            s=hs['parent'];off=s['off'];nd=s['nd'];nu=s['nu']
            ph=hs['phase']+hs['hb']@p[off:off+nd];cos=np.cos(ph).sum(axis=1)
            U=hs['utrace']@p[off+nd:off+nd+nu];fac=hs['ff']*np.exp(-2*np.pi**2*self.q*U)
            term=fac*cos;real+=term
            if jac:
                jr[:,off:off+nd]+=fac[:,None]*np.einsum('ns,nsj->nj',-np.sin(ph),hs['hb'])
                jr[:,off+nd:off+nd+nu]+=term[:,None]*(-2*np.pi**2*self.q[:,None])*hs['utrace'][None,:]
        scale=np.exp(p[-1]);ic=scale**2*(real**2+t*zr**2)
        if details:return ic,real,zr,scale,fp,t
        if jac:
            if self.dispersion:jr[:,-3]=zr
            J=2*scale**2*(real[:,None]*jr+t*zr[:,None]*jz)
            if self.dispersion:J[:,-2]=scale**2*zr**2
            J[:,-1]=2*ic
            return ic,J
        return ic

    def adp_restraints(self,p,jac=False):
        res=[];rows=[]
        for s in self.spec:
            off=s['off']+s['nd'];nu=s['nu'];ub=s['ub']
            U=np.einsum('t,tij->ij',p[off:off+nu],ub)
            vals,vecs=np.linalg.eigh(U)
            for i,v in enumerate(vals):
                r=1000*min(v-.002,0);res.append(r)
                if jac:
                    row=np.zeros(len(p))
                    if v<.002:row[off:off+nu]=1000*np.einsum('i,tij,j->t',vecs[:,i],ub,vecs[:,i])
                    rows.append(row)
        if jac:return np.array(rows)
        return np.array(res)

    def update(self,p):
        out=[]
        for s in self.spec:
            a=dict(s['a']);off=s['off'];nd=s['nd'];nu=s['nu']
            a['xyz']=((np.array(a['xyz'])+np.array(a['basis'])@p[off:off+nd])%1).tolist()
            U=np.einsum('t,tij->ij',p[off+nd:off+nd+nu],s['ub']);a['Ucart']=U.tolist();a['Uiso']=float(np.trace(U)/3)
            a['U_eigenvalues']=np.linalg.eigvalsh(U).tolist();out.append(a)
        return out

    def update_hydrogens(self,p):
        out=[]
        for ah in self.hydrogens:
            ah=dict(ah);s=self.spec[ah['parent_index']];off=s['off'];nd=s['nd'];nu=s['nu']
            ah['xyz']=((np.array(ah['xyz'])+np.array(s['a']['basis'])@p[off:off+nd])%1).tolist()
            U=np.einsum('t,tij->ij',p[off+nd:off+nd+nu],s['ub'])
            ah['Uiso']=float(np.trace(U)*.4);ah['construction']='riding aromatic C-H 0.95 A; Uiso=1.2 Ueq(parent)'
            out.append(ah)
        return out

def refine_aniso(inp,tag,dispersion=True,dmin=1.,dmax=20.,nfev=100,aweight=.06):
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q'];mask=(q<=1/dmin**2)&(q>=1/dmax**2)
    m=AnisoModel(inp['atoms'],h[mask],dispersion);p=m.p0.copy();p[-1]=np.log(inp['scale'])
    if dispersion:p[-3]=inp.get('Zr_fp',-4);p[-2]=inp.get('Zr_fpp_squared',9)
    im=I[mask];sm=sig[mask];w=1/np.sqrt(sm**2+(aweight*np.maximum(im,0))**2+.01)
    def fun(p):return np.r_[(m.calc(p)-im)*w,m.adp_restraints(p)]
    def jac(p):
        _,j=m.calc(p,True)
        return np.vstack([j*w[:,None],m.adp_restraints(p,True)])
    print('ANISO',tag,'nref',len(im),'npar',len(p),'dispersion',dispersion,flush=True)
    opt=least_squares(fun,p,jac=jac,bounds=(m.lo,m.hi),x_scale='jac',max_nfev=nfev,ftol=2e-8)
    atoms=m.update(opt.x);mall=AnisoModel(atoms,h,dispersion);pa=mall.p0.copy();pa[-1]=opt.x[-1]
    if dispersion:pa[-3:]=opt.x[-3:]
    ic,real,zr,scale,fp,t=mall.calc(pa,details=True)
    report={'atoms':atoms,'scale':float(scale),'Zr_fp':float(fp),'Zr_fpp_squared':float(t),'dispersion_refined':dispersion,
            'n_params':len(p),'dmin':dmin,'dmax':dmax,'metrics_fit':metrics(I,sig,np.sqrt(ic),mask),
            'metrics_all':metrics(I,sig,np.sqrt(ic)),'cost':float(opt.cost),'success':bool(opt.success),'message':opt.message,'nfev':opt.nfev}
    (OUT/(tag+'.json')).write_text(json.dumps(report,indent=2));np.savez_compressed(OUT/(tag+'_fc.npz'),h=h,ic=ic,real=real,zr=zr,scale=scale,fp=fp,t=t)
    print('RESULT',report['metrics_fit'],'scale',scale,'fp',fp,'fpp2',t,'cost',opt.cost,flush=True)
    for a in atoms:print(a['label'],np.round(a['xyz'],6),'Ueig',np.round(a['U_eigenvalues'],5),flush=True)
    # Covariance of weighted LS without rescaling, stored for later uncertainty checks.
    cov=np.linalg.pinv(opt.jac.T@opt.jac,rcond=1e-10)
    np.savez_compressed(OUT/(tag+'_covariance.npz'),cov=cov,jac=opt.jac,parameters=opt.x,residual=opt.fun)
    amp=np.sqrt(real**2+t*zr**2);fo=np.sqrt(np.maximum(I,0))/scale
    # Real projection of conventional complex difference coefficients in centrosymmetry.
    delta=(fo/np.maximum(amp,1e-8)-1)*real
    for typ,co in [('diff',delta),('2fo',real+2*delta)]:
        rho=map_fft(h[mask],co[mask],shape=(240,240,120))
        np.save(OUT/(tag+'_'+typ+'_map.npy'),rho.astype(np.float32))
        save_peaks(peaks(rho,120,.7),tag+'_'+typ+'_peaks.csv')
    inspect_peaks(tag+'_diff_peaks.csv',atoms,25)
    return report

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='framework_iso1');ap.add_argument('--tag',default='framework_aniso1');ap.add_argument('--no-dispersion',action='store_true');ap.add_argument('--dmin',type=float,default=1.);ap.add_argument('--dmax',type=float,default=20.)
    args=ap.parse_args();inp=json.loads((OUT/(args.input+'.json')).read_text())
    refine_aniso(inp,args.tag,not args.no_dispersion,args.dmin,args.dmax)
