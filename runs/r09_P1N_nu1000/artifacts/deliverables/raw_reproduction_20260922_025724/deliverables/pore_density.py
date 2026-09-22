"""Experimental non-atomic pore-density estimation. No solvent species assigned.
Alternating measured-amplitude / positivity / support projections, low-pass limited.
Cross-validation flags exclude test reflections from density reconstruction.
"""
from crystal import *
from scipy.fft import fftn,ifftn

class PoreDensity:
    def __init__(self,atoms,h,q,shape=(192,192,96),cutoff=2.0,radius_scale=1.):
        self.shape=np.array(shape);self.h=h;self.q=q;self.cutoff=cutoff;self.radius_scale=radius_scale
        self.mask=np.ones(shape,np.float32)
        for a in atoms:
            radius={'C':2.8,'O':2.65,'Zr':3.1,'H':2.0}[a['element']]*radius_scale
            radii=np.ceil(radius*np.linalg.norm(np.linalg.inv(A),axis=1)*self.shape).astype(int)+1
            for x in expand_xyz(a['xyz']):
                center=np.rint(x*self.shape).astype(int)
                axes=[np.arange(-r,r+1)+c for r,c in zip(radii,center)]
                coords=np.array(np.meshgrid(*axes,indexing='ij')).reshape(3,-1).T
                xyz=coords/self.shape-x
                ds=np.linalg.norm(xyz@A.T,axis=1)
                # Cosine taper over 0.3 A outside excluded radii.
                vals=np.clip((ds-radius)/.3,0,1);vals=.5-.5*np.cos(np.pi*vals)
                ij=tuple((coords%self.shape).T)
                self.mask[ij]=np.minimum(self.mask[ij],vals)
        self.rho=np.zeros(shape,np.float32)
        eq=equiv_h(h).reshape(-1,3);idx=np.repeat(np.arange(len(h)),len(R))
        eh,ii=np.unique(eq,axis=0,return_index=True)
        self.fullh=eh;self.source=idx[ii];self.grid_index=tuple((eh%self.shape).T)
        self.h_index=tuple((h%self.shape).T)
        self.volfrac=float(np.mean(self.mask))
        self.gridq=None
        print('Pore mask volume fraction',self.volfrac,flush=True)

    def coefficients(self):
        grid=ifftn(self.rho,workers=1)*V
        return grid[self.h_index].real

    def iterate(self,framework,fo,sig,scale,work,niter=60,step=.8):
        use=work&(self.q<=1/self.cutoff**2)
        conf=np.clip(fo**2/(fo**2+sig),0,1)
        fs=self.coefficients();history=[]
        for it in range(niter):
            total=framework+fs
            coef=np.zeros(len(fo));coef[use]=(fo[use]/scale-np.abs(total[use]))*np.sign(total[use])*conf[use]
            grid=np.zeros(self.shape,complex);grid[self.grid_index]=coef[self.source]
            drho=fftn(grid,workers=1).real/V
            self.rho=np.maximum((self.rho+step*drho)*self.mask,0).astype(np.float32)
            fs=self.coefficients()
            if it%10==0 or it==niter-1:
                pred=np.abs(framework+fs)*scale
                err=float(np.sum(abs(fo[use]-pred[use]))/np.sum(fo[use]))
                test=(~work)&(self.q<=1/self.cutoff**2)
                te=float(np.sum(abs(fo[test]-pred[test]))/np.sum(fo[test]))
                rec={'iteration':it,'R_amp_work_low':err,'R_amp_test_low':te,'pore_electrons':float(self.rho.mean()*V),'rho_max':float(self.rho.max())};history.append(rec)
                print('PORE',rec,flush=True)
        return fs,history

if __name__=='__main__':
    from anisotropic import AnisoModel
    r=json.loads((OUT/'framework_iso1.json').read_text());d=np.load(OUT/'merged_equal.npz');h=d['h'];q=d['q'];I=d['I'];sig=d['sig']
    m=AnisoModel(r['atoms'],h,False);p=m.p0;p[-1]=np.log(r['scale']);ic,real,zr,scale,fp,t=m.calc(p,details=True)
    work=np.random.default_rng(712711).random(len(h))>.07
    pore=PoreDensity(r['atoms'],h,q)
    fs,hist=pore.iterate(real,np.sqrt(np.maximum(I,0)),sig,scale,work,niter=150)
    np.savez_compressed(OUT/'pore_initial.npz',rho=pore.rho,mask=pore.mask,fs=fs,work=work)
    (OUT/'pore_initial_history.json').write_text(json.dumps(hist,indent=2))
