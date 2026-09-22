"""Independent centrosymmetric charge-flipping phase-extension diagnostic.
No atomic or reference structure library. Experimental; not automatically accepted.
"""
from crystal import *
from scipy.fft import fftn,ifftn
from anisotropic import AnisoModel
import argparse

def run(niter=800):
    d=np.load(OUT/'merged_equal.npz');h=d['h'];I=d['I'];sig=d['sig'];q=d['q']
    r=json.loads((OUT/'framework_h_highangle.json').read_text());m=AnisoModel(r['atoms'],h,True,hydrogens=True);p=m.p0;p[-3]=r['Zr_fp'];p[-2]=r['Zr_fpp_squared'];p[-1]=np.log(r['scale'])
    ic,fr,zr,scale,fp,t=m.calc(p,details=True)
    # For phase extension subtract the estimated imaginary Zr contribution.
    ie=I/scale**2-t*zr**2;se=sig/scale**2
    fo=np.sqrt(np.maximum(ie,0));flo=np.sqrt(np.maximum(ie-se,0));fhi=np.sqrt(np.maximum(ie+se,0))
    work=np.random.default_rng(712711).random(len(h))>.07
    used=(q<=1)&work;test=(q<=1)&~work;shape=np.array([128,128,64])
    eq=equiv_h(h[q<=1]).reshape(-1,3);source=np.repeat(np.where(q<=1)[0],len(R));eh,jj=np.unique(eq,axis=0,return_index=True);source=source[jj];ij=tuple((eh%shape).T)
    usefull=used[source];testfull=test[source];hj=tuple((h%shape).T)
    fgrid=np.zeros(shape,complex);fgrid[ij]=fr[source]
    rng=np.random.default_rng(62891);history=[]
    for trial,delta_ratio in enumerate([.35,.6,.85]):
        F=fgrid.copy();F[ij]=np.where(usefull,fo[source]*np.sign(fr[source]),F[ij]);F[0,0,0]=0
        best=None
        for it in range(niter):
            rho=fftn(F,workers=1).real/V
            delta=delta_ratio*np.std(rho)
            flipped=np.where(rho<delta,-rho,rho)
            G=ifftn(flipped,workers=1)*V
            vals=G[ij].real
            amp=np.abs(vals)
            # Project amplitudes to their one-sigma uncertainty intervals.
            amp[usefull]=np.clip(amp[usefull],flo[source[usefull]],fhi[source[usefull]])
            newvals=amp*np.sign(vals)
            F.fill(0);F[ij]=newvals;F[0,0,0]=G[0,0,0].real
            if it%40==0 or it==niter-1:
                pred=G[hj].real
                obs=(I>2*sig)&(q<=1)
                fco=np.sqrt(scale**2*(pred**2+t*zr**2));fo0=np.sqrt(np.maximum(I,0))
                rw=float(np.sum(abs(fco[obs&work]-fo0[obs&work]))/np.sum(fo0[obs&work]));rt=float(np.sum(abs(fco[obs&~work]-fo0[obs&~work]))/np.sum(fo0[obs&~work]))
                changed=int(np.sum(np.sign(F[hj][used].real)!=np.sign(fr[used])))
                rec={'trial':trial,'delta_ratio':delta_ratio,'iteration':it,'R1_work_modified_density':rw,'R1_test_modified_density':rt,'phases_changed':changed,'DC':float(F[0,0,0].real)};history.append(rec);print('FLIP',rec,flush=True)
                if best is None or rw<best['score']:
                    best={'score':rw,'F':F.copy(),'G':G.copy(),'record':rec}
        np.savez_compressed(OUT/f'chargeflip_trial{trial}.npz',F=best['F'].astype(np.complex64),G=best['G'].astype(np.complex64),shape=shape,h=h,work=work)
        rho=fftn(best['F']).real/V;save_peaks(peaks(rho,100,.7),f'chargeflip_trial{trial}_peaks.csv')
    (OUT/'chargeflip_history.json').write_text(json.dumps(history,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--iterations',type=int,default=800);x=a.parse_args();run(x.iterations)
