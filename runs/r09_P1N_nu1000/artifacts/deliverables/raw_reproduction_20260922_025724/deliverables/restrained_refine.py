from restraints import *
import argparse

def run(inp,tag,dispersion=False,dmin=1.,dmax=50.,epochs=4,aweight=.08,robust=False,solvent=None,all_data=False,hydrogens=False,geometry_strength=1.):
    data=np.load(OUT/'merged_equal.npz');h=data['h'];I=data['I'];sig=data['sig'];q=data['q']
    work=np.random.default_rng(712711).random(len(h))>.07
    fitrange=(q<=1/dmin**2)&(q>=1/dmax**2)
    fit=fitrange if all_data else fitrange&work
    seed=json.loads((OUT/'framework_iso1.json').read_text())['atoms']
    mall=AnisoModel(inp['atoms'],h,dispersion,solvent,hydrogens)
    m=AnisoModel(inp['atoms'],h[fit],dispersion,None if solvent is None else solvent[fit],hydrogens)
    geo=Geometry(m,seed,geometry_strength);simu=geo.adp_matrix()
    p=m.p0.copy();p[-1]=np.log(inp['scale'])
    if dispersion:p[-3]=inp.get('Zr_fp',-4);p[-2]=inp.get('Zr_fpp_squared',9)
    im=I[fit];sm=sig[fit];history=[]
    for epoch in range(epochs):
        ic=m.calc(p);P=(np.maximum(im,0)+2*ic)/3
        w=1/np.sqrt(sm**2+(aweight*P)**2+.01)
        def fun(p):return np.r_[(m.calc(p)-im)*w,geo.fun(p),simu@p,10*m.adp_restraints(p)]
        def jac(p):
            _,J=m.calc(p,True)
            return np.vstack([J*w[:,None],geo.jac(p),simu,10*m.adp_restraints(p,True)])
        print('RESTRAINED',tag,'epoch',epoch,'parameters',len(p),'fit',len(im),'restraints',len(fun(p))-len(im),flush=True)
        opt=least_squares(fun,p,jac=jac,bounds=(m.lo,m.hi),x_scale='jac',max_nfev=60,ftol=5e-8,loss='soft_l1' if robust else 'linear',f_scale=2.)
        p=opt.x;icfull=mall.calc(p)
        hist={'epoch':epoch,'work':metrics(I,sig,np.sqrt(icfull),fitrange&work),'test':metrics(I,sig,np.sqrt(icfull),fitrange&~work),'cost':float(opt.cost),'max_geom_z':float(max(abs(geo.fun(p)))),'nfev':int(opt.nfev)}
        history.append(hist);print(hist,flush=True)
    atoms=m.update(p);ic,real,zr,scale,fp,t=mall.calc(p,details=True)
    report={'atoms':atoms,'riding_hydrogens':m.update_hydrogens(p),'scale':float(scale),'Zr_fp':float(fp),'Zr_fpp_squared':float(t),'dispersion_refined':dispersion,'n_params':len(p),
            'dmin':dmin,'dmax':dmax,'all_data':all_data,'weight_a':aweight,'robust':robust,'hydrogens_included':hydrogens,'geometry_strength':geometry_strength,'metrics_fit':metrics(I,sig,np.sqrt(ic),fit),
            'metrics_work':metrics(I,sig,np.sqrt(ic),fitrange&work),'metrics_test':metrics(I,sig,np.sqrt(ic),fitrange&~work),
            'metrics_all_range':metrics(I,sig,np.sqrt(ic),fitrange),'metrics_all_recorded':metrics(I,sig,np.sqrt(ic)),
            'n_restraints':len(fun(p))-len(im),'geometry_restraint_rms_z':float(np.sqrt(np.mean(geo.fun(p)**2))),
            'history':history,'cost':float(opt.cost),'success':bool(opt.success),'message':opt.message}
    (OUT/(tag+'.json')).write_text(json.dumps(report,indent=2));np.savez_compressed(OUT/(tag+'_fc.npz'),h=h,ic=ic,real=real,zr=zr,scale=scale,fp=fp,t=t,work=work,fit=fit)
    cov=np.linalg.pinv(opt.jac.T@opt.jac,rcond=1e-10)
    np.savez_compressed(OUT/(tag+'_covariance.npz'),cov=cov,jac=opt.jac,parameters=p,residual=opt.fun)
    print('FINAL',tag,report['metrics_all_range'],'test',report['metrics_test'],'scale',scale,'fp',fp,'fpp2',t,flush=True)
    for a in atoms:print(a['label'],np.round(a['xyz'],6),'Ueig',np.round(a['U_eigenvalues'],5),flush=True)
    mask=fitrange;fo=np.sqrt(np.maximum(I,0))/scale;amp=np.sqrt(real**2+t*zr**2);delta=(fo/np.maximum(amp,1e-8)-1)*real
    rho=map_fft(h[mask],delta[mask],shape=(240,240,120));np.save(OUT/(tag+'_diff_map.npy'),rho.astype(np.float32));save_peaks(peaks(rho,120,.7),tag+'_diff_peaks.csv')
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='framework_iso1');ap.add_argument('--tag',default='framework_restrained');ap.add_argument('--dispersion',action='store_true');ap.add_argument('--robust',action='store_true');ap.add_argument('--dmin',type=float,default=1.);ap.add_argument('--dmax',type=float,default=50.);ap.add_argument('--epochs',type=int,default=4);ap.add_argument('--all-data',action='store_true');ap.add_argument('--hydrogens',action='store_true');ap.add_argument('--geometry-strength',type=float,default=1.)
    a=ap.parse_args();inp=json.loads((OUT/(a.input+'.json')).read_text())
    run(inp,a.tag,a.dispersion,a.dmin,a.dmax,a.epochs,robust=a.robust,all_data=a.all_data,hydrogens=a.hydrogens,geometry_strength=a.geometry_strength)
