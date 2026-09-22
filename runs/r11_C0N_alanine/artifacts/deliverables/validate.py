"""Independent numerical checks, conditional cross-validation, and geometry."""
from pathlib import Path
import hashlib
import json
import platform
import struct
import sys
import itertools
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares
from rawio import read_frame
from geometry import WAVELENGTH
from geometry_flexible import unpack,q_from_peaks
from refine import Model
from refine_aniso import AnisoModel,optimize,hydrogens,u_matrices,LABELS,H_LABELS
from solve import SIGNS,TRANS

OUT=Path(__file__).resolve().parent


def frame_checks():
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    errors=[]
    for rec in manifest:
        p=OUT.parent/f"inputs/alanine/frames/pgw240033_Mo_{rec['run']}_{rec['frame']}.rodhypix"
        b=p.read_bytes()
        assert hashlib.sha256(b).hexdigest()==rec["sha256"]
        errors.append(abs(struct.unpack_from("<d",b,2584)[0]*620000-rec["total"]))
    tests=[(1,1),(2,55),(3,73),(4,25),(5,80),(6,66)]
    for run,frame in tests:
        a,_=read_frame(OUT.parent/f"inputs/alanine/frames/pgw240033_Mo_{run}_{frame}.rodhypix")
        assert np.array_equal(a,np.load(OUT/f"cache/{run}_{frame}.npy"))
    return {"formal_frames":len(manifest),"all_raw_sha256_unchanged":True,
            "max_header_mean_total_discrepancy_counts":float(max(errors)),
            "total_integer_counts":sum(v["total"] for v in manifest),
            "maximum_pixel_count":max(v["maximum"] for v in manifest),
            "decoded_cache_roundtrip_tests":len(tests),"preexperiment_frames_not_used":30}


def geometry_jackknife():
    geom=json.loads((OUT/"geometry_flexible.json").read_text())
    z0=np.array(geom["parameter_vector"])
    g=np.loadtxt(OUT/"geometry_flexible_residuals.csv",delimiter=",")
    cells=[]
    residuals=[]
    for run in range(1,7):
        gg=g[g[:,0]!=run]
        h=gg[:,10:13]
        w=np.minimum(np.sqrt(gg[:,4]/300),2)
        def fun(z):
            return ((q_from_peaks(gg,z)-h@unpack(z)[0].T)*w[:,None]).ravel()
        span=np.r_[np.full(3,.03),np.full(3,.02),np.full(6,.7),np.full(6,.02),.5,.5,np.full(5,.3)]
        fit=least_squares(fun,z0,bounds=(z0-span,z0+span),loss="soft_l1",f_scale=.0003,
                          x_scale="jac",max_nfev=80,ftol=1e-10,gtol=1e-10,xtol=1e-10)
        cells.append(np.exp(fit.x[3:6]))
        test=g[g[:,0]==run]
        residuals.append(float(np.sqrt(np.mean((q_from_peaks(test,fit.x)-test[:,10:13]@unpack(fit.x)[0].T)**2))))
    cells=np.array(cells)
    jack=np.std(cells,axis=0)*np.sqrt(5)
    used=np.maximum(jack,geom["cell_esd"])
    return {"leave_one_run_out_cells":cells.tolist(),"jackknife_cell_esd":jack.tolist(),
            "conditional_cell_esd":geom["cell_esd"],"adopted_internal_cell_esd":used.tolist(),
            "heldout_q_rms_conditional_on_scan_zero":residuals,
            "warning":"Scan-zero initialization and fixed pixel pitch are shared. Not an external calibration."}


def completeness(cell):
    df=pd.read_csv(OUT/"merged_corrected.csv")
    h=np.array(np.meshgrid(np.arange(12),np.arange(12),np.arange(24),indexing="ij")).reshape(3,-1).T
    q=np.linalg.norm(h/cell,axis=1)
    allowed=(q>0)&((np.count_nonzero(h,axis=1)>1)|(h.sum(axis=1)%2==0))
    ds=np.divide(1,q,out=np.full_like(q,np.inf),where=q>0)
    rows=[]
    for lim in [.70,.73,.75,.80,.90,1.,1.2,1.5,2.]:
        expected=int(np.sum(allowed&(ds>=lim)))
        measured=int(np.sum(df.d>=lim))
        rows.append({"d_min":lim,"expected":expected,"measured":measured,"fraction":measured/expected})
    return rows


def crossvalidate(model):
    df=pd.read_csv(OUT/"merged_final.csv")
    z0=np.array(model["parameter_vector"])
    cell=np.array(model["cell"])
    rng=np.random.default_rng(170723)
    permutation=rng.permutation(len(df))
    folds=np.empty(len(df),int)
    folds[permutation]=np.arange(len(df))%5
    predicted=np.empty(len(df))
    metrics=[]
    for fold in range(5):
        train=df[folds!=fold]
        test=df[folds==fold]
        z,_,fit,sigma=optimize(train,cell,z0.copy(),cycles=2)
        fc=AnisoModel(test[["h","k","l"]].values,cell).calculate(z)[1]
        predicted[folds==fold]=fc
        fo=np.sqrt(np.maximum(test.I,0))
        r=float(np.sum(abs(fo-np.sqrt(fc)))/np.sum(fo))
        metrics.append({"fold":fold,"train_n":len(train),"test_n":len(test),"R1_test":r})
    fo=np.sqrt(np.maximum(df.I,0))
    robs=df.I>2*df.sigma
    total=float(sum(abs(fo-np.sqrt(predicted)))/sum(fo))
    observed=float(sum(abs(fo[robs]-np.sqrt(predicted[robs])))/sum(fo[robs]))
    result=df.copy()
    result["cv_fold"]=folds
    result["heldout_Fc2"]=predicted
    result.to_csv(OUT/"cross_validation_reflections.csv",index=False)
    return {"folds":metrics,"combined_R1_test_all":total,"combined_R1_test_observed":observed,
            "qualification":"Conditional parameter-refinement validation; initial phasing and SG used all data."}


def geometry_tables(model,cell_esd):
    z=np.array(model["parameter_vector"])
    cell=np.array(model["cell"])
    cov=np.load(OUT/"final_parameter_covariance.npy")
    ccov=np.diag(np.asarray(cell_esd)**2)
    combined=np.zeros((60,60))
    combined[:57,:57]=cov
    combined[57:,57:]=ccov
    zz=np.r_[z,cell]
    def uncertainty(fn):
        grad=np.zeros(60)
        for i in range(60):
            step=1e-5 if i>=57 else 1e-6
            x1,x2=zz.copy(),zz.copy()
            x1[i]+=step
            x2[i]-=step
            grad[i]=(fn(x1)-fn(x2))/(2*step)
        return float(np.sqrt(max(grad@combined@grad,0)))
    bonds=[]
    for i,j in [(0,3),(1,3),(3,4),(2,4),(4,5)]:
        fn=lambda a:np.linalg.norm((a[:18].reshape(6,3)[i]-a[:18].reshape(6,3)[j])*a[57:])
        bonds.append({"atom1":LABELS[i],"atom2":LABELS[j],"distance":float(fn(zz)),"esd":uncertainty(fn)})
    angles=[]
    for i,j,k in [(0,3,1),(0,3,4),(1,3,4),(2,4,3),(2,4,5),(3,4,5)]:
        def fn(a):
            x=a[:18].reshape(6,3)*a[57:]
            v,w=x[i]-x[j],x[k]-x[j]
            return np.degrees(np.arccos(np.clip(v@w/(np.linalg.norm(v)*np.linalg.norm(w)),-1,1)))
        angles.append({"atom1":LABELS[i],"atom2":LABELS[j],"atom3":LABELS[k],
                       "angle":float(fn(zz)),"esd":uncertainty(fn)})
    xyz=z[:18].reshape(6,3)
    hx=hydrogens(xyz,cell,z[55:57])
    hbonds=[]
    for hi in [1,2,3]:
        donor=xyz[2]
        hydrogen=hx[hi]
        for oi in [0,1]:
            for op,(sg,tr) in enumerate(zip(SIGNS,TRANS)):
                for shift in itertools.product([-2,-1,0,1,2],repeat=3):
                    acceptor=xyz[oi]*sg+tr+shift
                    ha=np.linalg.norm((acceptor-hydrogen)*cell)
                    da=np.linalg.norm((acceptor-donor)*cell)
                    v,w=(donor-hydrogen)*cell,(acceptor-hydrogen)*cell
                    ang=np.degrees(np.arccos(np.clip(v@w/(np.linalg.norm(v)*np.linalg.norm(w)),-1,1)))
                    if ha<2.6 and da<3.3 and ang>120:
                        hbonds.append({"donor":"N1","hydrogen":H_LABELS[hi],"acceptor":LABELS[oi],
                                       "symop":op+1,"translation":list(shift),"H_A":float(ha),
                                       "D_A":float(da),"D_H_A":float(ang)})
    pd.DataFrame(bonds).to_csv(OUT/"bond_lengths.csv",index=False)
    pd.DataFrame(angles).to_csv(OUT/"bond_angles.csv",index=False)
    pd.DataFrame(hbonds).to_csv(OUT/"hydrogen_bonds.csv",index=False)
    return {"bonds":bonds,"angles":angles,"hydrogen_bonds":hbonds}


def numerical_checks(model):
    z=np.array(model["parameter_vector"])
    cell=np.array(model["cell"])
    h=np.array([[1,2,3],[2,-3,1],[0,2,1],[2,0,0],[0,1,0],[3,0,0],[0,0,5]])
    m=AnisoModel(h,cell)
    f,ic=m.calculate(z)
    fm=AnisoModel(-h,cell).calculate(z)[0]
    assert np.max(abs(f-fm.conj()))<1e-10
    assert np.max(abs(f[4:]))<1e-10
    eig=np.linalg.eigvalsh(u_matrices(z[18:54].reshape(6,6)))
    assert eig.min()>0
    iso=Model(h,cell,model["types"])
    zi=np.r_[z[:18],z[18:54].reshape(6,6)[:,:3].mean(axis=1),z[54]]
    _,_,analytic=iso.calculate(zi)
    numeric=np.empty_like(analytic)
    for i in range(len(zi)):
        za,zb=zi.copy(),zi.copy()
        za[i]+=1e-7
        zb[i]-=1e-7
        numeric[:,i]=(iso.calculate(za,False)[1]-iso.calculate(zb,False)[1])/2e-7
    # Exclude exactly absent reflections, where the derivative of |F| is undefined.
    err=np.max(abs(numeric[:4]-analytic[:4]))
    assert err<2e-4
    return {"Friedel_conjugacy_max_error":float(max(abs(f-fm.conj()))),
            "screw_absence_max_amplitude":float(max(abs(f[4:]))),
            "isotropic_amplitude_jacobian_max_error":float(err),
            "minimum_U_eigenvalue":float(eig.min())}


if __name__=="__main__":
    model=json.loads((OUT/"model_final.json").read_text())
    print("Validating raw frames",flush=True)
    raw=frame_checks()
    print(raw,flush=True)
    print("Geometry leave-one-run-out",flush=True)
    gj=geometry_jackknife()
    print(gj,flush=True)
    print("Five-fold refinement validation",flush=True)
    cv=crossvalidate(model)
    report={"runtime":{"python":sys.version,"executable":sys.executable,"platform":platform.platform(),
                       "numpy":np.__version__,"scipy":scipy.__version__},
            "raw_data":raw,"geometry_stability":gj,"completeness":completeness(np.array(model["cell"])),
            "cross_validation":cv,"numerical_checks":numerical_checks(model),
            "molecular_geometry":geometry_tables(model,gj["adopted_internal_cell_esd"])}
    (OUT/"validation.json").write_text(json.dumps(report,indent=2))
    print("Validation complete",json.dumps(cv),flush=True)
