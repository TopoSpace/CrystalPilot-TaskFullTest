"""Refine independent detector module planes and a common orthorhombic metric."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.optimize import least_squares
from geometry import WAVELENGTH,rz,gonio

OUT=Path(__file__).resolve().parent


def unpack(z):
    ub=Rotation.from_rotvec(z[:3]).as_matrix()@np.diag(WAVELENGTH/np.exp(z[3:6]))
    origins=z[6:12].reshape(2,3)
    mats=Rotation.from_rotvec(z[12:18].reshape(2,3)).as_matrix()
    return ub,origins,mats


def q_from_peaks(g,z):
    ub,origins,mats=unpack(z)
    mod=(g[:,2]>400).astype(int)
    mx=np.where(mod,605.968,192.5)
    v=origins[mod]-(g[:,2]-mx)[:,None]*.1*mats[mod,:,1]-(g[:,3]-387.5)[:,None]*.1*mats[mod,:,2]
    v=np.einsum("hij,hj->hi",rz(-g[:,7]),v)
    s=v/np.linalg.norm(v,axis=1)[:,None]
    omega=g[:,6]-.36049+np.r_[0,z[20:25]][g[:,0].astype(int)-1]
    r=gonio(omega,g[:,8]+z[18],g[:,9],-1,-1,-1,0,-z[19])
    return np.einsum("hji,hj->hi",r,s-[1,0,0])


def project_flexible(s,theta,z):
    ub,origins,mats=unpack(z)
    v=np.asarray(s)@rz(theta).T
    result=[]
    for mod in [0,1]:
        origin=origins[mod]
        r=mats[mod]
        n=r[:,0]
        t=(origin@n)/(v@n)
        local=(v*t[...,None]-origin)@r
        x=[192.5,605.968][mod]-local[...,1]/.1
        y=387.5-local[...,2]/.1
        result.append(np.stack([x,y,t],axis=-1))
    return result


def initial_z(geom):
    ub=np.array(geom["ub"])
    u,s,v=np.linalg.svd(ub@np.diag(np.array(geom["cell"])/WAVELENGTH))
    rotation=u@v
    if np.linalg.det(rotation)<0:
        raise ValueError("Left-handed basis")
    rv=Rotation.from_matrix(rotation).as_rotvec()
    dd,cx,cy,tzero,alpha,kzero,ozero=geom["params"]
    origins=[]
    mats=[]
    for mod in [0,1]:
        r=rz([19,-19.016][mod])
        org=np.array([60,59.857][mod]*r[:,0])+[dd-60,(cx-400)*.1,(cy-387.5)*.1]
        org=rz(-tzero)@org
        r=rz(-tzero)@r
        origins.append(org)
        mats.append(Rotation.from_matrix(r).as_rotvec())
    return np.r_[rv,np.log(geom["cell"]),np.ravel(origins),np.ravel(mats),kzero,alpha,np.zeros(5)]


def run():
    geom=json.loads((OUT/"geometry.json").read_text())
    z0=initial_z(geom)
    g=np.loadtxt(OUT/"indexed_centroids.csv",delimiter=",",skiprows=1)
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    runlim={r:(min(f["angles"][0][0] for f in manifest if f["run"]==r),
               max(f["angles"][1][0] for f in manifest if f["run"]==r)) for r in range(1,7)}
    use=np.array([(row[6]>runlim[int(row[0])][0]+1.5 and
                   row[6]<runlim[int(row[0])][1]-1.5) for row in g])
    use&=g[:,4]>300
    g=g[use]
    h=g[:,10:13]
    weights=np.minimum(np.sqrt(g[:,4]/300),2)
    def res(z):
        ub=unpack(z)[0]
        return ((q_from_peaks(g,z)-h@ub.T)*weights[:,None]).ravel()
    lo=z0-np.r_[np.full(3,.05),np.full(3,.03),np.full(6,1.5),np.full(6,.03),1,1,np.full(5,.5)]
    hi=z0+np.r_[np.full(3,.05),np.full(3,.03),np.full(6,1.5),np.full(6,.03),1,1,np.full(5,.5)]
    fit=least_squares(res,z0,bounds=(lo,hi),loss="soft_l1",f_scale=.0003,
                      max_nfev=200,x_scale="jac",ftol=1e-11,xtol=1e-11,gtol=1e-11)
    z=fit.x
    ub,origins,mats=unpack(z)
    errors=q_from_peaks(g,z)-h@ub.T
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*np.sum(fit.fun**2)/(len(fit.fun)-len(z))
    esd=np.sqrt(np.diag(cov))
    cell=np.exp(z[3:6])
    print("N",len(g),"initial rms",np.sqrt(np.mean(res(z0)**2)),
          "final rms",np.sqrt(np.mean(errors**2)),"evaluations",fit.nfev,flush=True)
    print("cell",cell,"conditional ESD",cell*esd[3:6],"\norigins",origins,
          "\nrotations(deg)",np.rad2deg(z[12:18].reshape(2,3)),
          "\nkzero alpha",z[18:20],"\nrun omega",z[20:],flush=True)
    result=dict(parameter_vector=z.tolist(),parameter_esd=esd.tolist(),
                ub=ub.tolist(),cell=cell.tolist(),cell_esd=(cell*esd[3:6]).tolist(),
                origins=origins.tolist(),rotations=mats.tolist(),
                residual_q_rms=float(np.sqrt(np.mean(errors**2))),n_centroids=len(g),
                note="Conditional least-squares standard errors, not external calibration accuracy.")
    (OUT/"geometry_flexible.json").write_text(json.dumps(result,indent=2))
    np.savetxt(OUT/"geometry_flexible_residuals.csv",np.c_[g,errors],delimiter=",")
    return result


if __name__=="__main__":
    run()
