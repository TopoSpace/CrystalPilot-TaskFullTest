from pathlib import Path
import json
import numpy as np
from geometry import WAVELENGTH,rz,gonio
from geometry_flexible import unpack,project_flexible
from integrate import integrate
from shadow_mask import apply
from merge import merge

OUT=Path(__file__).resolve().parent


def predict(z,manifest,dmin=.70):
    ub=unpack(z)[0]
    cell=np.exp(z[3:6])
    limits=np.ceil(cell/dmin).astype(int)
    h=np.array(np.meshgrid(*[np.arange(-n,n+1) for n in limits],indexing="ij")).reshape(3,-1).T
    q=h@ub.T
    qnorm=np.linalg.norm(q,axis=1)
    use=(qnorm>0)&(qnorm<WAVELENGTH/dmin)
    h,q,qn=h[use],q[use],qnorm[use]
    allrows=[]
    for run in range(1,7):
        frames=[s for s in manifest if s["run"]==run]
        amin=min(s["angles"][0][0] for s in frames)
        amax=max(s["angles"][1][0] for s in frames)
        _,theta,kappa,phi,_=frames[0]["angles"][0]
        r=gonio(0,kappa+z[18],phi,-1,-1,-1,0,-z[19])
        v=q@r.T
        rho=np.linalg.norm(v[:,:2],axis=1)
        c=-qn**2/(2*rho)
        ok=np.abs(c)<1
        delta=np.arctan2(v[:,1],v[:,0])
        ozero=-.36049+np.r_[0,z[20:25]][run-1]
        for sign in [-1,1]:
            om=np.degrees(delta+sign*np.arccos(np.clip(c,-1,1)))-ozero
            om=(om+180)%360-180
            take=ok&(om>amin+.7)&(om<amax-.7)
            hh,vv,oo=h[take],v[take],om[take]
            qq=np.einsum("...ij,...j->...i",rz(-oo-ozero),vv)
            ss=qq+[1,0,0]
            for mod,pr in enumerate(project_flexible(ss,theta,z)):
                x,y=pr[:,:2].T
                xmin,xmax=[(9,375),(424,790)][mod]
                good=(x>xmin)&(x<xmax)&(y>9)&(y<765)&(pr[:,2]>0)&(np.abs(qq[:,1])>.025)
                for hi,oi,xi,yi,qi,si in zip(hh[good],oo[good],x[good],y[good],qq[good],ss[good]):
                    pol=.5*(1+si[0]**2)
                    lp=abs(qi[1])/pol
                    width=np.clip(.35+.18*np.linalg.norm(qi)/abs(qi[1]),.65,2.8)
                    if oi<amin+width+.15 or oi>amax-width-.15:
                        continue
                    allrows.append([run,*hi,oi,xi,yi,mod,lp,WAVELENGTH/np.linalg.norm(qi)])
    return np.array(allrows)


if __name__=="__main__":
    geom=json.loads((OUT/"geometry_flexible.json").read_text())
    z=np.array(geom["parameter_vector"])
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    pred=predict(z,manifest)
    np.savetxt(OUT/"predictions_flexible.csv",pred,delimiter=",",
               header="run,h,k,l,omega,x,y,module,lp,d",comments="")
    print("predictions",len(pred),flush=True)
    centroids=np.loadtxt(OUT/"indexed_centroids.csv",delimiter=",",skiprows=1)
    df,profiles=integrate(pred,centroids,manifest,adaptive=True)
    df.to_csv(OUT/"unmerged_flexible.csv",index=False)
    np.save(OUT/"rocking_profiles_flexible.npy",np.array(profiles))
    with np.load(OUT/"background_transmission.npz") as archive:
        masks={k:archive[k] for k in archive.files}
    apply(df,masks)
    m,obs=merge("unmerged_masked.csv","merged_flexible")
    selected=m[(~m.absent_212121)&(m.d>=.75)]
    selected.to_csv(OUT/"merged_final.csv",index=False)
    print("selected final",len(selected),"of",len(m),flush=True)
