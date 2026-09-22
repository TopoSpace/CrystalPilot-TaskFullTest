"""Prediction, 3-D summation integration, local background, rotation LP."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from geometry import *

OUT=Path(__file__).resolve().parent


def predictions(ub,params,manifest,dmin=.70):
    cell=WAVELENGTH*np.linalg.norm(np.linalg.inv(ub).T,axis=0)
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
        r=gonio(0,kappa+params[5],phi,-1,-1,-1,0,-params[4])
        v=q@r.T
        rho=np.linalg.norm(v[:,:2],axis=1)
        c=-qn**2/(2*rho)
        ok=np.abs(c)<1
        delta=np.arctan2(v[:,1],v[:,0])
        for sign in [-1,1]:
            om=np.degrees(delta+sign*np.arccos(np.clip(c,-1,1)))-params[6]
            om=(om+180)%360-180
            take=ok&(om>amin+.8)&(om<amax-.8)
            hh,vv,oo=h[take],v[take],om[take]
            qq=np.einsum("...ij,...j->...i",rz(-oo-params[6]),vv)
            ss=qq+[1,0,0]
            for mod,pr in enumerate(project(ss,theta,params)):
                x,y=pr[:,:2].T
                xmin,xmax=[(9,375),(424,790)][mod]
                good=(x>xmin)&(x<xmax)&(y>9)&(y<765)&(pr[:,2]>0)&(np.abs(qq[:,1])>.025)
                for hi,oi,xi,yi,qi,si in zip(hh[good],oo[good],x[good],y[good],qq[good],ss[good]):
                    pol=.5*(1+si[0]**2)
                    lp=abs(qi[1])/pol
                    allrows.append([run,*hi,oi,xi,yi,mod,lp,WAVELENGTH/np.linalg.norm(qi)])
    return np.array(allrows)


def integrate(pred,centroids,manifest,radius=4,rocking=.85,adaptive=False):
    lookup={tuple(map(int,g[[0,10,11,12]])):g for g in centroids}
    images={}
    for s in manifest:
        images[(s["run"],s["frame"])]=np.load(OUT/f"cache/{s['run']}_{s['frame']}.npy",mmap_mode="r")
    rows=[]
    profiles=[]
    for n,row in enumerate(pred):
        run,h,k,l,omega,x,y,mod,lp,d=row
        run=int(run)
        key=tuple(map(int,row[:4]))
        center=lookup.get(key)
        dx=dy=dw=0.
        if center is not None:
            dx,dy,dw=center[2]-x,center[3]-y,center[6]-omega
            if abs(dx)>3 or abs(dy)>3 or abs(dw)>.5:
                dx=dy=dw=0
        x+=dx
        y+=dy
        omega+=dw
        imin,imax=int(round(x))-8,int(round(x))+9
        jmin,jmax=int(round(y))-8,int(round(y))+9
        yy,xx=np.mgrid[jmin:jmax,imin:imax]
        r2=(xx-x)**2+(yy-y)**2
        signal=r2<=radius**2
        back=(r2>=36)&(r2<=64)
        ns,nb=signal.sum(),back.sum()
        total=variance=background=0.
        nfr=0
        width=rocking
        if adaptive:
            qnorm=WAVELENGTH/d
            pol=.5*(1+(1-qnorm*qnorm/2)**2)
            width=np.clip(.35+.18*qnorm/(lp*pol),.65,2.8)
        for s in manifest:
            if s["run"]!=run:
                continue
            om0,om1=s["angles"][0][0],s["angles"][1][0]
            if om1<omega-width or om0>omega+width:
                continue
            img=images[(run,s["frame"])]
            patch=img[jmin:jmax,imin:imax].astype(float)
            bc=patch[back]
            # Reject unrelated Bragg pixels from the background only.
            mean0=bc.mean()
            clean=bc[bc<max(4,mean0+5*np.sqrt(mean0+1))]
            bmean=clean.mean()
            raw=patch[signal].sum()
            value=raw-ns*bmean
            var=raw+ns**2*max(bmean,.02)/len(clean)
            total+=value
            variance+=var
            background+=ns*bmean
            nfr+=1
            profiles.append([n,run,s["frame"],.5*(om0+om1),raw,bmean,value,var])
        rows.append([*row,total,np.sqrt(max(variance,1)),total*lp,np.sqrt(max(variance,1))*lp,
                     background,nfr,dx,dy,dw])
        if n%500==0:
            print("integrated",n,len(pred),flush=True)
    cols=["run","h","k","l","omega","x","y","module","lp","d","counts","sigma_counts",
          "I","sigma","background","nframes","dx","dy","domega"]
    return pd.DataFrame(rows,columns=cols),profiles


if __name__=="__main__":
    geom=json.loads((OUT/"geometry.json").read_text())
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    pred=predictions(np.array(geom["ub"]),geom["params"],manifest)
    np.savetxt(OUT/"predictions.csv",pred,delimiter=",",
               header="run,h,k,l,omega,x,y,module,lp,d",comments="")
    print("predictions",len(pred))
    centroids=np.loadtxt(OUT/"indexed_centroids.csv",delimiter=",",skiprows=1)
    df,profiles=integrate(pred,centroids,manifest)
    df.to_csv(OUT/"unmerged.csv",index=False)
    np.save(OUT/"rocking_profiles.npy",np.array(profiles))
    print("positive >3sig",np.sum(df.I>3*df.sigma))
