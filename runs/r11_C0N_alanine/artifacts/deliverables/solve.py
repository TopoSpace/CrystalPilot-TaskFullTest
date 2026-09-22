"""Ab initio density modification in the experimentally supported P212121."""
from pathlib import Path
import itertools
import json
import time
import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter

OUT=Path(__file__).resolve().parent
SIGNS=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])
TRANS=np.array([[0,0,0],[.5,.5,0],[0,.5,.5],[.5,0,.5]])


def make_grid(df,shape=(32,32,64)):
    shape=np.array(shape)
    amp=np.zeros(shape)
    obs=np.zeros(shape,bool)
    # A smooth Wilson-style radial normalization removes atomic falloff.
    s2=1/(4*df.d.values**2)
    shell=pd.cut(s2,np.linspace(0,.55,16),labels=False)
    mx=[]
    for k in range(15):
        sel=shell==k
        if np.sum(sel)>10:
            mx.append([np.mean(s2[sel]),np.log(np.maximum(df.I.values[sel].mean(),1))])
    mx=np.array(mx)
    slope,intercept=np.polyfit(mx[:,0],mx[:,1],1)
    expected=np.exp(intercept+slope*s2)
    E=np.sqrt(np.maximum(df.I.values,0)/expected)
    E=np.minimum(E,4)
    for row,a in zip(df[["h","k","l"]].values.astype(int),E):
        for sg in itertools.product([-1,1],repeat=3):
            h=row*sg
            idx=tuple(h%shape)
            amp[idx]=a
            obs[idx]=True
    return amp,obs


def symmetry_indices(shape):
    grid=np.array(np.indices(shape)).reshape(3,-1).T
    return [np.ravel_multi_index(((grid*s+t*np.array(shape)).astype(int)%shape).T,shape)
            for s,t in zip(SIGNS,TRANS)]


def symmetrize(rho,indices):
    a=rho.ravel()
    return sum(a[i] for i in indices).reshape(rho.shape)/4


def peaks(rho,cell,n=6):
    shape=np.array(rho.shape)
    maxima=(rho==maximum_filter(rho,size=5,mode="wrap"))
    ij=np.argwhere(maxima)
    ij=ij[np.argsort(rho[tuple(ij.T)])[::-1]]
    chosen=[]
    vals=[]
    for i in ij:
        x=i/shape
        equiv=np.array([x*s+t for s,t in zip(SIGNS,TRANS)])%1
        if chosen:
            diff=equiv[:,None,:]-np.array(chosen)[None,:,:]
            diff-=np.rint(diff)
            if np.min(np.linalg.norm(diff*cell,axis=-1))<.85:
                continue
        chosen.append(x)
        vals.append(rho[tuple(i)])
        if len(chosen)==n:
            break
    return np.array(chosen),np.array(vals)


def model_resid(xyz,df,cell):
    h=df[["h","k","l"]].values
    atom=np.concatenate([xyz*s+t for s,t in zip(SIGNS,TRANS)])
    f=np.exp(2j*np.pi*h@atom.T).sum(axis=1)
    s2=1/(4*df.d.values**2)
    ic=abs(f)**2*np.exp(-4*s2)
    fo=np.sqrt(np.maximum(df.I.values,0))
    fc=np.sqrt(ic)
    scale=fo@fc/(fc@fc)
    return np.sum(abs(fo-scale*fc))/np.sum(fo)


if __name__=="__main__":
    source=OUT/"phasing_input.csv"
    df=pd.read_csv(source if source.exists() else OUT/"merged.csv")
    cell=np.array(json.loads((OUT/("geometry_flexible.json" if source.exists() else "geometry.json")).read_text())["cell"])
    amp,obs=make_grid(df)
    shape=amp.shape
    indices=symmetry_indices(shape)
    rng=np.random.default_rng(20260921)
    t=time.time()
    best=1.
    history=[]
    for trial in range(24):
        rho=symmetrize(rng.normal(size=shape),indices)
        f=np.fft.fftn(rho)
        f[obs]=amp[obs]*np.exp(1j*np.angle(f[obs]))
        f[~obs]=0
        f[0,0,0]=0
        rho=np.fft.ifftn(f).real
        delta=np.std(rho)*(.8+.05*(trial%5))
        for step in range(900):
            modified=np.where(rho<delta,-rho,rho)
            fc=np.fft.fftn(modified)
            error=np.sum(abs(abs(fc[obs])-amp[obs]))/amp[obs].sum()
            fc[obs]=amp[obs]*np.exp(1j*np.angle(fc[obs]))
            fc[~obs]=0
            fc[0,0,0]=0
            rho=np.fft.ifftn(fc).real
            if step%10==0:
                rho=symmetrize(rho,indices)
            if step%100==99:
                xyz,vals=peaks(rho,cell)
                score=model_resid(xyz,df,cell)
                history.append([trial,step,error,score])
                if score<best:
                    best=score
                    print("best",trial,step,"map",error,"model",score,"time",time.time()-t,flush=True)
                    np.savez(OUT/"solution_density.npz",rho=rho,xyz=xyz,values=vals,score=score)
                    (OUT/"solution_initial.json").write_text(json.dumps({"xyz":xyz.tolist(),
                        "peak_values":vals.tolist(),"score":float(score),"trial":trial,"iteration":step},indent=2))
                if best<.19:
                    break
        print("trial",trial,"best",best,flush=True)
        if best<.19:
            break
    np.savetxt(OUT/"phasing_history.csv",history,delimiter=",",
               header="trial,iteration,modulus_R,approximate_model_R",comments="")
