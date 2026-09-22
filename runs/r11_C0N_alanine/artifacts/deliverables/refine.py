"""Self-contained kinematic X-ray structure factors and least-squares refinement."""
from pathlib import Path
import itertools
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from solve import SIGNS,TRANS

OUT=Path(__file__).resolve().parent
# Conventional neutral-atom four-Gaussian coefficients, explicitly recorded.
# Numerical constants are not taken from any local crystallographic package.
SCATTER={
    "C":([2.31,1.02,1.5886,.865],[20.8439,10.2075,.5687,51.6512],.2156),
    "N":([12.2126,3.1322,2.0125,1.1663],[.0057,9.8933,28.9975,.5826],-11.529),
    "O":([3.0485,2.2868,1.5463,.867],[13.2771,5.7011,.3239,32.9089],.2508),
    "H":([.493002,.322912,.140191,.040810],[10.5109,26.1257,3.14236,57.7997],.003038)
}


def ff(element,s2):
    a,b,c=SCATTER[element]
    return (np.array(a)*np.exp(-np.asarray(s2)[:,None]*np.array(b))).sum(axis=1)+c


class Model:
    def __init__(self,h,cell,types):
        self.h=np.asarray(h,float)
        self.cell=np.array(cell)
        self.s2=np.sum((self.h/self.cell)**2,axis=1)/4
        self.types=types
        self.n=len(types)
        self.sf=np.array([ff(t,self.s2) for t in types]).T
        self.hr=self.h[:,None,:]*SIGNS[None,:,:]
        self.ht=self.h@TRANS.T

    def calculate(self,z,jac=True):
        n=self.n
        xyz=z[:3*n].reshape(n,3)
        u=z[3*n:4*n]
        scale=np.exp(z[-1])
        phase=np.exp(2j*np.pi*(np.einsum("hsd,ad->has",self.hr,xyz)+self.ht[:,None,:]))
        dw=np.exp(-8*np.pi**2*self.s2[:,None]*u)
        coeff=self.sf*dw
        atom=coeff*phase.sum(axis=2)
        f=atom.sum(axis=1)
        amp=np.sqrt(scale)*abs(f)
        if not jac:
            return f,amp
        dxyz=2j*np.pi*np.einsum("has,hsd->had",phase,self.hr)*coeff[:,:,None]
        du=-8*np.pi**2*self.s2[:,None]*atom
        deriv=np.c_[dxyz.reshape(len(f),-1),du]
        da=np.sqrt(scale)*(f.conj()[:,None]*deriv).real/np.maximum(abs(f)[:,None],1e-12)
        return f,amp,np.c_[da,amp/2]


def fit_model(df,cell,types,xyz,u=None,max_nfev=300,weighted=False):
    model=Model(df[["h","k","l"]].values,cell,types)
    fo=np.sqrt(np.maximum(df.I.values,0))
    if u is None:
        u=np.full(len(types),.02)
    z=np.r_[xyz.ravel(),u,0.]
    _,fc=model.calculate(z,False)
    z[-1]=2*np.log(fo@fc/(fc@fc))
    if weighted:
        sigma=np.sqrt((df.sigma.values**2+(0.03*np.maximum(df.I.values,0))**2)/
                      (4*np.maximum(df.I.values,4))+.2)
    else:
        sigma=np.ones(len(fo))
    bounds=(np.r_[np.full(3*len(types),-5),np.full(len(types),.001),-20],
            np.r_[np.full(3*len(types),5),np.full(len(types),.12),20])
    def fun(z):
        return (model.calculate(z,False)[1]-fo)/sigma
    def jac(z):
        return model.calculate(z,True)[2]/sigma[:,None]
    fit=least_squares(fun,z,jac=jac,bounds=bounds,max_nfev=max_nfev,x_scale="jac",
                      ftol=2e-9,xtol=2e-9,gtol=2e-8)
    _,fc=model.calculate(fit.x,False)
    r1=np.sum(abs(fo-fc))/sum(fo)
    return fit,model,float(r1)


def distances(xyz,cell):
    pairs=[]
    for i,x in enumerate(xyz):
        for j,y in enumerate(xyz):
            for op,(s,t) in enumerate(zip(SIGNS,TRANS)):
                delta=y*s+t-x
                shift=np.rint(delta)
                delta-=shift
                d=np.linalg.norm(delta*cell)
                if .2<d<1.9:
                    pairs.append((i,j,op,float(d),(-shift).astype(int).tolist()))
    return pairs


if __name__=="__main__":
    source=OUT/"phasing_input.csv"
    df=pd.read_csv(source if source.exists() else OUT/"merged.csv")
    cell=json.loads((OUT/("geometry_flexible.json" if source.exists() else "geometry.json")).read_text())["cell"]
    initial=json.loads((OUT/"solution_initial.json").read_text())
    xyz=np.array(initial["xyz"])
    fit,model,r=fit_model(df,cell,["C"]*6,xyz)
    xyz=fit.x[:18].reshape(6,3)%1
    print("all-C",r,"eval",fit.nfev,"xyz",xyz,"u",fit.x[18:24],sep="\n",flush=True)
    print("connectivity",distances(xyz,cell),flush=True)
    best=1
    for n in range(6):
        for oxy in itertools.combinations([i for i in range(6) if i!=n],2):
            types=["C"]*6
            types[n]="N"
            for i in oxy:
                types[i]="O"
            f,m,r=fit_model(df,cell,types,xyz,fit.x[18:24],max_nfev=80)
            if r<best:
                best=r
                result=dict(types=types,xyz=(f.x[:18].reshape(6,3)%1).tolist(),
                            uiso=f.x[18:24].tolist(),log_scale=float(f.x[-1]),R1=r,
                            cell=cell,nfev=f.nfev)
                (OUT/"model_isotropic.json").write_text(json.dumps(result,indent=2))
                print("best",types,r,flush=True)
    print("FINAL",json.dumps(result),flush=True)
