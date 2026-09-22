"""Anisotropic F^2 refinement with geometrically riding hydrogen atoms."""
from pathlib import Path
import itertools
import json
import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from refine import ff,SCATTER
from solve import SIGNS,TRANS

OUT=Path(__file__).resolve().parent
TYPES=["O","O","N","C","C","C"]
LABELS=["O1","O2","N1","C1","C2","C3"]
H_LABELS=["H2","H1N","H2N","H3N","H1C","H2C","H3C"]
H_PARENTS=np.array([4,2,2,2,5,5,5])
H_FACTORS=np.array([1.2,1.5,1.5,1.5,1.5,1.5,1.5])


def unit(v):
    return v/np.linalg.norm(v)


def nearest_image(y,x,cell):
    images=y*SIGNS+TRANS
    images-=np.rint(images-x)
    return images[np.argmin(np.linalg.norm((images-x)*cell,axis=1))]


def canonical_molecule(model):
    x=np.array(model["xyz"])
    cell=np.array(model["cell"])
    types=np.array(model["types"])
    carbons=np.where(types=="C")[0]
    oxygens=np.where(types=="O")[0]
    nitrogen=int(np.where(types=="N")[0][0])
    d=np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            d[i,j]=np.linalg.norm((nearest_image(x[j],x[i],cell)-x[i])*cell)
    carboxyl=int(carbons[np.argmin([sum(d[c,oxygens]) for c in carbons])])
    other=carbons[carbons!=carboxyl]
    alpha=int(other[np.argmin(d[other,nitrogen])])
    methyl=int(other[other!=alpha][0])
    order=[*oxygens,nitrogen,carboxyl,alpha,methyl]
    xyz=x[order].copy()
    for i in [2,3,5]:
        xyz[i]=nearest_image(xyz[i],xyz[4],cell)
    for i in [0,1]:
        xyz[i]=nearest_image(xyz[i],xyz[3],cell)
    assert all(np.linalg.norm((xyz[i]-xyz[4])*cell)<1.8 for i in [2,3,5])
    assert all(np.linalg.norm((xyz[i]-xyz[3])*cell)<1.5 for i in [0,1])
    return xyz,np.array(model["uiso"])[order]


def hydrogens(xyz,cell,torsions):
    xyz=xyz*np.asarray(cell)
    h=[]
    direction=-sum(unit(xyz[i]-xyz[4]) for i in [2,3,5])
    h.append(xyz[4]+unit(direction)*1.00)
    for parent,ref,length,torsion in [(2,3,.91,torsions[0]),(5,2,.98,torsions[1])]:
        axis=unit(xyz[4]-xyz[parent])
        v=xyz[ref]-xyz[4]
        e1=unit(v-axis*np.dot(axis,v))
        e2=np.cross(axis,e1)
        for angle in torsion+np.arange(3)*2*np.pi/3:
            direction=-axis/3+np.sqrt(8/9)*(np.cos(angle)*e1+np.sin(angle)*e2)
            h.append(xyz[parent]+length*direction)
    return np.array(h)/cell


def u_matrices(u):
    a=np.zeros((len(u),3,3))
    a[:,0,0],a[:,1,1],a[:,2,2]=u[:,:3].T
    a[:,1,2]=a[:,2,1]=u[:,3]
    a[:,0,2]=a[:,2,0]=u[:,4]
    a[:,0,1]=a[:,1,0]=u[:,5]
    return a


class AnisoModel:
    def __init__(self,h,cell,hydrogen=True):
        self.h=np.asarray(h,float)
        self.cell=np.array(cell)
        self.s2=np.sum((self.h/self.cell)**2,axis=1)/4
        self.hr=self.h[:,None,:]*SIGNS[None,:,:]
        self.hrmetric=self.hr/self.cell
        self.ht=self.h@TRANS.T
        self.hydrogen=hydrogen
        types=TYPES+(["H"]*7 if hydrogen else [])
        self.sf=np.array([ff(t,self.s2) for t in types]).T

    def calculate(self,z):
        xyz=z[:18].reshape(6,3)
        u=z[18:54].reshape(6,6)
        tensors=u_matrices(u)
        if self.hydrogen:
            hx=hydrogens(xyz,self.cell,z[55:57])
            xyz=np.concatenate([xyz,hx])
            hu=np.mean(u[:,:3],axis=1)[H_PARENTS]*H_FACTORS
            tensors=np.concatenate([tensors,hu[:,None,None]*np.eye(3)[None,:,:]])
        arg=np.einsum("hsd,ade,hse->has",self.hrmetric,tensors,self.hrmetric)
        phase=2j*np.pi*(np.einsum("hsd,ad->has",self.hr,xyz)+self.ht[:,None,:])
        f=(self.sf*np.exp(phase-2*np.pi**2*arg).sum(axis=2)).sum(axis=1)
        intensity=np.exp(z[54])*abs(f)**2
        return f,intensity


def initial_z(model):
    x,uiso=canonical_molecule(model)
    u=np.zeros((6,6))
    u[:,:3]=uiso[:,None]
    return np.r_[x.ravel(),u.ravel(),model["log_scale"],0.,0.]


def optimize(df,cell,z,cycles=3,hydrogen=True):
    mdl=AnisoModel(df[["h","k","l"]].values,cell,hydrogen)
    io=df.I.values
    sig=df.sigma.values
    lo=np.r_[np.full(18,-2),np.tile([.0005,.0005,.0005,-.08,-.08,-.08],6),-10,-30,-30]
    hi=np.r_[np.full(18,3),np.tile([.12,.12,.12,.08,.08,.08],6),10,30,30]
    for cycle in range(cycles):
        fc=mdl.calculate(z)[1]
        p=(np.maximum(io,0)+2*fc)/3
        sigma=np.sqrt(sig**2+(.03*p)**2)
        def res(z):
            return (mdl.calculate(z)[1]-io)/sigma
        fit=least_squares(res,z,bounds=(lo,hi),x_scale="jac",
                          max_nfev=180,ftol=3e-9,xtol=3e-9,gtol=1e-7)
        z=fit.x
        fo=np.sqrt(np.maximum(io,0))
        calc=np.sqrt(mdl.calculate(z)[1])
        r=sum(abs(fo-calc))/sum(fo)
        print("cycle",cycle,"R1",r,"gof",np.sqrt(sum(fit.fun**2)/(len(io)-len(z))),
              "nfev",fit.nfev,flush=True)
    return z,mdl,fit,sigma


def metrics(df,mdl,z,sigma,nparams=57):
    f,ic=mdl.calculate(z)
    io=df.I.values
    fo=np.sqrt(np.maximum(io,0))
    fc=np.sqrt(ic)
    sel=io>2*df.sigma.values
    r=sum(abs(fo-fc))/sum(fo)
    robs=sum(abs(fo[sel]-fc[sel]))/sum(fo[sel])
    chi=sum(((io-ic)/sigma)**2)
    return dict(R1_all=float(r),R1_gt2sigma=float(robs),
                wR2=float(np.sqrt(chi/sum((io/sigma)**2))),
                gof=float(np.sqrt(chi/(len(io)-nparams))),n_reflections=len(io),
                n_observed=int(sel.sum()),n_parameters=nparams,chi2=float(chi))


def difference_map(df,cell,z,prefix=""):
    full=[]
    iobs=[]
    for row in df.itertuples():
        for h in set(tuple(np.array([row.h,row.k,row.l])*s) for s in itertools.product([-1,1],repeat=3)):
            full.append(h)
            iobs.append(row.I)
    h=np.array(full)
    mdl=AnisoModel(h,cell)
    f,ic=mdl.calculate(z)
    fo=np.sqrt(np.maximum(iobs,0)/np.exp(z[54]))
    coeff=(fo-abs(f))*np.exp(1j*np.angle(f))
    shape=np.array([64,64,128])
    grid=np.zeros(shape,complex)
    grid[tuple((h%shape).T)]=coeff
    rho=np.fft.fftn(grid).real/np.prod(cell)
    np.save(OUT/(prefix+"difference_density.npy"),rho.astype(np.float32))
    result={"min_e_A3":float(rho.min()),"max_e_A3":float(rho.max()),
            "min_xyz":(np.array(np.unravel_index(np.argmin(rho),shape))/shape).tolist(),
            "max_xyz":(np.array(np.unravel_index(np.argmax(rho),shape))/shape).tolist(),
            "rms_e_A3":float(np.sqrt(np.mean(rho**2)))}
    print("difference",result,flush=True)
    return result


if __name__=="__main__":
    final="--final" in sys.argv
    prefix="final_" if final else ""
    df=pd.read_csv(OUT/("merged_final.csv" if final else "merged.csv"))
    start=json.loads((OUT/"model_isotropic.json").read_text())
    cell=np.array(json.loads((OUT/"geometry_flexible.json").read_text())["cell"] if final else start["cell"])
    z=initial_z(start)
    best=None
    # Rotational search covers both terminal tetrahedral hydrogen groups.
    for a,b in itertools.product([0,np.pi/3],repeat=2):
        zz=z.copy()
        zz[55:57]=[a,b]
        z1,mdl,fit,sigma=optimize(df,cell,zz,cycles=1)
        score=sum(fit.fun**2)
        if best is None or score<best[0]:
            best=(score,z1)
    z,mdl,fit,sigma=optimize(df,cell,best[1],cycles=3)
    stats=metrics(df,mdl,z,sigma)
    eig=np.linalg.eigvalsh(u_matrices(z[18:54].reshape(6,6)))
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*stats["gof"]**2
    result=dict(cell=cell.tolist(),types=TYPES,labels=LABELS,
                xyz=z[:18].reshape(6,3).tolist(),uij=z[18:54].reshape(6,6).tolist(),
                hydrogen_xyz=hydrogens(z[:18].reshape(6,3),cell,z[55:57]).tolist(),
                hydrogen_labels=H_LABELS,log_scale=float(z[54]),torsions=z[55:57].tolist(),
                parameter_vector=z.tolist(),parameter_esd=np.sqrt(np.diag(cov)).tolist(),
                statistics=stats,u_eigenvalues=eig.tolist(),
                difference_map=difference_map(df,cell,z,prefix))
    (OUT/("model_final.json" if final else "model_anisotropic.json")).write_text(json.dumps(result,indent=2))
    np.save(OUT/(prefix+"parameter_covariance.npy"),cov)
    f,ic=mdl.calculate(z)
    df["Fc2"]=ic
    df["weight"]=1/sigma**2
    df["phase_rad"]=np.angle(f)
    df.to_csv(OUT/(prefix+"refinement_reflections.csv"),index=False)
    print(json.dumps(stats),"\nU eigenvalues",eig,"\nTorsions",z[55:57],flush=True)
