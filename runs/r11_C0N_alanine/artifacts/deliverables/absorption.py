"""Uniform-volume ray absorption through the supplied convex face model."""
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import qmc
from scipy.spatial import ConvexHull
from geometry import gonio
from geometry_flexible import unpack

OUT=Path(__file__).resolve().parent


def shape():
    text=(OUT.parent/"inputs/alanine/pgw240033_Mo.CAP_shape").read_text()
    rows=[list(map(float,line.split()[2:])) for line in text.splitlines() if line.startswith("SHAPE FACE")]
    data=np.array(rows)
    distances=data[:,3]
    normals=data[:,4:7]
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    vertices=[]
    for ii in itertools.combinations(range(len(normals)),3):
        mat=normals[list(ii)]
        if abs(np.linalg.det(mat))<1e-7:
            continue
        v=np.linalg.solve(mat,distances[list(ii)])
        if np.all(normals@v<=distances+1e-8):
            vertices.append(v)
    vertices=np.unique(np.round(vertices,10),axis=0)
    lo,hi=vertices.min(axis=0),vertices.max(axis=0)
    grid=qmc.Sobol(3,scramble=True,seed=20260921).random_base2(16)
    points=lo+(hi-lo)*grid
    points=points[np.all(points@normals.T<=distances,axis=1)]
    return normals,distances,vertices,points


def transmission(normals,distances,points,incoming,outgoing):
    mu=.11620
    gaps=distances[None,:]-points@normals.T
    tin=np.min(np.where((-normals@incoming)[None,:]>1e-10,
                        gaps/np.maximum((-normals@incoming)[None,:],1e-20),np.inf),axis=1)
    tout=np.min(np.where((normals@outgoing)[None,:]>1e-10,
                         gaps/np.maximum((normals@outgoing)[None,:],1e-20),np.inf),axis=1)
    return np.mean(np.exp(-mu*(tin+tout)))


if __name__=="__main__":
    normals,distances,vertices,points=shape()
    geom=json.loads((OUT/"geometry_flexible.json").read_text())
    z=np.array(geom["parameter_vector"])
    ub=unpack(z)[0]
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    df=pd.read_csv(OUT/"unmerged_masked_all.csv")
    transmissions=[]
    # Two deterministic quadrature resolutions verify integration convergence.
    errors=[]
    for row in df.itertuples():
        fr=next(fr for fr in manifest if fr["run"]==int(row.run))
        ka,ph=fr["angles"][0][2:4]
        om=row.omega-.36049+np.r_[0,z[20:25]][int(row.run)-1]
        r=gonio(om,ka+z[18],ph,-1,-1,-1,0,-z[19])
        inc=r.T@np.array([1.,0,0])
        out=inc+ub@np.array([row.h,row.k,row.l])
        out/=np.linalg.norm(out)
        t=transmission(normals,distances,points,inc,out)
        transmissions.append(t)
        if len(transmissions)%100==0:
            small=transmission(normals,distances,points[:len(points)//2],inc,out)
            errors.append(abs(t-small))
    df["transmission"]=transmissions
    df["I_before_absorption"]=df.I
    df.I/=df.transmission
    df.sigma/=df.transmission
    df.to_csv(OUT/"unmerged_absorption_all.csv",index=False)
    df[~df.shadowed].to_csv(OUT/"unmerged_absorption.csv",index=False)
    report=dict(mu_mm_inverse=.11620,mu_source="Supplied CAP_shape header; not independently rederived.",
                model="Uniform illuminated convex polyhedron, numerical ray quadrature.",
                vertices_mm=vertices.tolist(),quadrature_points=len(points),
                volume_mm3=float(ConvexHull(vertices).volume),
                xyz_extents_mm=(vertices.max(axis=0)-vertices.min(axis=0)).tolist(),
                Tmin=float(min(transmissions)),Tmax=float(max(transmissions)),
                quadrature_max_half_sample_difference=float(max(errors)))
    (OUT/"absorption.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
