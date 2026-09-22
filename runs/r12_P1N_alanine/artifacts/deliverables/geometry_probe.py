from pathlib import Path
import json,re,struct,itertools
import numpy as np
from scipy.spatial.transform import Rotation
from raw_frames import OUT,INP

def rz(deg):
    t=np.deg2rad(deg);c=np.cos(t);s=np.sin(t)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def rotate_z(v,deg):
    t=np.deg2rad(deg);c=np.cos(t);s=np.sin(t)
    return np.column_stack((c*v[:,0]-s*v[:,1],s*v[:,0]+c*v[:,1],v[:,2]))

def rays(x,y,theta,sgnx=1,sgny=1,factor=1,offset=0,kind=0):
    mod=(x>=400).astype(int)
    origin=np.where(mod,605.968,192.5)
    ang=np.deg2rad(np.where(mod,19.016,-19))
    dx=(x-origin)*.1;dy=(y-387.5)*.1
    if kind==0:
        dd=np.where(mod,59.857,60)
        v=np.column_stack((dd*np.cos(ang)-dx*np.sin(ang)-13,sgnx*(dd*np.sin(ang)+dx*np.cos(ang)),sgny*dy))
    elif kind==1:
        v=np.column_stack((np.full(len(x),47),sgnx*(x-406.52)*.1,sgny*(y-380.45)*.1))
    elif kind==2:
        dd=np.where(mod,47.324,47.901)
        dx=(x-np.where(mod,569.88,241.275))*.1
        v=np.column_stack((dd*np.cos(ang)-dx*np.sin(ang),sgnx*(dd*np.sin(ang)+dx*np.cos(ang)),sgny*(y-380.45)*.1))
    v=rotate_z(v,theta*factor+offset)
    return v/np.linalg.norm(v,axis=1)[:,None]

def main():
    p=np.load(OUT/'peaks.npy');meta=json.loads((OUT/'frame_metadata.json').read_text())
    p=p[p[:,5]>100];p=p[::max(1,len(p)//900)]
    inds=p[:,2].astype(int)
    om=np.array([(meta[i]['omega0']+meta[i]['omega1'])/2 for i in inds])
    th=np.array([meta[i]['theta'] for i in inds]);ka=np.array([meta[i]['kappa'] for i in inds]);ph=np.array([meta[i]['phi'] for i in inds])
    ub=np.array([[-.0036164,.1093447,-.0225878],[.0317839,-.0440864,-.051439],[-.1181851,-.0152375,-.0131508]])
    uinv=np.linalg.inv(ub)
    best=[]
    # rotation options tested against observed integer indices, without atom coordinates.
    transforms=[]
    for os,ps,ks,plane in itertools.product([1,-1],[1,-1],[1,-1],[1,-1,2,-2]):
        axis=np.array([np.sign(plane)*np.sin(np.deg2rad(49.925)),0,np.cos(np.deg2rad(49.925))]) if abs(plane)==1 else np.array([0,np.sign(plane)*np.sin(np.deg2rad(49.925)),np.cos(np.deg2rad(49.925))])
        kval=ka*ks
        kr=Rotation.from_rotvec(np.deg2rad(kval)[:,None]*axis).as_matrix()
        # inverse Rz(om) Rk(ka) Rz(phi)
        ro=np.array([rz(-o*os) for o in om]);rp=np.array([rz(-pp*ps) for pp in ph])
        m=np.einsum('ab,nbc,ndc,nde->nae',uinv,rp,kr,ro)
        transforms.append(((os,ps,ks,plane),m))
    for sx,sy,tf,off,kind in itertools.product([1,-1],[1,-1],[1,-1,2,-2],[0,17,-17],[0,1,2]):
        q=rays(p[:,3],p[:,4],th,sx,sy,tf,off,kind)-[1,0,0]
        for tr,m in transforms:
            h=np.einsum('nij,nj->ni',m,q);err=np.linalg.norm(h-np.rint(h),axis=1)
            score=np.mean(err<.15)
            best.append((score,np.median(err),(sx,sy,tf,off,kind),tr))
    for x in sorted(best,reverse=True)[:20]:print(x,flush=True)
    (OUT/'geometry_trials.json').write_text(json.dumps(sorted(best,reverse=True)[:20],indent=2))
if __name__=='__main__':main()
