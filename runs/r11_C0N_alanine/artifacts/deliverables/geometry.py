import numpy as np
from scipy.spatial.transform import Rotation

UB0=np.array([[-.0036164,.1093447,-.0225878],
              [.0317839,-.0440864,-.0514390],
              [-.1181851,-.0152375,-.0131508]])
WAVELENGTH=.71073


def rot(axis,deg):
    return Rotation.from_rotvec(np.asarray(axis)*np.deg2rad(deg)).as_matrix()


def rz(deg):
    t=np.deg2rad(deg)
    c,s=np.cos(t),np.sin(t)
    r=np.zeros(np.shape(t)+(3,3))
    r[...,0,0]=r[...,1,1]=c
    r[...,0,1]=-s
    r[...,1,0]=s
    r[...,2,2]=1
    return r


def detector_vectors(x,y,theta, sx=1,sy=1,st=1,sm=1,dd=47,cx=406.52,cy=380.45):
    module=np.asarray(x)>=400
    mx=np.where(module,605.968,192.5)
    md=np.where(module,59.857,60)
    ang=np.where(module,19.016,-19)*sm
    local=np.stack([md,(np.asarray(x)-mx)*.1*sx,(np.asarray(y)-cy)*.1*sy],axis=-1)
    v=np.einsum("...ij,...j->...i",rz(ang),local)
    # Mechanical translation from the module calibration distance to the scan distance.
    v[...,0] += dd-60
    v[...,1] += (400-cx)*.1*sx
    v=np.einsum("...ij,...j->...i",rz(st*theta),v)
    return v/np.linalg.norm(v,axis=-1)[...,None]


def gonio(omega,kappa,phi,so=1,sk=1,sp=1,axis=0,alpha=49.92498):
    a=np.deg2rad(alpha)
    ax=np.array([np.sin(a) if axis==0 else 0,
                 np.sin(a) if axis==1 else 0,np.cos(a)])
    k=np.asarray(kappa)
    rk=Rotation.from_rotvec(k[...,None]*np.deg2rad(sk)*ax).as_matrix()
    return rz(so*omega)@rk@rz(sp*phi)


def calibrated_q(p, params):
    dd,cx,cy,tzero,alpha,kzero,ozero = params
    s=detector_vectors(p[:,2],p[:,3],p[:,7]+tzero,-1,-1,-1,-1,dd,cx,cy)
    r=gonio(p[:,6]+ozero,p[:,8]+kzero,p[:,9],-1,-1,-1,0,-alpha)
    return np.einsum("...ji,...j->...i",r,s-[1,0,0])


def project(s,theta,params):
    dd,cx,cy,tzero,alpha,kzero,ozero=params
    v=np.asarray(s)@rz(theta+tzero).T
    result=[]
    for mod in [0,1]:
        md=[60,59.857][mod]
        mx=[192.5,605.968][mod]
        ang=[19,-19.016][mod]
        r=rz(ang)
        n=r[:,0]
        origin=md*n+np.array([dd-60,(cx-400)*.1,0])
        t=(origin@n)/(v@n)
        local=(v*t[...,None]-origin)@r
        x=mx-local[...,1]/.1
        y=cy-local[...,2]/.1
        result.append(np.stack([x,y,t],axis=-1))
    return result
