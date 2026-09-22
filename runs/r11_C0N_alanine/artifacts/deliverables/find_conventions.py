from pathlib import Path
import itertools
import numpy as np
from geometry import UB0, detector_vectors, gonio

out=Path(__file__).resolve().parent
p=np.load(out/"peaks.npy")
p=p[p[:,4]>100]
if len(p)>1500:
    p=p[np.linspace(0,len(p)-1,1500).astype(int)]
inv=np.linalg.inv(UB0)
ans=[]
for sx,sy,st,sm,tm,to,bm in itertools.product([1,-1],[1,-1],[1,-1],[1,-1],[1,2],[0,17,-17],[1,-1]):
    s=detector_vectors(p[:,2],p[:,3],p[:,7]*tm+to,sx,sy,st,sm)
    q=(s-np.array([1,0,0]))*bm
    for so,sk,sp,axis,sa in itertools.product([1,-1],[1,-1],[1,-1],[0,1],[1,-1]):
        r=gonio(p[:,6]-.36049,p[:,8]+.29559,p[:,9],so,sk,sp,axis,49.92498*sa)
        q0=np.einsum("...ji,...j->...i",r,q)
        hf=q0@inv.T
        err=np.linalg.norm(hf-np.rint(hf),axis=1)
        score=np.mean(np.exp(-err**2/.01))
        ans.append((score,np.median(err),sx,sy,st,sm,tm,to,bm,so,sk,sp,axis,sa))
for row in sorted(ans,reverse=True)[:30]:
    print(row,flush=True)
