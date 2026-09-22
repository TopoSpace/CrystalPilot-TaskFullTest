import numpy as np,json
from raw_frames import OUT
from scipy.spatial import cKDTree
x=np.loadtxt(OUT/'scaled_observations.csv',delimiter=',',skiprows=1)
p=np.load(OUT/'indexed_peaks.npy')
print('rows',len(x),'peaks',len(p))
for j in np.argsort(abs(x[:,-1]))[::-1][:15]:
 r=x[j];print('\nOUTLIER',r[[0,1,2,3,4,5,6,7,8,11,12,17,21,22,23]])
 ids=np.where((np.abs(x[:,1:4])==np.abs(r[1:4])).all(axis=1))[0]
 print('equivs run,hkl,om,x,y,I,sig,Iscaled,z',x[ids][:,[0,1,2,3,4,5,6,11,12,21,23]])
 near=p[(p[:,0]==r[0])&(abs(p[:,9]-r[4])<2)]
 if len(near):
  dis=np.linalg.norm(near[:,3:5]-r[5:7],axis=1);ii=np.argsort(dis)[:8];print('nearest peak run/frame/x/y/count/om/hf',np.c_[dis[ii],near[ii][:,[0,1,3,4,5,9,13,14,15,16]]])
print('scale outliers by run',[(r,int(((x[:,0]==r)&(abs(x[:,-1])>5)).sum())) for r in range(1,7)])
