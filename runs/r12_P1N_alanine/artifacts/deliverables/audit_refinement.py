import numpy as np,json
from raw_frames import OUT
r=np.loadtxt(OUT/'refinement_reflections.csv',delimiter=',',skiprows=1)
u=np.loadtxt(OUT/'scaled_observations.csv',delimiter=',',skiprows=1)
rej=np.loadtxt(OUT/'rejected_observations.csv',delimiter=',',skiprows=1)
print('Worst F2 residuals: h k l Fo2 sigma Fc2 z, then equivalent observations')
for i in np.argsort(abs(r[:,-1]))[::-1][:25]:
 a=r[i];print('REF',' '.join(f'{v:.4g}' for v in a[[0,1,2,3,4,5,9]]))
 rows=u[(np.abs(u[:,1:4])==a[:3]).all(axis=1)]
 print('OBS [run omega x y I sigma z]')
 for z in rows:print(' '.join(f'{v:.4g}' for v in z[[0,4,5,6,21,22,23]]))
