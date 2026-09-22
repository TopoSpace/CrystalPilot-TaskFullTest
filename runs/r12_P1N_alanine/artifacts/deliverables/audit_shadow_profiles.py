"""Numeric audits of measured rocking profiles and raw backgrounds."""
import json
import numpy as np
from scipy.ndimage import uniform_filter
from raw_frames import OUT,INP,read_frame
from integrate import aperture,rock_width

def main():
 out=[]
 for run,key in [(6,[1,0,2]),(3,[0,3,1]),(4,[0,5,7]),(3,[2,3,1]),(2,[2,4,14])]:
  data=np.load(OUT/f'integrated_run{run}.npy');prof=np.load(OUT/f'profiles_run{run}.npy')
  occ=np.load(OUT/f'pixel_audit_run{run}.npz')['occupied']
  for i in np.where((abs(data[:,1:4])==key).all(axis=1))[0]:
   a=data[i];pr=prof[prof[:,0]==i];x,y=np.rint(a[5:7]).astype(int)
   if len(pr)==0:continue
   sm5=uniform_filter(occ,5);sm15=uniform_filter(occ,15);sm41=uniform_filter(occ,41)
   row=dict(run=run,hkl=a[1:4].tolist(),omega=float(a[4]),xy=a[5:7].tolist(),Linv=float(a[8]),I=float(a[11]),sigma=float(a[12]),background_mean=float(a[16]/max(a[17],1)),occupancy_5_15_41=[float(sm[y,x]) for sm in [sm5,sm15,sm41]],profile=pr.tolist())
   frame=int(pr[np.argmax(pr[:,3]),1]);im,_=read_frame(INP/'frames'/f'pgw240033_Mo_{run}_{frame}.rodhypix')
   row['max_frame']=frame;row['max_frame_apertures']={str(r):list(map(float,aperture(im,a[5],a[6],r)[:4])) for r in [3,4.5,6,8]}
   print(json.dumps(row),flush=True);out.append(row)
 (OUT/'shadow_profile_audit.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
