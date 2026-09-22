"""Audit raw-frame occupancy and invalid pixels without a structure model."""
import json,time
import numpy as np
from scipy.ndimage import uniform_filter,binary_dilation
from concurrent.futures import ProcessPoolExecutor
from raw_frames import OUT,INP,read_frame

def audit_run(run):
 meta=json.loads((OUT/'frame_metadata.json').read_text())
 meta=[m for m in meta if int(m['file'].split('_')[-2])==run]
 occupied=np.zeros((775,800),float);invalid=occupied.copy();clipped=occupied.copy();rows=[]
 for m in meta:
  im,_=read_frame(INP/'frames'/m['file'])
  occupied+=im>0;invalid+=im<0;clipped+=np.clip(im,0,2)
  rows.append([int(m['file'].split('_')[-1].split('.')[0]),int(np.sum(im<0)),int(np.sum(im>0))])
 occupied/=len(meta);invalid/=len(meta);clipped/=len(meta)
 np.savez_compressed(OUT/f'pixel_audit_run{run}.npz',occupied=occupied,invalid=invalid,clipped=clipped)
 points=[(246,381),(351,381),(293,381),(595,381),(549,418),(290,504),(290,259),(352,371),(200,381),(100,381)]
 info=[]
 for x,y in points:
  box=np.s_[y-5:y+6,x-5:x+6]
  info.append(dict(x=x,y=y,occupied=float(occupied[box].mean()),invalid=float(invalid[box].mean()),clipped=float(clipped[box].mean())))
 row=dict(run=run,frame_invalid_counts=rows,points=info)
 print(json.dumps(dict(run=run,invalid_range=[min(r[1] for r in rows),max(r[1] for r in rows)],points=info)),flush=True)
 return row

def main():
 if (OUT/'refined_model.json').exists() and not (OUT/'baseline_before_pixel_mask.json').exists():
  model=json.loads((OUT/'refined_model.json').read_text())
  (OUT/'baseline_before_pixel_mask.json').write_text(json.dumps(model,indent=2))
  np.savez_compressed(OUT/'baseline_before_pixel_mask.npz',merged=np.load(OUT/'merged.npy'),reflections=np.loadtxt(OUT/'refinement_reflections.csv',delimiter=',',skiprows=1))
 with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(audit_run,range(1,7)))
 (OUT/'pixel_audit.json').write_text(json.dumps(rows,indent=2))
 build_masks()

def build_masks():
 masks=[];summary=[]
 yy,xx=np.mgrid[-6:7,-6:7];disk=xx*xx+yy*yy<=6**2
 for run in range(1,7):
  a=np.load(OUT/f'pixel_audit_run{run}.npz')['occupied']
  sm=uniform_filter(a,5)
  broad=uniform_filter(a,81)
  dark=(sm<.005)|(sm<.45*broad)
  mask=binary_dilation(dark,structure=disk)
  masks.append(mask)
  summary.append(dict(run=run,dark_pixels=int(dark.sum()),dilated_pixels=int(mask.sum()),occupancy_quantiles=np.percentile(sm[:,10:375],[0,1,5,10,50,90,100]).tolist()))
 np.save(OUT/'static_masks.npy',np.array(masks))
 (OUT/'mask_settings.json').write_text(json.dumps(dict(method='5x5 mean raw photon occupancy < 0.005 OR <45 percent of local 81x81 mean; dilate by 6-pixel disk to reject partially shadowed apertures',structural_model_used=False,summary=summary),indent=2))
 print('Static raw-background masks',json.dumps(summary),flush=True)
if __name__=='__main__':main()
