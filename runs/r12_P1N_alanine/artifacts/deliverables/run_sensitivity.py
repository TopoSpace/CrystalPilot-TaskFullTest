"""Independent reductions varying aperture and rocking-window size."""
import os,sys,json,shutil,subprocess,time
import numpy as np
from raw_frames import OUT,ROOT

def main():
 scriptdir=ROOT/'deliverables';results=[]
 for name,radius,window in [('wider_rock',4.5,1.75),('larger_aperture',6.0,1.25)]:
  dest=OUT/'sensitivity'/name
  if dest.exists():raise FileExistsError(dest)
  dest.mkdir(parents=True)
  for fname in ['geometry.json','frame_metadata.json','molecular_solution.json','static_masks.npy']:
   shutil.copy2(OUT/fname,dest/fname)
  if radius>4.5:
   from scipy.ndimage import binary_dilation
   masks=np.load(dest/'static_masks.npy');yy,xx=np.mgrid[-2:3,-2:3]
   masks=np.array([binary_dilation(m,structure=xx*xx+yy*yy<=4) for m in masks])
   np.save(dest/'static_masks.npy',masks)
  env=os.environ.copy();env.update(ALANINE_OUTPUT=str(dest),ALANINE_RADIUS=str(radius),ALANINE_ROCK_MIN=str(window),OPENBLAS_NUM_THREADS='1')
  t=time.time()
  for script in ['integrate.py','merge.py','refine.py']:
   with (dest/(script+'.log')).open('w',encoding='utf8') as log:
    r=subprocess.run([sys.executable,'-u',str(scriptdir/script)],env=env,stdout=log,stderr=subprocess.STDOUT)
   if r.returncode:raise RuntimeError(f'{script} failed: {dest}')
  model=json.loads((dest/'refined_model.json').read_text());base=json.loads((OUT/'refined_model.json').read_text())
  dx=(np.array(model['coordinates'])[:6]-np.array(base['coordinates'])[:6])*np.array(model['cell'])
  # Shared initialization; no symmetry/origin matching is required here.
  row=dict(name=name,radius=radius,minimum_rock_halfwidth=window,statistics=model['statistics'],nonH_displacement_A=np.linalg.norm(dx,axis=1).tolist(),maximum_nonH_displacement_A=float(np.max(np.linalg.norm(dx,axis=1))),bonds=json.loads((dest/'bond_geometry.json').read_text()),seconds=time.time()-t)
  results.append(row);print(json.dumps(row),flush=True)
 (OUT/'sensitivity_summary.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
