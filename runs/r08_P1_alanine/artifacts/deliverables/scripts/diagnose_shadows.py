"""Trace suspicious intensity measurements to raw detector pixels, without model filtering."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import uniform_filter, binary_dilation
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from run_command import ROOT

work = ROOT / 'work/final_data'
exps = ExperimentListFactory.from_json_file(str(work / 'integrated.expt'), check_format=True)
refl = flex.reflection_table.from_file(str(work / 'integrated.refl'))
targets = {(2,1,4),(0,2,3),(0,2,2),(0,2,1),(0,2,6),(0,3,3),(1,0,2),(0,2,4),(0,3,1)}
rows=[]
for i in range(len(refl)):
    hkl=tuple(refl['miller_index'][i])
    if tuple(abs(v) for v in hkl) not in targets:
        continue
    eid=int(refl['id'][i]); panel=int(refl['panel'][i]); x,y,z=refl['xyzcal.px'][i]
    exp=exps[eid]; frame=min(max(int(z)-exp.scan.get_array_range()[0],0),len(exp.imageset)-1)
    image=exp.imageset.get_raw_data(frame)[panel].as_numpy_array()
    ix,iy=int(x),int(y)
    patch=image[max(0,iy-10):min(image.shape[0],iy+11),max(0,ix-10):min(image.shape[1],ix+11)]
    row={'row':i,'exp':eid,'hkl':str(hkl),'panel':panel,'x':round(x,2),'y':round(y,2),'frame':frame+1,'I_sum':round(refl['intensity.sum.value'][i],2),'sigma_sum':round(max(refl['intensity.sum.variance'][i],0)**.5,2),'I_prf':round(refl['intensity.prf.value'][i],2),'sigma_prf':round(max(refl['intensity.prf.variance'][i],0)**.5,2),'patch_sum':float(patch.sum()),'patch_nonzero':int(np.count_nonzero(patch)),'patch_n':patch.size}
    rows.append(row)
with (work/'suspicious_reflections_raw_pixels.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print('Suspicious group measurements:',len(rows))
for r in rows:
    if r['I_prf']<30:
        print(r)
# Aggregate raw pixels by detector geometry, never by agreement with a structure model.
groups={}
for eid,e in enumerate(exps):
    # The imported/refined geometry is shared within sets of scans.
    key=tuple(round(v,5) for p in e.detector for v in p.get_origin())
    groups.setdefault(key,[]).append(eid)
summary=[]
for gi,eids in enumerate(groups.values()):
    sums=[np.zeros((p.get_image_size()[1],p.get_image_size()[0]),dtype=np.float64) for p in exps[eids[0]].detector]
    nframes=0
    for eid in eids:
        for fi in range(len(exps[eid].imageset)):
            data=exps[eid].imageset.get_raw_data(fi)
            for pi,im in enumerate(data):
                sums[pi]+=np.maximum(im.as_numpy_array(),0)
            nframes+=1
    masks=[]
    for pi,s in enumerate(sums):
        # A 5x5 patch with exactly zero counts over an entire geometry group.
        zero_core=uniform_filter(s,size=5,mode='constant',cval=1e6)<1e-8
        # Record only at this stage; no data are removed by this diagnostic.
        candidate=binary_dilation(zero_core,iterations=3)
        masks.append(candidate)
        summary.append({'group':gi,'experiments':eids,'panel':pi,'frames':nframes,'total_pixels':s.size,'zero_core_pixels':int(zero_core.sum()),'candidate_masked_pixels':int(candidate.sum()),'mean_counts_per_pixel_per_frame':float(s.mean()/nframes)})
    np.savez_compressed(work/f'shadow_diagnostic_group_{gi}.npz',sum0=sums[0],sum1=sums[1],candidate_mask0=masks[0],candidate_mask1=masks[1])
(work/'shadow_diagnostic_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
