"""Apply conservative instrument-shadow masks defined from raw-pixel background.
No Fo-Fc threshold and no reflection-index exclusion are used.
See raw-pixel diagnostics and evaluate_shadow_background.py for the evidence.
"""
import json
import pickle
import sys
import numpy as np
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.format.image import ImageBool
from dials.array_family import flex
from run_command import ROOT, run

work=ROOT/'work/masked_data'
work.mkdir(parents=True,exist_ok=True)
experiments=ExperimentListFactory.from_json_file(str(ROOT/'work/dials/refined_joint_esd.expt'),check_format=False)
groups=[[0,1],[2],[3],[4,5]]
# Guard rectangles encompass the beamstop/support shadow and its partially
# attenuating boundary. Pixel intervals are half-open in each panel's coordinates.
rectangles={0:(120,360,340,445),1:(250,385,340,435),2:(120,360,340,445),3:(250,385,340,435)}
summary=[]
for g,eids in enumerate(groups):
    archive=np.load(ROOT/f'work/final_data/shadow_diagnostic_group_{g}.npz')
    invalid=[archive['candidate_mask0'].copy(),archive['candidate_mask1'].copy()]
    x0,x1,y0,y1=rectangles[g]
    invalid[0][y0:y1,x0:x1]=True
    masks=tuple(flex.bool(np.ascontiguousarray(~m)) for m in invalid)
    path=work/f'geometry_{g}.mask'
    with path.open('wb') as f:
        pickle.dump(masks,f,pickle.HIGHEST_PROTOCOL)
    for eid in eids:
        experiments[eid].imageset.external_lookup.mask.filename=str(path)
        experiments[eid].imageset.external_lookup.mask.data=ImageBool(masks)
    summary.append({'geometry':g,'experiments':eids,'panel0_guard_rectangle_x0_x1_y0_y1':rectangles[g],'masked_pixels_each_panel':[int(m.sum()) for m in invalid],'mask_fraction_each_panel':[float(m.mean()) for m in invalid],'source':'Raw accumulated counts and spatial background depression; includes persistent-zero mask plus conservative guard rectangle.'})
(work/'mask_definition.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
experiments.as_file(str(work/'masked_input.expt'))
print(json.dumps(summary,indent=2))
commands=[
('22_masked_integrate',['dials.integrate','masked_input.expt','../dials/refined_joint_esd.refl','nproc=4','prediction.d_min=0.65']),
('23_masked_symmetry',['dials.symmetry','integrated.expt','integrated.refl','output.html=None']),
('24_masked_scale',['dials.scale','symmetrized.expt','symmetrized.refl','../dials/scale.phil']),
('25_masked_export',['dials.export','scaled.expt','scaled.refl','format=shelx','intensity=scale','shelx.composition=C3H7NO2','shelx.hklout=alanine.hkl','shelx.ins=alanine.ins']),
('26_masked_audit',[sys.executable,'../../audit_intensities.py','.']),
]
for label,cmd in commands:
    result=run('work/masked_data',label,cmd)
    if result:
        raise SystemExit(result)
