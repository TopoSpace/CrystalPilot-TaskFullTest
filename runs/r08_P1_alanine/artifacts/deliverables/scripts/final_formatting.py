"""Prepare publication-compatible LIST 4 while preserving LIST 6 map output."""
import json
import re
import shutil
from run_command import ROOT,run
source=ROOT/'work/refine/07_final'
work=ROOT/'work/refine/10_publication'
work.mkdir(exist_ok=False)
text=(source/'alanine.ins').read_text()
text=re.sub(r'^TITL.*$','TITL alanine final LIST4 and geometry tables',text,count=1,flags=re.M)
text=text.replace('LIST 6','LIST 4\nCONF')
# Report three clear ammonium H-bonds, not incidental long/weak contacts.
text=text.replace('HTAB N1 O2_$2\n','').replace('EQIV $4 -x+1, y-1/2, -z+1/2\nHTAB C2 O2_$4\n','')
(work/'alanine.ins').write_text(text,encoding='ascii')
shutil.copy2(source/'alanine.hkl',work/'alanine.hkl')
code=run(str(work.relative_to(ROOT)),'10_publication',[r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe','alanine'])
assert code==0
assert 'REM wR2 = 0.0732' in (work/'alanine.res').read_text()
backup=ROOT/'work/validation/pre_LIST4_delivery'
backup.mkdir(exist_ok=False)
for ext in ['ins','res','hkl','fcf','lst','cif']:
    shutil.copy2(ROOT/f'deliverables/alanine.{ext}',backup/f'alanine.{ext}')
shutil.copy2(source/'alanine.fcf',ROOT/'deliverables/provenance/alanine_fourier_LIST6.fcf')
(work/'stage_provenance.json').write_text(json.dumps({'source':'work/refine/07_final/alanine.ins','changes':['LIST 4 instead of LIST 6','CONF torsion table','only three strong ammonium H-bonds reported'],'unchanged':'All reflection data, fitted parameters, restraints, weights and resolution'},indent=2),encoding='utf-8')
