"""Check ammonium H density and refine their coordinates/Uiso without restraints."""
import re
import shutil
import json
from run_command import ROOT, run

source=ROOT/'work/refine/03_masked_riding/alanine.res'
raw=source.read_text()
base=raw.split('HKLF')[0]
match=re.search(r'AFIX\s+137\n(H1A.*?H1C[^\n]*\n)AFIX\s+0\n',base,re.S)
assert match
h_lines=match.group(1)
for stage,omit in [('04_NH_omit',True),('05_NH_free',False)]:
    work=ROOT/'work/refine'/stage
    work.mkdir(exist_ok=False)
    text=re.sub(r'^TITL.*$',f'TITL alanine {stage}',base,count=1,flags=re.M)
    if omit:
        replacement=''
    else:
        replacement=h_lines.replace('-1.50000','0.02100')
    text=text.replace(match.group(0),replacement)
    if not omit:
        text=text.replace('BOND $H','BOND $H\nHTAB')
    text+='HKLF 4\nEND\n'
    (work/'alanine.ins').write_text(text,encoding='ascii')
    shutil.copy2(source.with_suffix('.hkl'),work/'alanine.hkl')
    (work/'stage_provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'purpose':'Ammonium H omit difference density' if omit else 'Unrestrained isotropic refinement of three ammonium H atoms; C-H remains riding','weight':'Unchanged WGHT 0.05','no_reflection_omission':True},indent=2),encoding='utf-8')
    code=run(str(work.relative_to(ROOT)),stage,[r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe','alanine'])
    if code:
        raise SystemExit(code)
