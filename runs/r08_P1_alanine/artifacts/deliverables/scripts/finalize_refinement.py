"""Refine recommended weights, explicit H-bonds, and resolution controls."""
import re
import shutil
import json
from run_command import ROOT, run
EXE=r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe'

def stage(source,name,cutoff=None,weight=True,hbonds=False):
    raw=source.read_text()
    text=raw.split('HKLF')[0]
    text=re.sub(r'^TITL.*$',f'TITL alanine {name}',text,count=1,flags=re.M)
    if weight:
        recommended=re.findall(r'^WGHT\s+(.+)$',raw,re.M)[-1]
        text=re.sub(r'^WGHT.*$','WGHT '+recommended,text,count=1,flags=re.M)
    if cutoff is not None:
        text=re.sub(r'^SHEL.*$',f'SHEL 999 {cutoff}',text,flags=re.M)
    if hbonds:
        instructions='\n'.join(line for line in raw.split('END',1)[1].splitlines() if line.startswith(('EQIV ','HTAB ')))
        assert instructions
        text=re.sub(r'^HTAB\s*$',instructions,text,count=1,flags=re.M)
    text+='HKLF 4\nEND\n'
    work=ROOT/'work/refine'/name
    work.mkdir(exist_ok=False)
    (work/'alanine.ins').write_text(text,encoding='ascii')
    shutil.copy2(source.with_suffix('.hkl'),work/'alanine.hkl')
    (work/'stage_provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'cutoff_override':cutoff,'recommended_weight':weight,'explicit_hbonds':hbonds},indent=2),encoding='utf-8')
    code=run(str(work.relative_to(ROOT)),name,[EXE,'alanine'])
    if code:
        raise SystemExit(code)
    assert 'REM wR2' in (work/'alanine.res').read_text()
    return work/'alanine.res'

source=stage(ROOT/'work/refine/05_NH_free/alanine.res','06_weighted',hbonds=True)
final=stage(source,'07_final')
stage(final,'08_cutoff_080',cutoff=0.80,weight=False)
stage(final,'09_full_measured',cutoff=0.70,weight=False)
