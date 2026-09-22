"""Create a fresh, inspectable SHELXL stage from a previous RES (never omit reflections)."""
import argparse
import json
from pathlib import Path
import re
import shutil
from run_command import ROOT, run

p=argparse.ArgumentParser()
p.add_argument('source')
p.add_argument('stage')
p.add_argument('--hkl')
p.add_argument('--add-h',action='store_true')
p.add_argument('--recommended-weight',action='store_true')
p.add_argument('--exti',action='store_true')
p.add_argument('--cutoff',type=float)
p.add_argument('--nitrogen-as-oxygen',action='store_true')
args=p.parse_args()
source=(ROOT/args.source).resolve()
raw=source.read_text()
text=raw.split('HKLF')[0]
text=re.sub(r'^TITL.*$', 'TITL alanine '+args.stage, text, count=1, flags=re.M)
text=re.sub(r'^L\.S\..*$', 'L.S. 20',text,count=1,flags=re.M)
if args.recommended_weight:
    weights=re.findall(r'^WGHT\s+(.+)$',raw,re.M)
    text=re.sub(r'^WGHT.*$', 'WGHT '+weights[-1],text,count=1,flags=re.M)
if args.cutoff is not None:
    text=re.sub(r'^SHEL.*$',f'SHEL 999 {args.cutoff}',text,flags=re.M)
if args.add_h:
    text=re.sub(r'^FVAR', 'HFIX 137 N1 C3\nHFIX 13 C2\nFVAR', text, count=1, flags=re.M)
if args.exti and not re.search(r'^EXTI',text,re.M):
    text+='EXTI 0.01\n'
if args.nitrogen_as_oxygen:
    text=re.sub(r'^(N1\s+)3(\s+)',r'\g<1>4\2',text,flags=re.M)
text+='HKLF 4\nEND\n'
work=ROOT/'work/refine'/args.stage
work.mkdir(parents=True,exist_ok=False)
(work/'alanine.ins').write_text(text,encoding='ascii')
hkl=(ROOT/args.hkl).resolve() if args.hkl else source.with_suffix('.hkl')
shutil.copy2(hkl,work/'alanine.hkl')
(work/'stage_provenance.json').write_text(json.dumps(vars(args),indent=2),encoding='utf-8')
code=run(str(work.relative_to(ROOT)),args.stage,[r'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe','alanine'])
if code or not (work/'alanine.res').exists():
    raise SystemExit(code or 1)
