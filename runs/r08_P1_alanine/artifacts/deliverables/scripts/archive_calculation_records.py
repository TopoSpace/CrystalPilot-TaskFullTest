"""Package this run only; never include or change raw inputs or other runs."""
import datetime
import hashlib
import json
import shutil
import zipfile
from pathlib import Path
import gemmi
from run_command import ROOT

out=ROOT/'deliverables'
assert out.is_dir()
for name in ['commands.jsonl','raw_frames_manifest.json','environment.json','PROCESSING_NOTES.md']:
    shutil.copy2(ROOT/name,out/'provenance'/name)
checks=ROOT/'work/validation_list4'
for name in ['alanine.chk','alanine.ckf','alanine.vrf','platon.out','console.log']:
    shutil.copy2(checks/'platon'/name,out/'validation'/name)
shutil.copy2(checks/'validation.json',out/'validation/validation.json')
# Keep script copies easy to inspect; zip also puts them at their original root.
scripts=out/'scripts'
scripts.mkdir(exist_ok=True)
for path in ROOT.glob('*.py'):
    shutil.copy2(path,scripts/path.name)
archive=out/'calculation_records.zip'
assert not archive.exists(),archive
sources=sorted(p for p in (ROOT/'work').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
sources+=sorted(ROOT.glob('*.py'))
sources += [ROOT/x for x in ['commands.jsonl','raw_frames_manifest.json','environment.json','PROCESSING_NOTES.md']]
with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for path in sources:
        z.write(path,path.relative_to(ROOT).as_posix())
with zipfile.ZipFile(archive) as z:
    corrupt=z.testzip()
    assert corrupt is None,corrupt
    members=len(z.infolist())
# Check delivered outputs, rather than assuming successful copies imply validity.
b=gemmi.cif.read_file(str(out/'alanine.cif')).sole_block()
f=gemmi.cif.read_file(str(out/'alanine.fcf')).sole_block()
assert int(f.find_value('_shelx_refln_list_code'))==4
assert len(f.find_values('_refln_index_h'))==int(b.find_value('_refine_ls_number_reflns'))==964
assert 'LIST 4' in (out/'alanine.ins').read_text()
assert not any(line.startswith('OMIT ') for line in (out/'alanine.ins').read_text().splitlines())
assert (out/'alanine.hkl').read_bytes()==(ROOT/'work/masked_data/alanine.hkl').read_bytes()
validation=json.loads((out/'validation/validation.json').read_text())
assert validation['all_adps_positive']
assert validation['raw_input_hashes']['checked']==646 and not validation['raw_input_hashes']['changed']
assert validation['clean_room_refinement_rerun']['exit_code']==0
assert validation['platon']['alert_counts'].get('B',0)==0
summary={'created_local':datetime.datetime.now().isoformat(),'archive_members':members,'archive_bytes':archive.stat().st_size,'archive_crc_test':'passed','structure_files_match_final_calculation':True,'final_reflections':964,'platon_alerts':validation['platon']['alert_counts'],'scientific_limitations':'See README_zh.md and validation/ALERT_RESPONSES.md; not an unconditional publication-ready or online checkCIF pass.'}
(out/'validation/package_check.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
manifest=[]
for path in sorted(out.rglob('*')):
    if path.is_file() and path.name!='SHA256SUMS.json':
        manifest.append({'path':path.relative_to(out).as_posix(),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(out/'SHA256SUMS.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
# Verify the just-built manifest once, including the zip.
for entry in manifest:
    assert hashlib.sha256((out/entry['path']).read_bytes()).hexdigest()==entry['sha256']
print(json.dumps({**summary,'deliverable_files_hashed':len(manifest),'hash_manifest_verification':'passed'},indent=2))
