"""Generate a final task-local integrity manifest; do not modify science outputs."""
import json,hashlib,datetime,sys
from raw_frames import OUT

PRIMARY=['README.md','REPORT_zh.md','METHODS_zh.md','alanine.cif','alanine.fcf','observed.hkl','refined_model.json','refinement_reflections.csv','validation.json','algorithm_tests.json','reproducibility_comparison.json','run_pipeline.py']

def main():
 for f in PRIMARY:
  if not (OUT/f).is_file():raise FileNotFoundError(f)
 checks=json.loads((OUT/'validation.json').read_text());tests=json.loads((OUT/'algorithm_tests.json').read_text());rep=json.loads((OUT/'reproducibility_comparison.json').read_text())
 assert checks['all_basic_checks_pass'] and checks['quality_targets_met'] and tests['all_pass'] and rep['passed']
 source=json.loads((OUT/'reproduction_verified'/'run_provenance.json').read_text())
 # Source files that implement the executed scientific pipeline must not have
 # changed since the verified run; later reporting-only scripts are separate.
 mismatches=[]
 dependencies=sorted({step['script'] for step in source['steps']}|{'solve.py','geometry_probe.py','run_pipeline.py'})
 for name in dependencies:
  sha=hashlib.sha256((OUT/name).read_bytes()).hexdigest()
  if sha!=source['source_sha256'][name]:mismatches.append(name)
 assert not mismatches,('Scientific source changed after verified execution',mismatches)
 entries=[]
 for p in sorted(OUT.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts or p.name in ['MANIFEST.sha256','delivery_manifest.json']:continue
  b=p.read_bytes();entries.append(dict(path=p.relative_to(OUT).as_posix(),bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
 manifest=dict(created=datetime.datetime.now().astimezone().isoformat(),primary_files=PRIMARY,files=entries,total_bytes=sum(e['bytes'] for e in entries),file_count=len(entries),excluded=['MANIFEST.sha256 itself','delivery_manifest.json itself','__pycache__'],final_model='alanine.cif',main_statistics=json.loads((OUT/'refined_model.json').read_text())['statistics'],basic_validation_passed=True,quality_targets_met=True,full_raw_data_reproduction_passed=True,verified_pipeline_sources_unchanged=True,independent_standard_crystallographic_validation=False,absolute_configuration='unknown',inputs_modified=False,raw_frame_hash_check=tests['tests']['all_616_raw_frame_hashes_unchanged'])
 (OUT/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
 (OUT/'MANIFEST.sha256').write_text('\n'.join(e['sha256']+'  '+e['path'] for e in entries)+'\n',encoding='utf8')
 print('Final delivery:',len(entries),'files;',manifest['total_bytes'],'bytes; scientific-source hashes match verified run.',flush=True)
 print('Primary structure:',OUT/'alanine.cif',flush=True)
if __name__=='__main__':main()
