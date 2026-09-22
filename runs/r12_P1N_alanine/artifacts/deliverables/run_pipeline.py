"""Reproduce this task from read-only raw inputs using only the allowed Python.
Usage (from task root): <environment.json Python> deliverables/run_pipeline.py
The default fresh destination is deliverables/reproduction. Existing outputs
are never overwritten by this driver. No installation, download or network.
"""
import os,sys,json,time,subprocess,platform,argparse,hashlib
from pathlib import Path
import numpy,scipy
ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/'deliverables'

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(SCRIPTS/'reproduction'));args=ap.parse_args()
 envspec=json.loads((ROOT/'environment.json').read_text(encoding='utf8'))
 if Path(sys.executable).resolve()!=Path(envspec['environment']['python']).resolve():raise RuntimeError('Use only the Python interpreter named in environment.json')
 dest=Path(args.output).resolve()
 if not dest.is_relative_to(SCRIPTS) or dest==SCRIPTS:raise ValueError('Choose a fresh subdirectory within deliverables')
 if dest.exists():raise FileExistsError('Preserving existing directory: '+str(dest))
 dest.mkdir(parents=True);logs=dest/'logs';logs.mkdir()
 env=os.environ.copy();env['ALANINE_OUTPUT']=str(dest);env['OPENBLAS_NUM_THREADS']='1'
 env.pop('ALANINE_RADIUS',None);env.pop('ALANINE_ROCK_MIN',None)
 provenance=dict(python=sys.executable,python_version=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),started=time.strftime('%Y-%m-%dT%H:%M:%S%z'),network_used=False,crystallographic_libraries_used=False,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SCRIPTS.glob('*.py')})
 steps=[('raw_frames.py',['--extract']),('geometry.py',[]),('audit_pixels.py',[]),('integrate.py',[]),('merge.py',[]),('search_molecule.py',[]),('refine.py',[]),('validate_geometry.py',[]),('test_algorithms.py',[]),('export_structure.py',[]),('validate_model.py',[])]
 history=[];start=time.time()
 for script,arguments in steps:
  t=time.time();print('START',script,flush=True)
  with (logs/(script+'.log')).open('w',encoding='utf8') as log:
   res=subprocess.run([sys.executable,'-u',str(SCRIPTS/script),*arguments],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
  history.append(dict(script=script,arguments=arguments,exit_code=res.returncode,seconds=time.time()-t))
  provenance.update(steps=history,elapsed_seconds=time.time()-start)
  (dest/'run_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf8')
  print('END',script,'return',res.returncode,'seconds',round(time.time()-t,2),flush=True)
  if res.returncode:raise RuntimeError('Failed '+script+'; see '+str(logs/(script+'.log')))
 print('REPRODUCTION COMPLETE',dest,flush=True)
 print(json.dumps(json.loads((dest/'refined_model.json').read_text())['statistics']),flush=True)
if __name__=='__main__':main()
