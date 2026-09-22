"""Reproduce using this interpreter only. No external crystallographic executables.
Default: independently validate the final text outputs without regenerating them.
--mode refit: repeat the exact final optimization from the saved CV checkpoint.
--mode from-inputs: regenerate a fresh raw-input solution/refinement in a NEW directory.
"""
from pathlib import Path
import argparse,subprocess,sys,os,json,datetime,hashlib,shutil

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'

def run(script,args=(),folder=OUT,log_prefix='reproduce'):
    cmd=[sys.executable,str(folder/script),*args]
    start=datetime.datetime.now().isoformat();print('RUN',cmd,flush=True)
    result=subprocess.run(cmd,cwd=folder.parent,capture_output=True,text=True,encoding='utf-8',errors='replace',env={**os.environ,'PYTHONIOENCODING':'utf-8'})
    stem=Path(script).stem;log=folder/(log_prefix+'_'+stem+'.log')
    log.write_text('Started: '+start+'\nCommand: '+json.dumps(cmd)+'\n'+result.stdout+'\nSTDERR:\n'+result.stderr+'\nExit: '+str(result.returncode)+'\n',encoding='utf-8')
    print(result.stdout,flush=True)
    if result.stderr:print(result.stderr,file=sys.stderr,flush=True)
    result.check_returncode()

def verify_inputs():
    provenance=json.loads((OUT/'provenance.json').read_text())
    for name,item in provenance['input_hashes'].items():
        assert hashlib.sha256((ROOT/'inputs'/name).read_bytes()).hexdigest()==item['sha256'],'Input hash mismatch: '+name

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['verify','finalize','refit','from-inputs'],default='verify');args=parser.parse_args()
    verify_inputs()
    if args.mode=='verify':run('verify_final.py');return
    if args.mode=='finalize':
        run('finalize.py');run('verify_final.py');return
    if args.mode=='refit':
        run('restrained_refine.py',['--input','framework_conservative_cv','--tag','reproduced_final','--dispersion','--hydrogens','--geometry-strength','4','--all-data','--epochs','5'])
        import numpy as np
        original=json.loads((OUT/'framework_final_unmasked.json').read_text());new=json.loads((OUT/'reproduced_final.json').read_text())
        delta=float(np.max(np.abs(np.array([a['xyz'] for a in original['atoms']])-np.array([a['xyz'] for a in new['atoms']]))))
        io=np.load(OUT/'framework_final_unmasked_fc.npz')['ic'];ir=np.load(OUT/'reproduced_final_fc.npz')['ic'];error=float(np.max(abs(io-ir))/max(io))
        result={'max_fractional_coordinate_difference':delta,'intensity_difference_normalized_to_max':error,'passed':delta<1e-7 and error<1e-6}
        (OUT/'reproduction_comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));assert result['passed'];return
    # This creates a new directory rather than overwriting the delivered history.
    dest=OUT/('raw_reproduction_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    dest.mkdir(exist_ok=False);(dest/'inputs').mkdir();(dest/'deliverables').mkdir()
    for name in ['start.hkl','start.ins','synthesis.png']:shutil.copy2(ROOT/'inputs'/name,dest/'inputs'/name)
    for source in OUT.glob('*.py'):shutil.copy2(source,dest/'deliverables'/source.name)
    target=dest/'deliverables'
    for script in ['crystal.py','diagnose.py','data_quality.py','solve_zr.py','refine.py','build_framework.py']:run(script,folder=target,log_prefix='raw')
    # A transparent compact path from the recovered Fourier framework. Nonlinear
    # histories differ from exploratory checkpoints: do not promise bitwise identity.
    stages=[('framework_iso1','raw_nodisp',False,False,1,5),('raw_nodisp','raw_disp',True,False,1,5),
            ('raw_disp','raw_h',True,True,4,5),('raw_h','framework_conservative_cv',True,True,4,5)]
    for source,tag,disp,H,strength,epochs in stages:
        opt=['--input',source,'--tag',tag,'--geometry-strength',str(strength),'--epochs',str(epochs)]
        if disp:opt+=['--dispersion']
        if H:opt+=['--hydrogens']
        run('restrained_refine.py',opt,folder=target,log_prefix=tag)
    run('restrained_refine.py',['--input','framework_conservative_cv','--tag','framework_h_highangle','--dispersion','--hydrogens','--geometry-strength','4','--dmax','3','--epochs','4'],folder=target,log_prefix='highangle')
    run('restrained_refine.py',['--input','framework_conservative_cv','--tag','framework_final_unmasked','--dispersion','--hydrogens','--geometry-strength','4','--all-data','--epochs','5'],folder=target,log_prefix='final')
    run('finalize.py',folder=target,log_prefix='raw');run('verify_final.py',folder=target,log_prefix='raw')
    print('Fresh raw-input reconstruction:',dest)

if __name__=='__main__':main()
