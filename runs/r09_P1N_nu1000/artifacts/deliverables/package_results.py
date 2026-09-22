"""Package this run's source, inputs and checkable outputs; verify member hashes."""
from pathlib import Path
import hashlib,json,zipfile,datetime
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent

def main():
    raw=OUT/'raw_reproduction_20260922_025724'/'deliverables'
    if not raw.is_dir():raw=OUT/'raw_reconstruction_records'
    mainmodel=json.loads((OUT/'final_model.json').read_text());other=json.loads((raw/'final_model.json').read_text())
    cell=np.array([[39.19,-39.19/2,0],[0,39.19*np.sqrt(3)/2,0],[0,0,16.61]])
    shifts=[]
    for a,b in zip(mainmodel['atoms'],other['atoms']):
        assert a['label']==b['label'];d=np.array(a['xyz'])-b['xyz'];d-=np.rint(d);shifts.append(float(np.linalg.norm(cell@d)))
    ri=mainmodel['final_statistics']['refinement']['R1_obs'];rr=other['final_statistics']['refinement']['R1_obs']
    valid=json.loads((raw/'final_independent_validation.json').read_text());geo=json.loads((raw/'final_validation.json').read_text())
    comp={'source_directory':str(raw.relative_to(OUT)),'reconstructed_R1_obs':rr,'main_R1_obs':ri,'absolute_R1_difference':abs(rr-ri),
        'max_non_H_coordinate_difference_A':max(shifts),'rms_non_H_coordinate_difference_A':float(np.sqrt(np.mean(np.square(shifts)))),
        'reconstructed_inventory':geo['inventory'],'reconstructed_node_sizes':geo['node_sizes'],'reconstructed_ligand_inventories':geo['ligand_inventories'],
        'reconstructed_periodic_rank':geo['periodic_connectivity_rank'],'independent_validation_passed':valid['all_passed'],
        'passed':max(shifts)<.001 and abs(rr-ri)<1e-5 and valid['all_passed'] and geo['periodic_connectivity_rank']==3}
    assert comp['passed'];(OUT/'raw_reconstruction_comparison.json').write_text(json.dumps(comp,indent=2))
    files=[]
    for ext in ['*.py','*.json','*.md','*.cif','*.xyz','*.fcf','*.csv','*.log']:
        files.extend((p,'deliverables/'+p.name) for p in OUT.glob(ext) if p.name not in ['package_manifest.json','package_integrity.json'])
    needed=['merged_equal.npz','merged.npz','framework_final_unmasked_fc.npz','framework_final_unmasked_covariance.npz',
        'framework_final_unmasked_diff_map.npy','framework_conservative_cv_fc.npz','framework_h_highangle_fc.npz','final_covariance.npz']
    for name in needed:files.append((OUT/name,'deliverables/'+name))
    for name in ['start.hkl','start.ins','synthesis.png']:files.append((ROOT/'inputs'/name,'inputs/'+name))
    # Keep evidence of the fresh end-to-end test without duplicating its large maps.
    for ext in ['*.log','final_*.json']:
        for p in raw.glob(ext):files.append((p,'deliverables/raw_reconstruction_records/'+p.name))
    files.sort(key=lambda item:item[1]);assert len({name for p,name in files})==len(files)
    manifest={'created_at':datetime.datetime.now().isoformat(),'members':[],'note':'Excluded large exploratory density arrays remain in the original deliverables directory. Manifest does not hash itself.'}
    for p,name in files:
        content=p.read_bytes();manifest['members'].append({'path':name,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()})
    mp=OUT/'package_manifest.json';mp.write_text(json.dumps(manifest,indent=2));files.append((mp,'deliverables/package_manifest.json'))
    archive=OUT/'nu1000_framework_bundle.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p,name in files:z.write(p,name)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for item in manifest['members']:assert hashlib.sha256(z.read(item['path'])).hexdigest()==item['sha256'],item['path']
    info={'archive':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'member_count':len(files),
        'CRC_check_passed':True,'all_member_SHA256_checks_passed':True,'raw_reconstruction_comparison':comp,'finished_at':datetime.datetime.now().isoformat()}
    (OUT/'package_integrity.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))

if __name__=='__main__':main()
