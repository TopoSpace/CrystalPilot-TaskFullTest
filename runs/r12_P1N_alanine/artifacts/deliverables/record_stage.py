"""Preserve immutable checkpoints of this task's own computed outputs."""
import sys,shutil,json,datetime,hashlib
from raw_frames import OUT
NAMES=['refined_model.json','merged.npy','merged_all.csv','refinement_reflections.csv','refinement_covariance.npy','merging_statistics.json','difference_statistics.json','bond_geometry.json','integration_settings.json','mask_settings.json','refinement_stages.json']
def record(name):
 assert name.replace('_','').replace('-','').isalnum()
 dest=OUT/'history'/name
 if dest.exists():raise FileExistsError(dest)
 dest.mkdir(parents=True)
 hashes={}
 for filename in NAMES:
  p=OUT/filename
  if p.exists():shutil.copy2(p,dest/filename);hashes[filename]=hashlib.sha256(p.read_bytes()).hexdigest()
 (dest/'checkpoint.json').write_text(json.dumps(dict(created=datetime.datetime.now().isoformat(),sha256=hashes),indent=2))
 print('Recorded',dest,flush=True)
if __name__=='__main__':record(sys.argv[1])
