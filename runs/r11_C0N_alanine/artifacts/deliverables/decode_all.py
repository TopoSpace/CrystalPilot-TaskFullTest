from pathlib import Path
import concurrent.futures
import hashlib
import json
import sys
import time
import numpy as np
from scipy.ndimage import maximum_filter, uniform_filter
from rawio import read_frame

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
CACHE = OUT / "cache"


def work(path):
    path = Path(path)
    a, meta = read_frame(path)
    run, frame = map(int,path.stem.split("_")[-2:])
    assert a.min() >= 0
    np.save(CACHE/f"{run}_{frame}.npy",a.astype(np.int16))
    angles = np.array([meta["angles_start_steps"],meta["angles_end_steps"]])
    angles = angles / np.array([12800,12800,6400,6400,12800])
    sm = uniform_filter(a.astype(float),3)*9
    sel = (sm==maximum_filter(sm,9)) & (sm>25)
    sel[:8] = False
    sel[-8:] = False
    sel[:,:8] = False
    sel[:,-8:] = False
    sel[:,379:422] = False
    peaks=[]
    for y,x in np.argwhere(sel):
        cut=a[y-4:y+5,x-4:x+5].astype(float)
        bg = (cut.sum()-cut[2:7,2:7].sum())/56
        sig = np.maximum(cut[2:7,2:7]-bg,0)
        s=sig.sum()
        if s < 30:
            continue
        yy,xx=np.mgrid[y-2:y+3,x-2:x+3]
        peaks.append([run,frame,(xx*sig).sum()/s,(yy*sig).sum()/s,s,float(a[y,x]),
                      *angles.mean(axis=0)[:4]])
    return dict(run=run,frame=frame,angles=angles.tolist(),total=int(a.sum()),
                maximum=int(a.max()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                peak_count=len(peaks)),peaks


if __name__ == "__main__":
    CACHE.mkdir(exist_ok=True)
    files=sorted((ROOT/"inputs/alanine/frames").glob("pgw*.rodhypix"),
                 key=lambda p:tuple(map(int,p.stem.split("_")[-2:])))
    t=time.time()
    stats,peaks=[],[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
        for i,(s,p) in enumerate(pool.map(work,files)):
            stats.append(s)
            peaks.extend(p)
            if i%30==0:
                print(i+1,len(files),round(time.time()-t,1),"sec",flush=True)
    (OUT/"frame_manifest.json").write_text(json.dumps(stats,indent=2))
    np.save(OUT/"peaks.npy",np.array(peaks))
    print("DONE",len(stats),len(peaks),time.time()-t)
