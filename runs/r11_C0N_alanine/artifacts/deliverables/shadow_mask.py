"""Estimate detector obscuration from diffuse counts, independently of Fc."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from geometry import detector_vectors

OUT=Path(__file__).resolve().parent


def create():
    manifest=json.loads((OUT/"frame_manifest.json").read_text())
    yy,xx=np.mgrid[:775,:800]
    results={}
    for run in range(1,7):
        frames=[s for s in manifest if s["run"]==run]
        theta=frames[0]["angles"][0][1]
        s=detector_vectors(xx.ravel(),yy.ravel(),theta,-1,-1,-1,-1)
        angle=np.arccos(np.clip(s[:,0],-1,1)).reshape(xx.shape)
        bins=np.minimum((angle/.025).astype(int),69)
        nbin=70
        maps=[]
        times=[]
        for subset in np.array_split(frames,4):
            total=np.zeros(xx.shape)
            for fr in subset:
                a=np.load(OUT/f"cache/{run}_{fr['frame']}.npy",mmap_mode="r")
                total+=np.minimum(a,2)
            mean=gaussian_filter(total/len(subset),7)
            valid=(xx>15)&(xx<784)&(yy>15)&(yy<759)&((xx<370)|(xx>429))
            med=np.array([np.median(mean[(bins==k)&valid]) if np.sum((bins==k)&valid)>50
                          else np.nan for k in range(nbin)])
            k=np.arange(nbin)
            med=np.interp(k,k[np.isfinite(med)],med[np.isfinite(med)])
            med=gaussian_filter1d(med,1)
            ratio=mean/np.maximum(med[bins],.005)
            maps.append(ratio)
            times.append(np.mean([np.mean(np.array(fr["angles"])[:,0]) for fr in subset]))
        results[f"run{run}"]=np.array(maps,dtype=np.float32)
        results[f"omega{run}"]=np.array(times)
        print("run",run,"masked fraction",np.mean(results[f"run{run}"]<.5),flush=True)
    np.savez_compressed(OUT/"background_transmission.npz",**results)
    return results


def apply(df,results):
    ratios=[]
    for row in df.itertuples():
        run=int(row.run)
        j=np.argmin(abs(results[f"omega{run}"]-row.omega))
        y,x=int(round(row.y)),int(round(row.x))
        # A small footprint protects the aperture from mask-boundary crossings.
        ratio=np.min(results[f"run{run}"][j,y-4:y+5,x-4:x+5])
        ratios.append(ratio)
    df=df.copy()
    df["background_transmission"]=ratios
    df["shadowed"]=df.background_transmission<.55
    df.to_csv(OUT/"unmerged_masked_all.csv",index=False)
    df[~df.shadowed].to_csv(OUT/"unmerged_masked.csv",index=False)
    print("masked",df.shadowed.sum(),"of",len(df),flush=True)
    return df


if __name__=="__main__":
    results=create()
    apply(pd.read_csv(OUT/"unmerged.csv"),results)
