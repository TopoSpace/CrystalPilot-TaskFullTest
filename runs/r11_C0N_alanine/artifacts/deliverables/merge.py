from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

OUT=Path(__file__).resolve().parent


def merge(file="unmerged.csv",output="merged",scan=False,validate_scales=False):
    df=pd.read_csv(OUT/file)
    for name in ["h","k","l"]:
        df[name]=df[name].round().astype(int)
        df["a"+name]=df[name].abs()
    keys=df[["ah","ak","al"]]
    codes,unique=pd.factorize(pd.MultiIndex.from_frame(keys),sort=True)
    I=df.I.values
    sig=df.sigma.values
    runs=df.run.values.astype(int)-1
    scales=np.ones(6)
    design=np.zeros((len(df),17 if scan else 5))
    for run in range(6):
        sel=runs==run
        if run:
            design[sel,run-1]=1
        if scan:
            om=df.omega.values[sel]
            t=2*(om-om.min())/(om.max()-om.min())-1
            design[sel,5+2*run]=t
            design[sel,6+2*run]=t*t-1/3
    beta=np.zeros(design.shape[1])
    rng=np.random.default_rng(6192026)
    validation_codes=rng.random(len(unique))<.2
    training=~validation_codes[codes] if validate_scales else np.ones(len(df),bool)
    robust=np.ones(len(I))
    means=np.array([np.median(I[codes==i]) for i in range(len(unique))])
    for cycle in range(15):
        factors=np.exp(design@beta)
        corr=I/factors
        ss=sig/factors
        z=np.abs(corr-means[codes])/np.sqrt(ss**2+(.04*np.maximum(means[codes],0))**2+1)
        robust=np.minimum(1,4/np.maximum(z,1e-9))**2
        w=robust/(ss**2+(.04*np.maximum(means[codes],0))**2+1)
        means=np.bincount(codes,w*corr)/np.bincount(codes,w)
        use=(means[codes]>10)&(robust>.5)&training
        f=least_squares(lambda z:(I[use]-np.exp(design[use]@z)*means[codes[use]])/
                        np.sqrt(sig[use]**2+(.04*I[use])**2+4),
                        beta,loss="soft_l1",bounds=(-.6,.6),max_nfev=40)
        beta=f.x
        z=np.abs(corr-means[codes])/np.sqrt(ss**2+(.04*np.maximum(means[codes],0))**2+1)
        robust=np.minimum(1,4/np.maximum(z,1e-9))**2
    factors=np.exp(design@beta)
    scales=np.array([np.median(factors[runs==r]) for r in range(6)])
    corr=I/factors
    ss=sig/factors
    merged=[]
    for i,hkl in enumerate(unique):
        ii=np.where(codes==i)[0]
        wi=robust[ii]/(ss[ii]**2+(.04*np.maximum(means[i],0))**2+1)
        mean=np.sum(wi*corr[ii])/wi.sum()
        neff=wi.sum()**2/(wi**2).sum()
        esd=max(np.sqrt(1/wi.sum()),np.sqrt(np.sum(wi*(corr[ii]-mean)**2)/wi.sum()/max(neff,1)))
        merged.append([*hkl,mean,esd,len(ii),df.d.values[ii].mean()])
    m=pd.DataFrame(merged,columns=["h","k","l","I","sigma","multiplicity","d"])
    absences=np.zeros(len(m),bool)
    for axis in range(3):
        hh=m[["h","k","l"]].values
        absent=(np.count_nonzero(hh,axis=1)==1)&(hh[:,axis]%2==1)
        absences|=absent
    m["absent_212121"]=absences
    m.to_csv(OUT/f"{output}_all.csv",index=False)
    m[~absences].to_csv(OUT/f"{output}.csv",index=False)
    df["I_scaled"]=corr
    df["sigma_scaled"]=ss
    df["robust_weight"]=robust
    df["merge_mean"]=means[codes]
    df["scale_factor"]=factors
    df.to_csv(OUT/f"{output}_observations.csv",index=False)
    ri=np.abs(corr-means[codes]).sum()/np.maximum(corr,0).sum()
    accepted=robust>.5
    ri_clean=np.abs(corr[accepted]-means[codes[accepted]]).sum()/np.maximum(corr[accepted],0).sum()
    rp=np.sqrt(np.sum((corr-means[codes])**2)/np.sum(corr**2))
    print("scales",scales,"Rint",ri,"Rint_filtered",ri_clean,"rejected",np.sum(~accepted),
          "Rrms",rp,"unique",len(m),"observed",np.sum(m.I>2*m.sigma))
    print("absences",m[absences].to_string(index=False))
    report={"scales":scales.tolist(),"Rint_all":float(ri),"Rint_filtered":float(ri_clean),
            "rejected_robust":int(np.sum(~accepted)),"Rrms":float(rp),"observations":len(df),
            "unique":len(m),"unique_allowed":int((~absences).sum()),
            "I_gt_2sigma":int((m.I>2*m.sigma).sum())}
    report["scan_polynomial"]=scan
    report["scale_coefficients"]=beta.tolist()
    if validate_scales:
        valid=validation_codes[codes]
        vc=valid&accepted
        report["scale_validation_unique"]=int(validation_codes.sum())
        report["scale_validation_Rint_all"]=float(np.abs(corr[valid]-means[codes[valid]]).sum()/np.maximum(corr[valid],0).sum())
        report["scale_validation_Rint_filtered"]=float(np.abs(corr[vc]-means[codes[vc]]).sum()/np.maximum(corr[vc],0).sum())
        print("scale validation",report["scale_validation_Rint_all"],report["scale_validation_Rint_filtered"],flush=True)
    (OUT/f"{output}_statistics.json").write_text(json.dumps(report,indent=2))
    return m,df


if __name__=="__main__":
    merge()
