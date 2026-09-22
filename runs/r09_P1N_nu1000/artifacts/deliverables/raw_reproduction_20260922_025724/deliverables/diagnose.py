from crystal import *
h,I,sig=read_raw()
hu,idx,n=np.unique(h,axis=0,return_inverse=True,return_counts=True)
iu=np.bincount(idx,weights=I)/n
su=np.sqrt(np.bincount(idx,weights=sig**2))/n
D={tuple(x):i for i,x in enumerate(hu)}
res=[]
for r in R:
    eq=hu@r
    a=[];b=[]
    for i,x in enumerate(eq):
        j=D.get(tuple(x))
        if j is not None:a.append(i);b.append(j)
    a=np.array(a);b=np.array(b)
    strong=(iu[a]>5*su[a])&(iu[b]>5*su[b])
    q=np.einsum('ni,ij,nj->n',hu[a],GI,hu[a])
    low=q<.5**2
    rs=float(np.sum(abs(iu[a[strong]]-iu[b[strong]]))/np.sum(abs(iu[a[strong]]+iu[b[strong]])))
    rl=float(np.sum(abs(iu[a[low]]-iu[b[low]]))/np.sum(abs(iu[a[low]]+iu[b[low]])))
    res.append({'r':r.tolist(),'pairs':len(a),'strong_n':int(strong.sum()),'corr_strong':float(np.corrcoef(iu[a[strong]],iu[b[strong]])[0,1]),'R_strong':rs,'R_low':rl})
print(json.dumps(res,indent=1))
# Equal-weight means avoid observation-dependent variance bias.
hc=canonical(h)
hm,ix,n=np.unique(hc,axis=0,return_inverse=True,return_counts=True)
im=np.bincount(ix,weights=I)/n
sigmean=np.sqrt(np.bincount(ix,weights=sig**2))/n
scatter=np.sqrt(np.bincount(ix,weights=(I-im[ix])**2)/np.maximum(n-1,1)/n)
se=np.maximum(sigmean,scatter)
q=np.einsum('ni,ij,nj->n',hm,GI,hm)
np.savez_compressed(OUT/'merged_equal.npz',h=hm,I=im,sig=se,multiplicity=n,q=q)
report={'symmetry':res,'bins':[]}
for dl,dh in [(0,0.1),(.1,.2),(.2,.4),(.4,.6),(.6,.8),(.8,1),(1,1.2),(1.2,1.5)]:
    m=(q>dl**2)&(q<=dh**2); mr=m[ix]
    rr={'inv_d':[dl,dh],'n':int(m.sum()),'mean_I_sig':float(np.mean(im[m]/se[m])),'Rint':float(np.sum(abs(I[mr]-im[ix[mr]]))/np.sum(abs(I[mr])))}
    report['bins'].append(rr)
print('EQUAL MERGE',json.dumps(report['bins'],indent=2))
(OUT/'diagnostics.json').write_text(json.dumps(report,indent=2))
