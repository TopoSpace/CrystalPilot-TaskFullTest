from crystal import *
h,I,sig=read_raw(); hc=canonical(h)
hu,ix,n=np.unique(hc,axis=0,return_inverse=True,return_counts=True)
rng=np.random.default_rng(20260922); split=rng.integers(0,2,len(I))
ni=[];ii=[]
for j in [0,1]:
    m=split==j
    nn=np.bincount(ix[m],minlength=len(hu)); ss=np.bincount(ix[m],weights=I[m],minlength=len(hu))
    ni.append(nn);ii.append(ss/np.maximum(nn,1))
D=np.load(OUT/'merged_equal.npz'); im=D['I'];se=D['sig'];q=D['q']; invd=np.sqrt(q)
lim=float(max(invd)); limh=int(np.ceil(39.19*lim))+1;liml=int(np.ceil(16.61*lim))+1
ha=np.array(np.meshgrid(np.arange(-limh,limh+1),np.arange(-limh,limh+1),np.arange(-liml,liml+1),indexing='ij')).reshape(3,-1).T
qa=np.einsum('ni,ij,nj->n',ha,GI,ha)
ha=ha[(qa>0)&(qa<=lim**2+1e-9)]
ht=np.unique(canonical(ha),axis=0);qt=np.einsum('ni,ij,nj->n',ht,GI,ht)
print('total completeness',len(hu),len(ht),len(hu)/len(ht),flush=True)
rows=[]
edges=[0,1/6,1/4,1/3,1/2.5,1/2,1/1.8,1/1.5,1/1.3,1/1.2,1/1.1,1,1/.9,1/.8,1/.75,lim+1e-8]
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(invd>lo)&(invd<=hi);mm=m&(ni[0]>=2)&(ni[1]>=2);mr=m[ix]
    cc=float(np.corrcoef(ii[0][mm],ii[1][mm])[0,1])
    nt=int(np.sum((qt>lo**2)&(qt<=hi**2)))
    row={'d_low':float(1/lo) if lo else None,'d_high':float(1/hi),'observed':int(m.sum()),'possible':nt,'completeness':float(m.sum()/nt),'I_sigma':float(np.mean(im[m]/se[m])),'CC_half':cc,'Rint':float(np.sum(abs(I[mr]-im[ix[mr]]))/np.sum(abs(I[mr]))),'n_half_pairs':int(mm.sum())}
    rows.append(row)
    print(row,flush=True)
(OUT/'data_quality.json').write_text(json.dumps({'total_observed':len(hu),'total_possible':len(ht),'total_completeness':len(hu)/len(ht),'shells':rows},indent=2))
