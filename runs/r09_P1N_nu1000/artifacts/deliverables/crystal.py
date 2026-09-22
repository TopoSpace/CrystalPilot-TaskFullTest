"""First-principles crystallographic utilities; only Python stdlib, numpy, scipy.
Written specifically for this run, no crystallographic software dependencies.
Fractional coordinates are rows, direct Cartesian vectors are cell matrix columns.
"""
from pathlib import Path
import json, re, time
import numpy as np
from scipy.ndimage import maximum_filter

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'deliverables'
A=np.array([[39.19,-39.19/2,0],[0,39.19*np.sqrt(3)/2,0],[0,0,16.61]])
G=A.T@A
GI=np.linalg.inv(G)
V=abs(np.linalg.det(A))
# Cromer-Mann neutral-atom four-Gaussian coefficients (analytical approximation).
CM={'C':([2.31,1.02,1.5886,.865],[20.8439,10.2075,.5687,51.6512],.2156),
    'O':([3.0485,2.2868,1.5463,.867],[13.2771,5.7011,.3239,32.9089],.2508),
    'Zr':([17.8765,10.948,5.41732,3.65721],[1.27618,11.916,.117622,87.6627],2.06929),
    'H':([.493002,.322912,.140191,.04081],[10.5109,26.1257,3.14236,57.7997],.003038)}

def rotations():
    lines=(ROOT/'inputs'/'start.ins').read_text().splitlines()
    expr=['x,y,z']+[s[5:].strip().lower() for s in lines if s.startswith('SYMM ')]
    mats=[]
    for e in expr:
        rows=[]
        for term in e.split(','):
            row=[]
            for v in 'xyz':
                m=re.search(r'([+-]?)'+v,term)
                row.append(0 if m is None else (-1 if m[1]=='-' else 1))
            rows.append(row)
        r=np.array(rows,int)
        mats += [r,-r]
    return np.unique(np.array(mats),axis=0)
R=rotations()

def equiv_h(h,rots=R):
    return np.einsum('ni,sij->nsj',h,rots)

def canonical(h,rots=R):
    eq=equiv_h(h,rots)
    # lexicographic maximum, with generously large integer radix
    key=(eq[:,:,0]+512)*1048576+(eq[:,:,1]+512)*1024+(eq[:,:,2]+512)
    return eq[np.arange(len(h)),np.argmax(key,axis=1)]

def expand_xyz(x,rots=R):
    e=np.einsum('sij,j->si',rots,np.array(x))%1
    return np.unique(np.round(e,9),axis=0)

def distvec(x,y):
    d=np.array(x)-np.array(y)
    d-=np.round(d)
    return d@A.T

def mindist(x,ys):
    # neighbors sufficient for hexagonal metric with wrapped differences
    d=np.array(ys)-np.array(x)
    d-=np.round(d)
    vals=[]
    for u,v in [(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]:
        dd=d+np.array([u,v,0])
        vals.append(np.linalg.norm(dd@A.T,axis=-1))
    return np.min(vals,axis=0)

def formfactor(el,s2):
    a,b,c=CM[el]
    return np.exp(-np.array(s2)[:,None]*np.array(b))@np.array(a)+c

def read_raw():
    data=[]
    for line in (ROOT/'inputs'/'start.hkl').read_text().splitlines():
        if len(line)<28: continue
        try: row=[int(line[:4]),int(line[4:8]),int(line[8:12]),float(line[12:20]),float(line[20:28])]
        except ValueError: continue
        if row[:3]==[0,0,0]: break
        data.append(row)
    d=np.array(data)
    return d[:,:3].astype(int),d[:,3],d[:,4]

def merge(h,I,sig,rots=R):
    hc=canonical(h,rots)
    hu,idx,n=np.unique(hc,axis=0,return_inverse=True,return_counts=True)
    w=1/np.maximum(sig,1e-4)**2
    sw=np.bincount(idx,weights=w)
    im=np.bincount(idx,weights=I*w)/sw
    se=1/np.sqrt(sw)
    # Internal scatter uncertainty of weighted mean, not lower than counting estimate.
    chi=np.bincount(idx,weights=w*(I-im[idx])**2)
    se*=np.sqrt(np.maximum(1,chi/np.maximum(n-1,1)))
    rint=np.sum(abs(I-im[idx]))/np.sum(abs(I))
    return hu,im,se,n,idx,rint

def full_reflections(h,val):
    e=equiv_h(h).reshape(-1,3)
    v=np.repeat(val,len(R))
    eh,ix=np.unique(e,axis=0,return_index=True)
    return eh,v[ix]

def map_fft(h,val,shape=None):
    eh,ev=full_reflections(h,val)
    if shape is None:
        from scipy.fft import next_fast_len
        shape=tuple(next_fast_len(int(2*m+2)) for m in np.max(abs(eh),axis=0))
    grid=np.zeros(shape,complex)
    ij=tuple((eh%np.array(shape)).T)
    grid[ij]=ev
    rho=np.fft.fftn(grid).real/V
    return rho

def peaks(rho,n=100,min_dist=0.8,unique=True):
    shape=np.array(rho.shape)
    loc=np.argwhere((rho==maximum_filter(rho,size=3,mode='wrap'))&(rho>0))
    order=np.argsort(rho[tuple(loc.T)])[::-1]
    ans=[]; expanded=[]
    for ij in loc[order]:
        x=ij/shape
        if expanded and np.min(mindist(x,np.array(expanded)))<min_dist: continue
        ans.append((x,float(rho[tuple(ij)])))
        expanded.extend(expand_xyz(x) if unique else [x])
        if len(ans)>=n: break
    return ans

def save_peaks(p,name):
    (OUT/name).write_text('x,y,z,height,multiplicity\n'+'\n'.join(','.join(map(str,[*x,v,len(expand_xyz(x))])) for x,v in p))

def load_data():
    d=np.load(OUT/'merged.npz')
    return d['h'],d['I'],d['sig']

def audit():
    h,I,sig=read_raw()
    hu,im,se,n,idx,rint=merge(h,I,sig)
    q=np.einsum('ni,ij,nj->n',hu,GI,hu)
    np.savez_compressed(OUT/'merged.npz',h=hu,I=im,sig=se,multiplicity=n,q=q)
    report={'cell':[39.19,39.19,16.61,90,90,120],'volume_A3':V,'symmetry_operators':R.tolist(),
        'raw_reflections':len(h),'unique_P6mmm':len(hu),'r_int_weighted_center':rint,
        'index_min':h.min(axis=0).tolist(),'index_max':h.max(axis=0).tolist(),
        'd_min':float(1/np.sqrt(max(q))),'I_range':[float(min(I)),float(max(I))],
        'negative_fraction':float(np.mean(I<0)),'unique_I_gt_2sig':int(np.sum(im>2*se)),
        'redundancy':float(np.mean(n))}
    # Compare symmetry components by merge residual and correlation.
    rep=[]
    for r in R:
        rr=np.array([np.eye(3,dtype=int),r])
        *_,ri=merge(h,I,sig,rr)
        rep.append({'op':r.tolist(),'rint':ri})
    report['individual_symmetry_residuals']=rep
    bins=np.linspace(0,np.sqrt(max(q)),11)
    report['resolution_bins']=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        mask=(np.sqrt(q)>lo)&(np.sqrt(q)<=hi)
        report['resolution_bins'].append({'inv_d_lo':lo,'inv_d_hi':hi,'n':int(mask.sum()),'mean_I_sig':float(np.mean(im[mask]/se[mask])),'fraction_2sig':float(np.mean(im[mask]>2*se[mask]))})
    (OUT/'audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['symmetry_operators','individual_symmetry_residuals']},indent=2))
    print('Symmetry Ri:',[(x['op'],round(x['rint'],4)) for x in rep])
    rho=map_fft(hu,np.maximum(im,0),shape=(192,192,96))
    save_peaks(peaks(rho,n=80,min_dist=.7),'patterson_peaks.csv')
    print('Top Patterson peaks:')
    for x,v in peaks(rho,n=25,min_dist=.7):print(np.round(x,6),round(v,4))

if __name__=='__main__':audit()
