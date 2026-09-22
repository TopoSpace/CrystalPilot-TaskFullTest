"""Task-local, independently written Oxford TY6 decoder and peak finder.
Format inferred from supplied byte streams, not copied from a format library.
Each row: one literal difference, 49 packed blocks of 16, 15 literal
values. Each packed block control contains two widths (low nibble first).
Width w codes eight biased signed differences; width 8 uses 254/255
escapes whose int16/int32 values follow that block's two packed halves.
"""
from pathlib import Path
import re, struct, json, hashlib, time, os
import numpy as np
from scipy import ndimage
ROOT=Path(__file__).resolve().parents[1]
INP=ROOT/'inputs'/'alanine'; OUT=Path(os.environ.get('ALANINE_OUTPUT',str(ROOT/'deliverables'))).resolve()
assert OUT.is_relative_to(ROOT/'deliverables'), 'Outputs must stay within this task deliverables'
TAB1=[tuple((v>>j)&1 for j in range(8)) for v in range(256)]
TAB2=[tuple(((v>>(2*j))&3)-1 for j in range(8)) for v in range(65536)]

def read_frame(path, pixels=True):
    b=Path(path).read_bytes()
    asciihead=b[:256].decode('ascii',errors='replace')
    nx,ny,oi,ol=map(int,re.search(r'NX=\s*(\d+) NY=\s*(\d+) OI=\s*(\d+) OL=\s*(\d+)',asciihead).groups())
    nh=int(re.search(r'NHEADER=\s*(\d+)',asciihead)[1])
    mot0=np.array(struct.unpack_from('<5i',b,1820))/[12800,12800,6400,6400,12800]
    mot1=np.array(struct.unpack_from('<5i',b,1860))/[12800,12800,6400,6400,12800]
    meta=dict(file=Path(path).name,nx=nx,ny=ny,omega0=mot0[0],omega1=mot1[0],theta=mot0[1],kappa=mot0[2],phi=mot0[3],distance=mot0[4],exposure=2.0,wavelength=struct.unpack_from('<d',b,2120)[0])
    if not pixels:return meta
    assert 'COMPRESSION=TY6' in asciihead
    n=struct.unpack_from('<I',b,nh)[0];pos=nh+4;end=pos+n
    arr=np.empty((ny,nx),np.int32);overflow16=overflow32=0
    def literal(pos):
        nonlocal overflow16,overflow32
        c=b[pos];pos+=1
        if c==254:
            c=struct.unpack_from('<h',b,pos)[0];pos+=2;overflow16+=1
        elif c==255:
            c=struct.unpack_from('<i',b,pos)[0];pos+=4;overflow32+=1
        else:c-=127
        return c,pos
    for y in range(ny):
        first,pos=literal(pos);vals=[first]
        for block in range((nx-1)//16):
            control=b[pos];pos+=1;v16=[]
            for w in (control&15,control>>4):
                assert 1<=w<=8,(path,y,pos,w)
                if w==1:v8=TAB1[b[pos]]
                elif w==2:v8=TAB2[b[pos]|(b[pos+1]<<8)]
                else:
                    v=int.from_bytes(b[pos:pos+w],'little');mask=(1<<w)-1;bias=(1<<(w-1))-1
                    v8=tuple(((v>>(j*w))&mask)-bias for j in range(8))
                pos+=w;v16.extend(v8)
            for j,c in enumerate(v16):
                if c==127:
                    v16[j]=struct.unpack_from('<h',b,pos)[0];pos+=2;overflow16+=1
                elif c==128:
                    v16[j]=struct.unpack_from('<i',b,pos)[0];pos+=4;overflow32+=1
            vals.extend(v16)
        for j in range((nx-1)%16):
            c,pos=literal(pos);vals.append(c)
        arr[y]=np.cumsum(vals)
    assert pos==end,(path,pos,end)
    assert overflow16==oi and overflow32==ol,(path,overflow16,oi,overflow32,ol)
    assert arr.min()>=-1,(path,'negative pixels',arr.min())
    meta.update(max_pixel=int(arr.max()),sum_pixel=int(arr.sum()),overflow16=oi,overflow32=ol)
    return arr,meta

def peaks_from_image(im, threshold=12):
    a=np.maximum(im,0).astype(float)
    # Local 3x3 signal, compact maxima; weighted centroids in a 7x7 box.
    sm=ndimage.uniform_filter(a,3)*9
    local=ndimage.maximum_filter(sm,size=7)
    yx=np.argwhere((sm==local)&(sm>threshold)&(a>1))
    out=[]
    for y,x in yx:
        if x<6 or x>793 or y<6 or y>768 or 379<x<421:continue
        box=a[y-4:y+5,x-4:x+5];bg=np.median(box)
        yy,xx=np.mgrid[y-4:y+5,x-4:x+5]
        rr=(yy-y)**2+(xx-x)**2
        z=np.maximum(box-bg,0)*(rr<=12)
        total=z.sum()
        if total<25:continue
        cy=(z*yy).sum()/total;cx=(z*xx).sum()/total
        varx=(z*(xx-cx)**2).sum()/total;vary=(z*(yy-cy)**2).sum()/total
        if max(varx,vary)>7:continue
        out.append([cx,cy,total,varx,vary,float(bg)])
    return np.array(out).reshape(-1,6)

def process_frame(path):
    im,meta=read_frame(path)
    peaks=peaks_from_image(im)
    return meta,peaks,hashlib.sha256(Path(path).read_bytes()).hexdigest()

def extract():
    from concurrent.futures import ProcessPoolExecutor
    files=sorted((INP/'frames').glob('pgw240033_Mo_[1-6]_*.rodhypix'),key=lambda p:tuple(map(int,p.stem.split('_')[-2:])))
    t=time.time(); rows=[];meta=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for idx,(m,p,sha) in enumerate(pool.map(process_frame,files,chunksize=4)):
            m['sha256']=sha;m['npeaks']=len(p);meta.append(m)
            run,frame=map(int,files[idx].stem.split('_')[-2:])
            for pk in p:rows.append([run,frame,idx,*pk])
            if idx%50==0:print('decoded',idx+1,'peaks',len(rows),'seconds',round(time.time()-t,1),flush=True)
    np.save(OUT/'peaks.npy',np.array(rows))
    (OUT/'frame_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
    print('DONE',len(files),len(rows),'sec',time.time()-t,flush=True)

if __name__=='__main__':
    import sys
    if '--extract' in sys.argv:extract()
    else:
        for run,frame in [(1,1),(1,92),(2,1),(3,1),(4,1),(5,1),(6,1)]:
            p=INP/'frames'/f'pgw240033_Mo_{run}_{frame}.rodhypix'
            im,m=read_frame(p);print(m);print('quantiles',np.percentile(im,[0,50,90,99,99.9,100]),'peaks',peaks_from_image(im)[:10])
