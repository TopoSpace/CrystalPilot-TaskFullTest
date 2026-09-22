from pathlib import Path
import struct, collections, numpy as np
R=Path(__file__).resolve().parents[1]; O=R/'deliverables'
b=(R/'inputs/alanine/frames/pgw240033_Mo_1_1.rodhypix').read_bytes(); n=struct.unpack_from('<I',b,6576)[0]
print('length',len(b),'blocklen',n,'remaining',len(b)-6576-4-n)
for off in [6576+n-16,6576+n,6576+n+4,len(b)-128]:
 print(off,b[off:off+128].hex(' '),np.frombuffer(b[off:off+128],'<i4'))
s=b[6580:6580+n]
for maxw in [3,4,5,6]:
 vals=[]; i=0; ctrl=collections.Counter()
 while i<len(s):
  c=s[i];i+=1;lo=c&15;hi=c>>4
  if 1<=lo<=maxw and 1<=hi<=maxw:
   ctrl[c]+=1
   for w in [lo,hi]:
    v=int.from_bytes(s[i:i+w],'little'); i+=w
    vals.extend([(v>>(j*w)&((1<<w)-1))-((1<<(w-1))-1) for j in range(8)])
  else: vals.append(c-127)
 a=np.array(vals,dtype='int64');v=np.cumsum(a)
 print('maxw',maxw,'consumed',i,'pixels',len(a),'ctrl',ctrl,'sumdiff',a.sum(),'range',v.min(),v.max(),'first100',v[:100])
 if len(a)>=620000:
  im=v[:620000].reshape(775,800); print('quantiles',np.percentile(im,[0,25,50,90,99,100]),'row0',im[0,370:430]);np.save(O/f'decode_candidate_{maxw}.npy',im)
print('Header motor area')
for off in range(1536,1880,4):
 print(off,struct.unpack_from('<i',b,off)[0])
