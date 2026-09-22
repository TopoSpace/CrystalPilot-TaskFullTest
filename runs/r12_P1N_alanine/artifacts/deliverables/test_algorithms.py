"""Task-local synthetic tests and raw-input integrity checks."""
import json,struct,hashlib
import numpy as np
from raw_frames import OUT,INP,read_frame
from geometry import indices
from integrate import aperture

# Independent simple encoder, only for synthetic decoder tests.
def encode_test(image):
 stream=bytearray();oi=ol=0;widths=set()
 def escaped(v):
  nonlocal oi,ol
  if -127<=v<=126:return v+127,b''
  if -32768<=v<=32767:oi+=1;return 254,struct.pack('<h',v)
  ol+=1;return 255,struct.pack('<i',v)
 def literal(v):
  c,b=escaped(v);stream.append(c);stream.extend(b)
 for row in image:
  diff=np.diff(np.r_[0,row]).astype(int);literal(int(diff[0]));pos=1
  for _ in range((len(row)-1)//16):
   block=diff[pos:pos+16];pos+=16;ws=[];bits=[];ext=[]
   for half in [block[:8],block[8:]]:
    w=8
    for test in range(1,8):
     if half.min()>=-(2**(test-1)-1) and half.max()<=2**(test-1):w=test;break
    widths.add(w);ws.append(w)
    if w==8:
     codes=[]
     for v in half:
      code,extra=escaped(int(v));codes.append(code);ext.append(extra)
     bits.append(bytes(codes))
    else:
     value=sum((int(v)+(2**(w-1)-1))<<(w*j) for j,v in enumerate(half));bits.append(value.to_bytes(w,'little'))
   stream.append(ws[0]|(ws[1]<<4));stream.extend(b''.join(bits));stream.extend(b''.join(ext))
  for v in diff[pos:]:literal(int(v))
 ny,nx=image.shape;header=bytearray(6576)
 text=f'OD SAPPHIRE  4.0\r\nCOMPRESSION=TY6(test)\r\nNX={nx} NY={ny} OI={oi} OL={ol}\r\nNHEADER=6576\r\n'.encode()
 header[:len(text)]=text;struct.pack_into('<d',header,2120,.71073)
 return bytes(header)+struct.pack('<I',len(stream))+stream,widths,oi,ol

def main():
 rng=np.random.default_rng(1276);im=rng.poisson(.08,size=(12,800)).astype(np.int32)
 diff=[0]
 for block in range(49):
  for half in range(2):
   w=(block*2+half)%8+1;diff.extend([2**(w-1),-(2**(w-1)-1)]*4)
 diff.extend([0]*15);im[0]=np.cumsum(diff)
 im[1,11]=1000000;im[2,710]=32760;im[3,100]=-1
 blob,widths,oi,ol=encode_test(im);p=OUT/'synthetic_decoder_fixture.rodhypix';p.write_bytes(blob)
 decoded,meta=read_frame(p)
 tests=dict(TY6_synthetic_roundtrip=bool(np.array_equal(decoded,im)),all_eight_bit_widths_exercised=widths==set(range(1,9)),both_escape_types_exercised=bool(oi>0 and ol>0))
 flat=np.ones((775,800),np.int32);val=aperture(flat,200.3,250.7)
 tests['flat_background_gives_zero_net']=abs(val[0])<1e-12
 flat[251,200]+=500;val=aperture(flat,200.3,250.7)
 tests['known_signal_recovered']=abs(val[0]-500)<1e-12
 par=np.array(json.loads((OUT/'geometry.json').read_text())['parameters']);a=np.concatenate([np.load(OUT/f'integrated_run{r}.npy') for r in range(1,7)])
 p=np.zeros((len(a),9));p[:,0]=a[:,0];p[:,3:5]=a[:,5:7]
 meta=json.loads((OUT/'frame_metadata.json').read_text());th={int(m['file'].split('_')[-2]):m['theta'] for m in meta}
 ang=np.zeros((len(a),4));ang[:,0]=a[:,4];ang[:,1]=[th[int(r)] for r in a[:,0]]
 hi=indices(par,p,ang);error=float(np.max(abs(hi-a[:,1:4])))
 tests['prediction_geometry_inverse']=error<1e-9
 changed=[]
 for m in meta:
  sha=hashlib.sha256((INP/'frames'/m['file']).read_bytes()).hexdigest()
  if sha!=m['sha256']:changed.append(m['file'])
 tests['all_616_raw_frame_hashes_unchanged']=not changed and len(meta)==616
 tests={key:bool(value) for key,value in tests.items()}
 result=dict(tests=tests,all_pass=all(tests.values()),n_raw_frames=len(meta),changed_inputs=changed,synthetic_widths=sorted(widths),synthetic_int16_escapes=oi,synthetic_int32_escapes=ol,max_prediction_index_error=error)
 (OUT/'algorithm_tests.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
 if not all(tests.values()):raise RuntimeError('Algorithm test failed')
if __name__=='__main__':main()
