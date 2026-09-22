"""Inspect ONLY this task's input files; no third-party format reader used."""
from pathlib import Path
import struct, collections, json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
INP=ROOT/'inputs'/'alanine'
OUT=ROOT/'deliverables'

def main():
    paths=list(INP.rglob('*'))
    counts=collections.Counter(p.suffix for p in paths if p.is_file())
    print('file counts',dict(counts))
    files=sorted((INP/'frames').glob('*.rodhypix'))
    print('frames',len(files))
    p=INP/'frames'/'pgw240033_Mo_1_1.rodhypix'
    b=p.read_bytes()
    print('FRAME',p.name,'size',len(b))
    print('FIRST 4096 repr',repr(b[:4096]))
    print('FIRST 4096 ascii',b[:4096].decode('latin1',errors='replace'))
    for off in range(0,min(8192,len(b)),16):
        chunk=b[off:off+16]
        print(f'{off:05d}',chunk.hex(' '), repr(chunk), struct.unpack('<4i',chunk),struct.unpack('<2d',chunk))
    for dtype in ['<i4','<u2','<f4']:
        if len(b)>=800*775*np.dtype(dtype).itemsize:
            arr=np.frombuffer(b[-800*775*np.dtype(dtype).itemsize:],dtype=dtype)
            print('tail',dtype,np.percentile(arr,[0,25,50,90,99,99.9,100]),'nonzero',np.count_nonzero(arr))
    inv=[{'path':str(p.relative_to(ROOT)), 'bytes':p.stat().st_size} for p in paths if p.is_file()]
    (OUT/'input_inventory.json').write_text(json.dumps(inv,indent=2),encoding='utf8')

if __name__=='__main__':
    import contextlib
    with (OUT/'probe_inputs.log').open('w',encoding='utf8') as f,contextlib.redirect_stdout(f):main()
