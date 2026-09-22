"""Self-written decoder inferred and checked against the supplied TY6 frames."""
from pathlib import Path
import re
import struct
import numpy as np


def read_frame(path):
    b = Path(path).read_bytes()
    text = b[:256].decode("ascii", errors="replace")
    nx, ny = map(int, re.search(r"NX=\s*(\d+) NY=\s*(\d+)", text).groups())
    nh = int(re.search(r"NHEADER=\s*(\d+)", text)[1])
    ncomp = struct.unpack_from("<I", b, nh)[0]
    start = nh + 4
    ends = np.frombuffer(b, "<u4", ny, start+ncomp)
    image = np.empty((ny,nx),np.int32)
    widths_hist = np.zeros(16,int)
    for y in range(ny):
        p = start + int(ends[y])
        ds = []
        # Each row: one ordinary difference, blocks of 16, remaining 15.
        v = b[p]
        p += 1
        ds.append(v-127)
        for _ in range((nx-1)//16):
            control = b[p]
            p += 1
            for bits in (control & 15, control >> 4):
                widths_hist[bits] += 1
                if not 1 <= bits <= 8:
                    raise ValueError(f"Unsupported width {bits}, row {y}, offset {p}")
                w = int.from_bytes(b[p:p+bits],"little")
                p += bits
                mask = (1<<bits)-1
                bias = (1<<(bits-1))-1
                ds.extend(((w>>(bits*j))&mask)-bias for j in range(8))
            for j in range(len(ds)-16,len(ds)):
                if ds[j] in (127,128):
                    n = 2 if ds[j] == 127 else 4
                    ds[j] = int.from_bytes(b[p:p+n],"little",signed=True)
                    p += n
        tailstart = len(ds)
        for _ in range((nx-1)%16):
            ds.append(b[p]-127)
            p += 1
        for j in range(tailstart,len(ds)):
            if ds[j] in (127,128):
                n = 2 if ds[j] == 127 else 4
                ds[j] = int.from_bytes(b[p:p+n],"little",signed=True)
                p += n
        expected = start + (int(ends[y+1]) if y+1<ny else ncomp)
        if p != expected:
            raise ValueError(f"Row {y}: consumed to {p}, expected {expected}")
        image[y] = np.cumsum(ds)
    return image, {"text":text,"widths":widths_hist.tolist(),
                   "angles_start_steps":struct.unpack_from("<5i",b,1820),
                   "angles_end_steps":struct.unpack_from("<5i",b,1860)}


if __name__ == "__main__":
    import sys
    a, meta = read_frame(sys.argv[1])
    print(meta)
    print(a.min(),a.max(),a.sum(),np.quantile(a,[0,.1,.5,.9,.99,.999,1]))
    print("row values",a[0])
