"""Inspect only the supplied raw data; no third-party format code is used."""
from pathlib import Path
import struct
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
p = ROOT / "inputs/alanine/frames/pgw240033_Mo_1_1.rodhypix"
b = p.read_bytes()
for start, stop in [(6400, 6740), (len(b)-400, len(b))]:
    for i in range(start, stop, 16):
        print(f"{i:7d}: {b[i:i+16].hex(' ')}")
print("Header nonzero double values (8-byte aligned within sections)")
for start, stop in [(256,768),(768,1536),(1536,2560),(2560,3072)]:
    for i in range(start, stop, 8):
        v = struct.unpack_from("<d", b, i)[0]
        if 1e-6 < abs(v) < 1e7:
            print(i, v)
print("At data start int32", struct.unpack_from("<10I",b,6576))
print("Length",len(b))
