"""Quantify detector-background attenuation near the suspect measurements."""
import csv
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter, uniform_filter
from run_command import ROOT
work=ROOT/'work/final_data'
rows=list(csv.DictReader((work/'suspicious_reflections_raw_pixels.csv').open()))
gmap={0:0,1:0,2:1,3:2,4:3,5:3}
arrays={g:np.load(work/f'shadow_diagnostic_group_{g}.npz') for g in range(4)}
for row in rows:
    if float(row['I_prf'])>30 or int(row['panel'])!=0:
        continue
    group=gmap[int(row['exp'])];s=arrays[group]['sum0'];x=int(float(row['x']));y=int(float(row['y']))
    vals=[]
    for dy in [-70,-40,-20,0,20,40,70]:
        yy=y+dy
        a=s[max(yy-5,0):yy+6,max(x-5,0):x+6]
        vals.append(round(float(np.median(a)),2))
    print('group',group,'xy',x,y,'hkl',row['hkl'],'median accumulated counts [y-70,-40,-20,y,+20,+40,+70]:',vals)
for g,a in arrays.items():
    s=a['sum0']
    print('\nGroup',g,'median accumulated counts in 10x10 patches; rows y, columns x=50,100,150,200,250,300,325,350,375')
    for y in range(320,451,10):
        print(y,[round(float(np.median(s[y-5:y+5,x-5:min(x+5,385)])),1) for x in [50,100,150,200,250,300,325,350,375]])
