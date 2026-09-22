from refine import *

z=json.loads((OUT/'zr_refined.json').read_text());atoms=z['atoms']
pe=np.loadtxt(OUT/'zr_refined_diff_peaks.csv',delimiter=',',skiprows=1)
for k,j in enumerate(range(2,8)):
    xyz=representatives_near(pe[j,:3])
    atoms.append(atom('O'+str(k+1),'O',xyz,.055,tolerance=.14))
for k,j in enumerate(range(8,20)):
    xyz=representatives_near(pe[j,:3])
    atoms.append(atom('C'+str(k+1),'C',xyz,.055,tolerance=.14))
(OUT/'framework_seed.json').write_text(json.dumps({'atoms':atoms,'scale':z['scale'],'description':'Zr positions plus six O and twelve C peaks from heavy-atom difference Fourier; assignments checked against connectivity.'},indent=2))
r,fc=refine(atoms,z['scale'],tag='framework_iso1',dmin=1.,dmax=20.,nfev=150,aweight=.08)
maps(r,fc,'framework_iso1')
inspect_peaks('framework_iso1_diff_peaks.csv',r['atoms'],40)
# All unique atom-orbit contacts at the representative atom.
contacts=[]
for a in r['atoms']:
    near=[]
    for b in r['atoms']:
        eq=expand_xyz(b['xyz']);ds=mindist(a['xyz'],eq)
        for j in np.where((ds>.2)&(ds<(2.7 if 'Zr' in [a['element'],b['element']] else 1.85)))[0]:
            near.append((b['label'],float(ds[j]),eq[j].tolist()))
    near.sort(key=lambda x:x[1]);contacts.append({'atom':a['label'],'neighbors':near})
    print('BONDS',a['label'],[(b,round(d,3)) for b,d,x in near],flush=True)
(OUT/'framework_iso1_contacts.json').write_text(json.dumps(contacts,indent=2))
