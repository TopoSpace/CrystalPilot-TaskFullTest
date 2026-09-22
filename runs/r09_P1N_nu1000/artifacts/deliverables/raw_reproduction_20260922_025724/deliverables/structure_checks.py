"""Geometry/connectivity audits, without crystallographic libraries."""
from anisotropic import *
from restraints import Geometry
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from collections import Counter
import sys

def check(report,tag):
    atoms=report['atoms'];xs=[];labels=[];els=[];owners=[]
    for i,a in enumerate(atoms):
        ex=expand_xyz(a['xyz'])
        xs.extend(ex);labels.extend([a['label']]*len(ex));els.extend([a['element']]*len(ex));owners.extend([i]*len(ex))
    xs=np.array(xs);els=np.array(els);N=len(xs)
    pairs=[];zrpair=[];contacts=[]
    shifts=np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1],indexing='ij')).reshape(3,-1).T
    for i in range(N):
        d=xs-xs[i]
        dd=d[:,None,:]+shifts[None,:,:]
        ds=np.linalg.norm(dd@A.T,axis=-1);jj=np.argmin(ds,axis=1);dist=ds[np.arange(N),jj]
        for j in np.where((dist>.2)&(dist<4.1))[0]:
            if j<i:continue
            e1,e2=els[i],els[j]
            if e1==e2=='Zr':zrpair.append((i,j,float(dist[j]),shifts[jj[j]].tolist()))
            threshold=2.7 if 'Zr' in [e1,e2] else (1.7 if 'O' in [e1,e2] else 1.8)
            if dist[j]<threshold:pairs.append((i,j,float(dist[j]),shifts[jj[j]].tolist()))
    def components(edges,which):
        ix=np.where(which)[0];lookup={int(j):k for k,j in enumerate(ix)};row=[];col=[]
        for i,j,*rest in edges:
            if i in lookup and j in lookup:row.extend([lookup[i],lookup[j]]);col.extend([lookup[j],lookup[i]])
        mat=coo_matrix((np.ones(len(row)),(row,col)),shape=(len(ix),len(ix)))
        n,lab=connected_components(mat,directed=False)
        return [ix[lab==l].tolist() for l in range(n)]
    zrclusters=components(zrpair,els=='Zr')
    organic=components(pairs,els!='Zr')
    ligand=[g for g in organic if any(els[j]=='C' for j in g)]
    print('Expanded inventory',Counter(els),'Zr clusters',[len(g) for g in zrclusters],flush=True)
    print('Ligands',[dict(Counter(els[g])) for g in ligand],flush=True)
    # Periodic dimensionality: cycle translations of covalent/coordination graph.
    adj=[[] for _ in range(N)]
    for i,j,dd,t in pairs:adj[i].append((j,np.array(t)));adj[j].append((i,-np.array(t)))
    assigned={};cycles=[]
    for root in range(N):
        if root in assigned:continue
        assigned[root]=np.zeros(3,dtype=int);todo=[root]
        while todo:
            i=todo.pop()
            for j,t in adj[i]:
                tr=assigned[i]+t
                if j not in assigned:assigned[j]=tr;todo.append(j)
                else:
                    c=tr-assigned[j]
                    if np.any(c):cycles.append(c.tolist())
    rank=int(np.linalg.matrix_rank(cycles)) if cycles else 0
    print('Periodic framework graph translation rank',rank,flush=True)
    # Asymmetric unit contacts, and local atom valences.
    for i,a in enumerate(atoms):
        n=owners.index(i);ns=[(labels[j],float(np.linalg.norm(A@(xs[j]+t-xs[n])))) for j,t in adj[n]]
        ns.sort(key=lambda t:t[1]);contacts.append({'label':a['label'],'element':a['element'],'neighbors':ns,'Uiso':a['Uiso'],'U_eigenvalues':a.get('U_eigenvalues')})
        print(a['label'],[(l,round(d,4)) for l,d in ns],flush=True)
    seed=json.loads((OUT/'framework_iso1.json').read_text())['atoms'];dummy=AnisoModel(atoms,np.array([[1,0,0]]),False);geom=Geometry(dummy,seed);p=dummy.p0
    distances=[]
    xx=geom.xyz(p)
    for i,j,r,t,target,sd in geom.bonds:
        dist=float(np.linalg.norm(A@(r@xx[j]+t-xx[i])));distances.append({'a':atoms[i]['label'],'b':atoms[j]['label'],'distance':dist,'target':target,'sigma_restraint':sd,'z':(dist-target)/sd})
    planes=[]
    for refs in geom.planes:
        xyz=np.array([A@(r@xx[j]+t) for j,r,t in refs]);xyz-=xyz.mean(axis=0);_,_,vh=np.linalg.svd(xyz,full_matrices=False)
        dev=xyz@vh[-1];planes.append({'atom_count':len(refs),'rms_A':float(np.sqrt(np.mean(dev**2))),'max_A':float(max(abs(dev)))})
    print('Plane deviations',planes,flush=True)
    print('Bond restraints',[(r['a'],r['b'],round(r['distance'],4),round(r['z'],2)) for r in distances],flush=True)
    reportout={'inventory':dict(Counter(els.tolist())),'node_sizes':[len(g) for g in zrclusters],'ligand_inventories':[dict(Counter(els[g].tolist())) for g in ligand],
        'periodic_connectivity_rank':rank,'contacts':contacts,'restrained_bonds':distances,'plane_deviations':planes,
        'min_ADP_eigenvalue':float(min(min(a.get('U_eigenvalues',[a['Uiso']])) for a in atoms)),
        'max_ADP_eigenvalue':float(max(max(a.get('U_eigenvalues',[a['Uiso']])) for a in atoms))}
    (OUT/(tag+'_validation.json')).write_text(json.dumps(reportout,indent=2))
    return reportout

if __name__=='__main__':
    tag=sys.argv[1] if len(sys.argv)>1 else 'framework_restrained_disp'
    check(json.loads((OUT/(tag+'.json')).read_text()),tag)
