"""Transparent chemical restraints derived from the solved atom connectivity.
Targets are generic aromatic/carboxylate chemistry, not a reference structure.
"""
from anisotropic import *
from itertools import combinations

class Geometry:
    def __init__(self,model,seed_atoms,strength=1.):
        self.model=model;self.seed=seed_atoms;self.bonds=[];self.angles=[];self.planes=[];self.strength=strength
        self.indices=[j for s in model.spec for j in range(s['off'],s['off']+s['nd'])]
        self.lab={a['label']:i for i,a in enumerate(seed_atoms)}
        self.ns=[]
        for i,a in enumerate(seed_atoms):
            ns=[]
            for j,b in enumerate(seed_atoms):
                if b['element']=='Zr':continue
                for r in np.array(b['reps']):
                    x=r@np.array(b['xyz']);shift=np.rint(np.array(a['xyz'])-x);x+=shift
                    d=np.linalg.norm(A@(x-a['xyz']))
                    if .3<d<1.8:ns.append((j,r,shift,d))
            self.ns.append(ns)
            if a['element']!='C':continue
            for j,r,t,d in ns:
                b=seed_atoms[j]
                if i>j and b['element']=='C':continue
                la,lb=a['label'],b['label']
                if b['element']=='O':target=1.26;sd=.025
                elif set([la,lb])=={'C10','C12'}:target=1.49;sd=.03
                elif set([la,lb])=={'C5','C6'}:target=1.49;sd=.035
                elif la in ['C2','C3','C4','C6','C8','C10'] and lb in ['C2','C3','C4','C6','C8','C10']:target=1.395;sd=.025
                else:target=1.41;sd=.04
                self.bonds.append((i,j,r,t,target,sd))
            for n1,n2 in combinations(ns,2):
                if a['label']=='C12':
                    target=np.cos(np.deg2rad(125 if seed_atoms[n1[0]]['element']=='O' and seed_atoms[n2[0]]['element']=='O' else 117.5))
                else:target=-.5
                self.angles.append((i,n1,n2,target,.045))
        # Individual phenyl and carboxylate planes, representative only.
        for labels in [['C2','C4','C10','C8','C3','C6','C12','C5'],['O4','O6','C12','C10']]:
            self.planes.append([(self.lab[l],np.eye(3),np.zeros(3)) for l in labels])
        # Symmetry-generated 16-member pyrene carbon skeleton near its center.
        refs=[];center=np.array([.5,.25,.5])
        for l in ['C1','C5','C7','C9','C11']:
            j=self.lab[l];a=seed_atoms[j]
            for r in np.array(a['reps']):
                x=r@np.array(a['xyz']);t=np.rint(center-x)
                if np.linalg.norm(A@(x+t-center))<4.1:refs.append((j,r,t))
        if len(refs)!=16:raise ValueError('Expected 16 solved pyrene carbons, got '+str(len(refs)))
        self.planes.append(refs)

    def xyz(self,p):
        xs=[]
        for s in self.model.spec:
            a=s['a'];off=s['off'];nd=s['nd']
            xs.append(np.array(a['xyz'])+np.array(a['basis'])@p[off:off+nd])
        return np.array(xs)

    def fun(self,p):
        xs=self.xyz(p);res=[]
        for i,j,r,t,target,sd in self.bonds:
            d=A@(r@xs[j]+t-xs[i]);res.append((np.linalg.norm(d)-target)/sd)
        for i,n1,n2,target,sd in self.angles:
            j,r,t,_=n1;k,s,u,_=n2
            d=A@(r@xs[j]+t-xs[i]);e=A@(s@xs[k]+u-xs[i])
            co=d@e/(np.linalg.norm(d)*np.linalg.norm(e));res.append((co-target)/sd)
        for refs in self.planes:
            xyz=np.array([A@(r@xs[j]+t) for j,r,t in refs]);xyz-=xyz.mean(axis=0)
            _,_,vh=np.linalg.svd(xyz,full_matrices=False);n=vh[-1]
            if n[np.argmax(abs(n))]<0:n=-n
            res.extend((xyz@n/.035).tolist())
        return np.array(res)*self.strength

    def jac(self,p):
        n=len(self.fun(p));J=np.zeros((n,len(p)));eps=1e-5
        for j in self.indices:
            a=p.copy();b=p.copy();a[j]+=eps;b[j]-=eps
            J[:,j]=(self.fun(a)-self.fun(b))/(2*eps)
        return J

    def adp_matrix(self):
        rows=[]
        for i,j,r,t,target,sd in self.bonds:
            si=self.model.spec[i];sj=self.model.spec[j];Q=A@r@AI
            xi=np.array(si['a']['xyz']);xj=r@np.array(sj['a']['xyz'])+t
            bond=A@(xj-xi);n=bond/np.linalg.norm(bond)
            ubi=si['ub'];ubj=np.array([Q@u@Q.T for u in sj['ub']])
            oi=si['off']+si['nd'];oj=sj['off']+sj['nd']
            # Rigid bond: mean square displacement parallel to bond similar.
            row=np.zeros(len(self.model.p0));row[oi:oi+si['nu']]+=np.einsum('i,tij,j->t',n,ubi,n)/.008;row[oj:oj+sj['nu']]-=np.einsum('i,tij,j->t',n,ubj,n)/.008;rows.append(row)
            # Similarity of directly bonded light-atom Cartesian ADPs.
            for e in E6:
                row=np.zeros(len(self.model.p0));row[oi:oi+si['nu']]+=np.einsum('ij,tij->t',e,ubi)/.06;row[oj:oj+sj['nu']]-=np.einsum('ij,tij->t',e,ubj)/.06;rows.append(row)
        return np.array(rows)
