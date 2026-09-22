from anisotropic import *
import json

out={}
# Verify group, metric invariance and inverse membership.
out['symmetry_order']=len(R)
out['group_closed']=all(any(np.array_equal(a@b,c) for c in R) for a in R for b in R)
out['metric_error']=float(max(np.max(abs(r.T@G@r-G)) for r in R))
report=json.loads((OUT/'framework_iso1.json').read_text());d=np.load(OUT/'merged_equal.npz')
rng=np.random.default_rng(23432);h=d['h'][rng.choice(len(d['h']),100,replace=False)]
m=AnisoModel(report['atoms'],h,True);p=m.p0;p[-1]=np.log(report['scale']);p[-3]=-3;p[-2]=9
ic,jac=m.calc(p,True);errors=[]
for j in range(len(p)):
    dp=1e-6;pa=p.copy();pb=p.copy();pa[j]+=dp;pb[j]-=dp
    numeric=(m.calc(pa)-m.calc(pb))/(2*dp)
    error=np.max(abs(numeric-jac[:,j]))/max(1,np.max(abs(numeric)))
    errors.append(error)
out['analytic_jacobian_max_relative_error']=float(max(errors))
# Explicit complex summation over all unit-cell atoms; no reduced real expression.
F=np.zeros(len(h),complex);q=np.einsum('ni,ij,nj->n',h,GI,h)
for a in report['atoms']:
    xyz=np.einsum('sij,j->si',a['reps'],a['xyz'])%1
    ff=formfactor(a['element'],q/4).astype(complex)
    if a['element']=='Zr':ff+=-3+3j
    F+=ff*np.exp(-2*np.pi**2*q*a['Uiso'])*a['occ']*np.exp(2j*np.pi*(h@xyz.T)).sum(axis=1)
iref=report['scale']**2*abs(F)**2
out['direct_complex_F2_relative_error']=float(np.max(abs(iref-ic))/np.max(iref))
# Intensity invariance across group rotations for non-isotropic test tensors.
for a in report['atoms']:
    ub=tensor_basis(a['stabilizer']);uc=np.eye(3)*a['Uiso']+.003*ub[0];a['Ucart']=uc.tolist()
m=AnisoModel(report['atoms'],h,True);p=m.p0;p[-1]=-3;p[-3]=-3;p[-2]=9
ref=m.calc(p);errs=[]
for r in R:
    mm=AnisoModel(report['atoms'],h@r,True);pp=mm.p0;pp[-3:]=p[-3:]
    errs.append(np.max(abs(mm.calc(pp)-ref))/max(ref))
out['anisotropic_symmetry_invariance_relative_error']=float(max(errs))
out['all_passed']=bool(out['group_closed'] and out['metric_error']<1e-8 and out['analytic_jacobian_max_relative_error']<1e-5 and out['direct_complex_F2_relative_error']<1e-9 and out['anisotropic_symmetry_invariance_relative_error']<1e-9)
print(json.dumps(out,indent=2));(OUT/'algorithm_validation.json').write_text(json.dumps(out,indent=2))
assert out['all_passed']
