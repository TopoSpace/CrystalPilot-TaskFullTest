"""Independent CIF read/complex summation plus analytic-derivative regression test.
Only Python stdlib/numpy; the derivative test separately imports our in-run model.
"""
from pathlib import Path
import json,re,shlex
from collections import Counter
import numpy as np

OUT=Path(__file__).resolve().parent

def read_cif(path):
    tokens=[];block=None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith(';'):
            if block is None:block=[]
            else:tokens.append('\n'.join(block));block=None
        elif block is not None:block.append(line)
        else:tokens.extend(shlex.split(line,comments=True))
    assert block is None,'Unterminated semicolon text'
    values={};loops=[];i=0
    def control(t):return t=='loop_' or t.startswith('_') or t.startswith('data_')
    while i<len(tokens):
        tok=tokens[i];i+=1
        if tok=='loop_':
            names=[]
            while i<len(tokens) and tokens[i].startswith('_'):names.append(tokens[i]);i+=1
            rows=[]
            while i<len(tokens) and not control(tokens[i]):
                rows.append(tokens[i:i+len(names)]);i+=len(names)
            assert all(len(row)==len(names) for row in rows)
            loops.append([dict(zip(names,row)) for row in rows])
        elif tok.startswith('_'):values[tok]=tokens[i];i+=1
        elif not tok.startswith('data_'):raise ValueError('Unexpected CIF token '+tok)
    return values,loops

def number(s):return float(s.split('(')[0])
def loop_with(loops,key):return next(loop for loop in loops if loop and key in loop[0])
def rotation(expr):
    rows=[]
    for term in expr.split(','):
        row=[]
        assert re.fullmatch('[xyz+\\-]+',term)
        for v in 'xyz':
            matches=re.findall('([+-]?)'+v,term);row.append(sum(-1 if sign=='-' else 1 for sign in matches))
        rows.append(row)
    return np.array(rows,int)

def main():
    values,loops=read_cif(OUT/'nu1000_framework.cif');report=json.loads((OUT/'final_model.json').read_text())
    fcfv,fcfl=read_cif(OUT/'final_reflections.fcf');refs=loop_with(fcfl,'_refln_index_h')
    h=np.array([[int(a['_refln_index_'+j]) for j in 'hkl'] for a in refs]);stored=np.array([number(a['_refln_F_squared_calc']) for a in refs])
    a,b,c=[number(values['_cell_length_'+j]) for j in 'abc'];alpha,beta,gamma=np.deg2rad([number(values['_cell_angle_'+j]) for j in ['alpha','beta','gamma']])
    cx=c*np.cos(beta);cy=c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
    cell=np.array([[a,b*np.cos(gamma),cx],[0,b*np.sin(gamma),cy],[0,0,np.sqrt(c*c-cx*cx-cy*cy)]])
    inverse=np.linalg.inv(cell);reciprocal=h@inverse;q=np.sum(reciprocal**2,axis=1);astar=np.linalg.norm(inverse,axis=1)
    rots=[rotation(row['_space_group_symop_operation_xyz']) for row in loop_with(loops,'_space_group_symop_id')]
    tensors={}
    for row in loop_with(loops,'_atom_site_aniso_label'):
        u=np.zeros((3,3))
        for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]:u[i,j]=u[j,i]=number(row[f'_atom_site_aniso_U_{i+1}{j+1}'])
        tensors[row['_atom_site_aniso_label']]=cell@(u*np.outer(astar,astar))@cell.T
    disp={row['_atom_type_symbol']:complex(number(row['_atom_type_scat_dispersion_real']),number(row['_atom_type_scat_dispersion_imag'])) for row in loop_with(loops,'_atom_type_symbol')}
    # Same published analytical approximation constants, but independent evaluator.
    coeff={'C':([2.31,1.02,1.5886,.865],[20.8439,10.2075,.5687,51.6512],.2156),
        'O':([3.0485,2.2868,1.5463,.867],[13.2771,5.7011,.3239,32.9089],.2508),
        'Zr':([17.8765,10.948,5.41732,3.65721],[1.27618,11.916,.117622,87.6627],2.06929),
        'H':([.493002,.322912,.140191,.04081],[10.5109,26.1257,3.14236,57.7997],.003038)}
    F=np.zeros(len(h),complex);inventory=Counter();multiplicity_ok=True;minu=1.;positions={}
    for row in loop_with(loops,'_atom_site_label'):
        label=row['_atom_site_label'];el=row['_atom_site_type_symbol'];occ=number(row['_atom_site_occupancy'])
        x=np.array([number(row['_atom_site_fract_'+j]) for j in 'xyz']);positions[label]=x
        U=tensors.get(label,np.eye(3)*number(row['_atom_site_U_iso_or_equiv']));minu=min(minu,float(np.linalg.eigvalsh(U).min()))
        ca,cb,cc=coeff[el];ff=np.full(len(h),cc,dtype=complex)+disp[el]
        for aa,bb in zip(ca,cb):ff+=aa*np.exp(-bb*q/4)
        orbit=[]
        for rot in rots:
            xx=(rot@x)%1
            if any(np.linalg.norm((xx-other-np.rint(xx-other))@cell.T)<1e-5 for other in orbit):continue
            orbit.append(xx);cartrot=cell@rot@inverse;uc=cartrot@U@cartrot.T
            damp=np.exp(-2*np.pi**2*np.sum((reciprocal@uc)*reciprocal,axis=1))
            F+=occ*ff*damp*np.exp(2j*np.pi*(h@xx))
        multiplicity_ok &= len(orbit)==int(row['_atom_site_symmetry_multiplicity']);inventory[el]+=len(orbit)*occ
    calc=report['scale']**2*abs(F)**2
    normalized_error=float(np.max(abs(calc-stored))/np.max(stored));relative_l2=float(np.linalg.norm(calc-stored)/np.linalg.norm(stored))
    # Check the CIF symmetry codes in every exported bond.
    max_bond_error=0.
    for row in loop_with(loops,'_geom_bond_atom_site_label_1'):
        x=positions[row['_geom_bond_atom_site_label_1']];y=positions[row['_geom_bond_atom_site_label_2']];code=row['_geom_bond_site_symmetry_2']
        if code!='.':
            ri,shift=code.split('_');y=rots[int(ri)-1]@y+np.array([int(s)-5 for s in shift])
        distance=np.linalg.norm(cell@(y-x));max_bond_error=max(max_bond_error,abs(distance-number(row['_geom_bond_distance'])))
    # Test all columns of the fitted analytical derivative, including the riding H terms.
    from anisotropic import AnisoModel
    data=np.load(OUT/'merged_equal.npz');saved=np.load(OUT/'framework_final_unmasked_fc.npz')
    model=AnisoModel(report['atoms'],data['h'],True,hydrogens=report['riding_hydrogens']);p=model.p0.copy();p[-3]=report['Zr_fp'];p[-2]=report['Zr_fpp_squared'];p[-1]=np.log(report['scale'])
    exact=model.calc(p);saved_error=float(np.max(abs(exact-saved['ic']))/max(saved['ic']))
    rng=np.random.default_rng(68199);indices=np.unique(np.r_[np.argsort(data['q'])[:40],rng.choice(len(h),160,replace=False)])
    model=AnisoModel(report['atoms'],data['h'][indices],True,hydrogens=report['riding_hydrogens']);p0=model.p0.copy();p0[-3:]=p[-3:]
    _,J=model.calc(p0,True);errors=[]
    for j in range(len(p0)):
        pa=p0.copy();pb=p0.copy();eps=1e-6;pa[j]+=eps;pb[j]-=eps;numerical=(model.calc(pa)-model.calc(pb))/(2*eps)
        errors.append(float(np.max(abs(numerical-J[:,j]))/max(1.,np.max(abs(numerical)))))
    out={'CIF_loop_parse_passed':True,'unit_cell_inventory':dict(inventory),'multiplicities_match':bool(multiplicity_ok),
        'minimum_ADP_eigenvalue_A2':minu,'CIF_vs_FCF_max_difference_normalized_to_max_Ic':normalized_error,
        'CIF_vs_FCF_relative_L2_error':relative_l2,'CIF_bond_symmetry_max_distance_error_A':float(max_bond_error),
        'saved_JSON_with_riding_H_vs_original_Ic_relative_error':saved_error,'riding_H_analytic_jacobian_max_relative_error':max(errors),
        'note':'Independent CIF reader/direct complex summation verifies encoding and implemented equations, not model completeness or experimental validity.'}
    out['all_passed']=bool(multiplicity_ok and inventory=={'Zr':18,'O':96,'C':264,'H':132} and minu>0 and normalized_error<1e-5 and relative_l2<1e-5 and max_bond_error<1e-4 and saved_error<1e-10 and max(errors)<1e-5)
    (OUT/'final_independent_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
    assert out['all_passed'],'Final verification failed; inspect report'

if __name__=='__main__':main()
