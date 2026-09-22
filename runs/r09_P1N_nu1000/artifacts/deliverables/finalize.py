"""Create final checkable structure, formal uncertainties and exact refinement statistics."""
from anisotropic import *
from export_structure import export,op_string
from structure_checks import check
from restraints import Geometry
import hashlib,csv,datetime

TAG='framework_final_unmasked'

def actual_stats(I,sig,ic,mask,npar=0,a=.08):
    ii=I[mask];ss=sig[mask];cc=ic[mask];fo=np.sqrt(np.maximum(ii,0));fc=np.sqrt(np.maximum(cc,0));obs=ii>2*ss
    P=(np.maximum(ii,0)+2*cc)/3;w=1/(ss**2+(a*P)**2+.01);chi=float(np.sum(w*(ii-cc)**2))
    return {'n':len(ii),'n_obs':int(obs.sum()),'R1_all':float(np.sum(abs(fo-fc))/np.sum(fo)),
        'R1_obs':float(np.sum(abs(fo[obs]-fc[obs]))/np.sum(fo[obs])),
        'wR2_actual_weights':float(np.sqrt(chi/np.sum(w*ii**2))),'weighted_chi_square':chi,
        'GoF':float(np.sqrt(chi/max(len(ii)-npar,1))),'n_parameters_subtracted':npar}

def main():
    r=json.loads((OUT/(TAG+'.json')).read_text());d=np.load(OUT/'merged_equal.npz');f=np.load(OUT/(TAG+'_fc.npz'));c=np.load(OUT/(TAG+'_covariance.npz'))
    h=d['h'];I=d['I'];sig=d['sig'];q=d['q'];fit=f['fit'];ic=f['ic'];real=f['real'];zr=f['zr'];scale=float(f['scale']);fp=float(f['fp']);t=float(f['t'])
    st=actual_stats(I,sig,ic,fit,r['n_params']);cov=c['cov']*st['GoF']**2
    m=AnisoModel(r['atoms'],h[:1],True,hydrogens=r['riding_hydrogens'])
    # Formal uncertainties conditional on the chosen scattering/density/geometry model.
    for a,s in zip(r['atoms'],m.spec):
        off=s['off'];nd=s['nd'];nu=s['nu'];basis=np.array(a['basis']);cp=cov[off:off+nd,off:off+nd]
        cf=basis@cp@basis.T;a['xyz_su']=np.sqrt(np.maximum(np.diag(cf),0)).tolist()
        cu=cov[off+nd:off+nd+nu,off+nd:off+nd+nu];ub=s['ub'];tr=np.trace(ub,axis1=1,axis2=2)/3
        a['Uiso_su']=float(np.sqrt(max(tr@cu@tr,0)))
        astar=np.sqrt(np.diag(GI));bu=np.array([AI@u@AI.T/np.outer(astar,astar) for u in ub])
        a['Uij_su']=np.sqrt(np.maximum(np.einsum('tij,tu,uij->ij',bu,cu,bu),0)).tolist()
    r['effective_dispersion_uncertainties']={'fp_su':float(np.sqrt(cov[-3,-3])),'fpp':float(np.sqrt(t)),'fpp_su':float(np.sqrt(cov[-2,-2])/(2*np.sqrt(t)))}
    # Atom displacement between independent resolution treatments quantifies model sensitivity.
    hi=json.loads((OUT/'framework_h_highangle.json').read_text());comparisons=[]
    for a,b in zip(r['atoms'],hi['atoms']):comparisons.append({'atom':a['label'],'delta_A_all_vs_highangle':float(np.linalg.norm(distvec(a['xyz'],b['xyz'])))})
    # Difference-map statistics with their own explicit resolution/projection convention.
    rho=np.load(OUT/(TAG+'_diff_map.npy'));diff={'maximum':float(rho.max()),'minimum':float(rho.min()),'rms':float(np.std(rho)),
        'max_fractional':(np.array(np.unravel_index(np.argmax(rho),rho.shape))/rho.shape).tolist(),
        'min_fractional':(np.array(np.unravel_index(np.argmin(rho),rho.shape))/rho.shape).tolist(),
        'coefficient':'Re[(Fo/|Fc|-1)Fc], scaled to electron units; F000=0; 1.0 A resolution; negative I gives Fo=0'}
    rawh,rawI,rawsig=read_raw();hc=canonical(rawh);qh=np.einsum('ni,ij,nj->n',rawh,GI,rawh)
    center={tuple(x):y for x,y in zip(h,I)};mean=np.array([center[tuple(x)] for x in hc]);rm=qh<=1.
    rint=float(np.sum(abs(rawI[rm]-mean[rm]))/np.sum(abs(rawI[rm])))
    cv=json.loads((OUT/'framework_conservative_cv.json').read_text())
    cvf=np.load(OUT/'framework_conservative_cv_fc.npz');cvm=cvf['work']&(q<=1.);cvt=(~cvf['work'])&(q<=1.)
    stats={'refinement':st,'all_recorded':actual_stats(I,sig,ic,np.ones(len(I),bool)),
        'mid_high_angle_3_to_1_A_with_final_model':actual_stats(I,sig,ic,(q<=1.)&(q>=1/9)),
        'cv_work':actual_stats(I,sig,cvf['ic'],cvm,cv['n_params']),'cv_test':actual_stats(I,sig,cvf['ic'],cvt),
        'cv_caveat':'Reflections excluded from refinement; initial heavy-atom/Fourier solution used all data. Diagnostic holdout, not fully blind Rfree.',
        'raw_observations':len(rawI),'merged_unique_all':len(I),'measured_d_min':float(1/np.sqrt(q.max())),'measured_d_max':float(1/np.sqrt(q.min())),
        'Rint_equal_mean_to_1_A':rint,'completeness_to_1_A':4395/4401,'possible_unique_to_1_A':4401,
        'difference_map':diff,'model_sensitivity':comparisons,'weighting':'w=1/[sigma(I)^2+(0.08P)^2+0.01], P=(max(I,0)+2Ic)/3',
        'uncertainty_caveat':'Formal covariance = inverse(JtJ) times GoF^2, with restraints and fixed cell. Does NOT include model bias, near-edge dispersion calibration, cell errors or pore scattering.'}
    r['final_statistics']=stats;r['solvent_correction']='none';r['status']='framework solved; not publication validated; residual pore scattering and effective dispersion unresolved'
    # Geometry bonds with formal covariance and CIF symmetry notation.
    seed=json.loads((OUT/'framework_iso1.json').read_text())['atoms'];geo=Geometry(m,seed);bonds=[]
    rots=sorted(R.tolist(),key=lambda rr:(not np.array_equal(rr,np.eye(3)),op_string(rr)))
    for i,j,rot,trans,target,sd in geo.bonds:
        ai=r['atoms'][i];aj=r['atoms'][j];si=m.spec[i];sj=m.spec[j]
        dv=A@(rot@np.array(aj['xyz'])+trans-np.array(ai['xyz']));dist=np.linalg.norm(dv);g=np.zeros(r['n_params'])
        oi=si['off'];oj=sj['off'];di=si['nd'];dj=sj['nd'];u=dv/dist
        g[oi:oi+di]-=u@A@np.array(ai['basis']);g[oj:oj+dj]+=u@A@rot@np.array(aj['basis'])
        su=float(np.sqrt(max(g@cov@g,0)));ri=next(k+1 for k,rr in enumerate(rots) if np.array_equal(rr,rot));tr=np.rint(trans).astype(int)
        code='.' if ri==1 and np.all(tr==0) else str(ri)+'_'+''.join(str(int(v)+5) for v in tr)
        bonds.append({'a':ai['label'],'b':aj['label'],'distance':float(dist),'su':su,'symmetry_2':code,'target':target,'restraint_sigma':sd/r['geometry_strength']})
    r['geometry_bonds']=bonds
    r['n_effective_restraint_rows']=len(geo.fun(m.p0))+len(geo.adp_matrix())
    r['restraint_count_details']={'geometry_rows':len(geo.fun(m.p0)),'ADP_similarity_rows':len(geo.adp_matrix()),'inactive_ADP_eigenvalue_lower_bound_rows':60,'total_optimizer_rows':r['n_restraints'],'note':'Rows are not all statistically independent.'}
    ev=np.linalg.eigvalsh(c['jac'].T@c['jac']);stats['normal_matrix_rank_at_1e_10_cutoff']=int(np.sum(ev>ev[-1]*1e-10));stats['normal_matrix_condition_number']=float(ev[-1]/ev[0])
    xs=[];labels=[]
    for a in r['atoms']:
        ex=expand_xyz(a['xyz']);xs.extend(ex);labels.extend([a['label']]*len(ex))
    for kind in ['max','min']:
        dd=mindist(diff[kind+'_fractional'],xs);j=int(np.argmin(dd));diff[kind+'_nearest_framework_atom']={'label':labels[j],'distance_A':float(dd[j])}
    (OUT/'final_model.json').write_text(json.dumps(r,indent=2));(OUT/'final_statistics.json').write_text(json.dumps(stats,indent=2));np.savez_compressed(OUT/'final_covariance.npz',cov=cov,unscaled_cov=c['cov'],parameters=c['parameters'])
    description='''Framework-only constrained/ restrained refinement in P 6/m m m.
18 Zr, 264 C and 96 O atoms per unit cell, plus 132 idealized aromatic H atoms.
No solvent or guest atomic sites, and NO pore-density correction, are in this model.
The sixfold formula-unit inventory is not a complete chemical composition:
node O/OH/H2O protonation, terminal-site chemistry and pore contents are unknown.
Aromatic C-H = 0.95 A riding, Uiso(H)=1.2Ueq(C); node hydrogen atoms are omitted.
Twenty non-H asymmetric sites were refined anisotropically, with site symmetry.
Generic aromatic/carboxylate distances, angles, planarity, and rigid-bond/similar
ADP restraints stabilize the weak light-atom data. There are 188 geometry/ADP
residual rows, not all independent, plus 60 inactive ADP-positivity penalty rows.
Restraint targets and strengths are in the accompanying source and JSON, not
from a reference structure. Coordinates have extra decimals for reproducibility.
Resolution limited to 1.0 A using signal diagnostics; ALL measured intensities
inside that limit, including negative values, enter the final F-squared refinement.
No reflection-specific outlier removal or low-angle omission is applied.
Zr real and imaginary dispersion terms are EFFECTIVE parameters fitted to these
data, not independently calibrated anomalous-scattering values. C,H,O corrections
were approximated as zero. The supplied wavelength is 0.68883 A.
Formal standard uncertainties are conditional on this model and the fixed cell;
they do not cover unmodeled pore scattering or the uncertain dispersion correction.
This structure is a checkable framework result, NOT a publication-validated model.
Low-angle residuals remain significant. No external crystallographic validation
program, crystallographic library, structure database, or downloaded data was used.'''
    export(r,OUT/'nu1000_framework.cif',description)
    validation=check(r,'final');(OUT/'final_validation.json').write_text(json.dumps(validation,indent=2))
    # Explicit structure-factor record: every original unique reflection retained.
    phase=np.degrees(np.arctan2(np.sqrt(t)*zr,real))
    with (OUT/'final_reflections.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['h','k','l','Iobs','sigma_I','Icalc','Fc_real_e','Fc_imag_e','d_A','in_final_refinement','cv_work_flag'])
        for j in range(len(h)):w.writerow([*h[j],f'{I[j]:.10g}',f'{sig[j]:.10g}',f'{ic[j]:.10g}',f'{real[j]:.10g}',f'{np.sqrt(t)*zr[j]:.10g}',f'{1/np.sqrt(q[j]):.8f}',int(fit[j]),int(cvf['work'][j])])
    lines=['data_NU1000_reflections',"_audit_creation_method 'Independent numpy/scipy implementation'",'_refln_special_details',';',
        'All merged unique reflections, including those outside the chosen refinement limit.',
        'F_squared_calc is scaled into the SAME arbitrary intensity units as F_squared_meas.',
        'Status o: included and I>2sigma; <: included but not above threshold; x: outside resolution cutoff.',
        'Calculated complex phase includes the empirically fitted Zr imaginary dispersion term.',';',
        'loop_','_refln_index_h','_refln_index_k','_refln_index_l','_refln_F_squared_calc','_refln_F_squared_meas','_refln_F_squared_sigma','_refln_phase_calc','_refln_observed_status']
    for j,hh in enumerate(h):
        status=('o' if I[j]>2*sig[j] else '<') if fit[j] else 'x'
        lines.append(f'{hh[0]} {hh[1]} {hh[2]} {ic[j]:.8f} {I[j]:.8f} {sig[j]:.8f} {phase[j]:.5f} {status}')
    (OUT/'final_reflections.fcf').write_text('\n'.join(lines)+'\n')
    with (OUT/'asymmetric_atoms.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['label','element','x','y','z','sx_formal','sy_formal','sz_formal','Ueq','sUeq_formal','occupancy','multiplicity','treatment'])
        for a in r['atoms']+r['riding_hydrogens']:w.writerow([a['label'],a['element'],*a['xyz'],*a.get('xyz_su',['']*3),a['Uiso'],a.get('Uiso_su',''),a['occ'],len(a['reps']),'riding' if a['element']=='H' else 'anisotropic restrained'])
    manifest={'generated_at':datetime.datetime.now().isoformat(),'input_hashes':{},'python':sys.executable,'numpy_version':np.__version__,'methods':'All crystallographic algorithms written in this run; numpy/scipy only; no network/install/external crystallographic code.'}
    for name in ['start.hkl','start.ins','synthesis.png']:
        raw=(ROOT/'inputs'/name).read_bytes();manifest['input_hashes'][name]={'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    (OUT/'provenance.json').write_text(json.dumps(manifest,indent=2))
    print('FINAL STATISTICS',json.dumps(stats,indent=2),flush=True)
    print('EFFECTIVE DISPERSION',r['Zr_fp'],np.sqrt(r['Zr_fpp_squared']),r['effective_dispersion_uncertainties'],flush=True)

if __name__=='__main__':main()
