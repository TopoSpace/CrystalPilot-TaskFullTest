"""Plain-text CIF/XYZ export, implemented without crystallographic libraries."""
from crystal import *
from collections import Counter
import argparse

def op_string(r):
    out=[]
    for row in r:
        s=''
        for c,v in zip(row,'xyz'):
            if c:s+=('+' if c>0 and s else '')+('-' if c<0 else '')+v
        out.append(s or '0')
    return ','.join(out)

def su_number(value,su=None,minimum_decimals=0):
    if su is None or su<=1e-12:return f'{value:.8f}'
    nd=max(minimum_decimals,min(8,int(-np.floor(np.log10(su)))+1))
    return f'{value:.{nd}f}({int(round(su*10**nd))})'

def export(report,path,description=None):
    atoms=report['atoms']+report.get('riding_hydrogens',[]);counts=Counter()
    for a in atoms:counts[a['element']]+=len(a['reps'])*a['occ']
    lines=['data_NU1000_framework','',"_audit_creation_method 'Independent Python/numpy/scipy implementation'",'_audit_creation_date 2026-09-22',
        '_audit_block_doi ?', '_chemical_name_common NU-1000',
        '_chemical_formula_sum ?', '_chemical_formula_weight ?', '_cell_formula_units_Z ?',
        '_cell_length_a 39.1900','_cell_length_b 39.1900','_cell_length_c 16.6100','_cell_angle_alpha 90','_cell_angle_beta 90','_cell_angle_gamma 120',f'_cell_volume {V:.4f}',
        "_space_group_name_H-M_alt 'P 6/m m m'", "_space_group_name_Hall '-P 6 2'", '_space_group_IT_number 191',
        '_diffrn_radiation_wavelength 0.68883','_diffrn_radiation_type x-ray','_diffrn_ambient_temperature ?', '_diffrn_radiation_source ?',
        '_exptl_crystal_description ?', '_exptl_crystal_colour ?', '_exptl_crystal_size_max ?', '_exptl_crystal_size_mid ?', '_exptl_crystal_size_min ?',
        '_refine_special_details',';',(description or 'Framework-only checkpoint. No solvent/guest atoms. Hydrogen positions not assigned. Do not treat as publication-validated.'),
        'The space group is supported by the input symmetry and strong-reflection equivalence.','The cell and wavelength are copied from start.ins; cell uncertainties and collection temperature are unknown.',
        'Unit-cell atom inventory: '+', '.join(f'{el} {n:g}' for el,n in counts.items()),
        'The formula, formula weight and Z are left unknown because node protonation and pore content are not determined.',';',
        'loop_','_space_group_symop_id','_space_group_symop_operation_xyz']
    rots=sorted(R.tolist(),key=lambda r:(not np.array_equal(r,np.eye(3)),op_string(r)))
    lines.extend(f"{i+1} '{op_string(r)}'" for i,r in enumerate(rots))
    if 'final_statistics' in report:
        st=report['final_statistics'];fit=st['refinement']
        lines += ['', '_refine_ls_structure_factor_coef Fsqd', '_refine_ls_matrix_type full', '_refine_ls_hydrogen_treatment constr',
            f"_refine_ls_number_reflns {fit['n']}",f"_refine_ls_number_parameters {report['n_params']}",f"_refine_ls_number_restraints {report.get('n_effective_restraint_rows',report.get('n_restraints','?'))}",
            f"_refine_ls_R_factor_gt {fit['R1_obs']:.6f}",f"_refine_ls_R_factor_all {fit['R1_all']:.6f}",f"_refine_ls_wR_factor_ref {fit['wR2_actual_weights']:.6f}",f"_refine_ls_goodness_of_fit_ref {fit['GoF']:.5f}",
            "_reflns_threshold_expression 'I > 2 sigma(I)'", f"_reflns_number_total {fit['n']}",f"_reflns_number_gt {fit['n_obs']}", f"_reflns_d_resolution_high {report['dmin']:.4f}",f"_reflns_d_resolution_low {st['measured_d_max']:.5f}",
            '_refine_ls_weighting_scheme calc',"_refine_ls_weighting_details 'w=1/[sigma(I)^2+(0.08P)^2+0.01]; P=(max(I,0)+2Ic)/3'",
            f"_refine_diff_density_max {st['difference_map']['maximum']:.5f}",f"_refine_diff_density_min {st['difference_map']['minimum']:.5f}",f"_refine_diff_density_rms {st['difference_map']['rms']:.5f}",
            'loop_','_atom_type_symbol','_atom_type_scat_dispersion_real','_atom_type_scat_dispersion_imag',
            'C 0 0','H 0 0','O 0 0',f"Zr {report['Zr_fp']:.6f} {np.sqrt(report['Zr_fpp_squared']):.6f}"]
    lines += ['','loop_','_atom_site_label','_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z','_atom_site_U_iso_or_equiv','_atom_site_adp_type','_atom_site_occupancy','_atom_site_symmetry_multiplicity','_atom_site_calc_flag','_atom_site_refinement_flags_posn']
    for a in atoms:
        x=np.array(a['xyz'])%1;adp='Uani' if 'Ucart' in a else 'Uiso';s=a.get('xyz_su',[None]*3)
        coord=' '.join(su_number(v,u,8) for v,u in zip(x,s));uiso=su_number(a['Uiso'],a.get('Uiso_su'),8)
        flags='calc R' if a['element']=='H' else 'd .'
        lines.append(f"{a['label']} {a['element']} {coord} {uiso} {adp} {a['occ']:.6f} {len(a['reps'])} {flags}")
    ani=[a for a in atoms if 'Ucart' in a]
    if ani:
        lines+=['','loop_','_atom_site_aniso_label','_atom_site_aniso_U_11','_atom_site_aniso_U_22','_atom_site_aniso_U_33','_atom_site_aniso_U_12','_atom_site_aniso_U_13','_atom_site_aniso_U_23']
        astar=np.sqrt(np.diag(GI));AI=np.linalg.inv(A)
        for a in ani:
            us=AI@np.array(a['Ucart'])@AI.T;u=us/np.outer(astar,astar)
            su=np.array(a.get('Uij_su',np.zeros((3,3))))
            pairs=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
            lines.append(a['label']+' '+' '.join(su_number(u[i,j],su[i,j],8) for i,j in pairs))
    if report.get('geometry_bonds'):
        lines += ['', 'loop_','_geom_bond_atom_site_label_1','_geom_bond_atom_site_label_2','_geom_bond_distance','_geom_bond_site_symmetry_2']
        for b in report['geometry_bonds']:lines.append(f"{b['a']} {b['b']} {su_number(b['distance'],b['su'],4)} {b['symmetry_2']}")
    Path(path).write_text('\n'.join(lines)+'\n',encoding='utf-8')
    xyz=[]
    for a in atoms:
        for x in expand_xyz(a['xyz']):xyz.append((a['element'],A@x,a['label']))
    Path(path).with_suffix('.xyz').write_text(str(len(xyz))+'\n'+'Framework unit cell, Cartesian Angstrom. Cell columns: '+str(A.T.tolist())+'\n'+'\n'.join(f'{e} {x[0]:.8f} {x[1]:.8f} {x[2]:.8f}' for e,x,l in xyz)+'\n')
    print('Exported',path,'expanded atom count',len(xyz),'inventory',dict(counts))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model');p.add_argument('output');a=p.parse_args()
    report=json.loads((OUT/a.model).read_text());export(report,OUT/a.output)
