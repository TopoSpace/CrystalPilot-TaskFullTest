"""Write task-local CIF/FCF without a crystallographic library."""
import json,datetime,itertools
import numpy as np
from raw_frames import OUT
from refine import unpack

OPS=['x,y,z','x+1/2,-y+1/2,-z','-x,y+1/2,-z+1/2','-x+1/2,-y,z+1/2']

def with_su(v,s):
 if not np.isfinite(s) or s<=0:return f'{v:.7f}'
 places=max(0,int(-np.floor(np.log10(s)))+1)
 return f'{v:.{places}f}({int(round(s*10**places))})'

def main():
 m=json.loads((OUT/'refined_model.json').read_text());stats=m['statistics']
 merge=json.loads((OUT/'merging_statistics.json').read_text());diff=json.loads((OUT/'difference_statistics.json').read_text())
 cell=np.array(m['cell']);x=np.array(m['coordinates']);u=np.array(m['U_cartesian']);cov=np.load(OUT/'refinement_covariance.npy');p=np.array(m['parameters'])
 cell_su=np.full(3,np.nan)
 gp=OUT/'geometry_validation.json'
 if gp.exists():cell_su=np.array(json.loads(gp.read_text()).get('cell_standard_uncertainty',[np.nan]*3))
 vol=float(np.prod(cell));density=4*89.093/(.602214076*vol)
 r=np.loadtxt(OUT/'refinement_reflections.csv',delimiter=',',skiprows=1)
 h=r[:,:3];d=1/np.linalg.norm(h/cell,axis=1);theta=np.rad2deg(np.arcsin(.71073/(2*d)))
 hall=np.array(list(itertools.product(*[range(int(np.ceil(c/d.min()))+1) for c in cell])))
 gn=np.linalg.norm(hall/cell,axis=1);possible=(gn>0)&(gn<=1/d.min()+1e-9)
 for axis in range(3):
  others=[j for j in range(3) if j!=axis]
  possible&=~((hall[:,others]==0).all(axis=1)&(hall[:,axis]%2==1))
 actual_completeness=len(h)/int(possible.sum())
 gcv=OUT/'geometry_covariance.npy'
 vol_su=np.nan
 if gcv.exists():
  cv=np.load(gcv)[:3,:3];grad=vol/cell;vol_su=float(np.sqrt(grad@cv@grad))
 ue=np.trace(u,axis1=1,axis2=2)/3
 # Propagate U uncertainties through the Cholesky parameterization.
 def uv(p):
  _,us,_,_=unpack(p,'aniso',cell,True)
  return np.c_[np.trace(us,axis1=1,axis2=2)/3,us[:,0,0],us[:,1,1],us[:,2,2],us[:,0,1],us[:,0,2],us[:,1,2]].ravel()
 j=np.zeros((13*7,len(p)))
 for a in range(len(p)):
  pp=p.copy();pm=p.copy();pp[a]+=1e-6;pm[a]-=1e-6;j[:,a]=(uv(pp)-uv(pm))/2e-6
 us=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',j,cov,j),0)).reshape(13,7)
 lines=['data_alanine_task_local',"_audit_creation_date '2026-09-22'","_audit_creation_method 'Independent task-local Python algorithms; numpy and scipy only'",'_audit_update_record',";",'Raw-frame decoding, geometry, integration, symmetry merging, direct-space', 'solution and least-squares refinement independently implemented in this task.', 'No crystallographic software/library or external reference structure used.', 'Provisional research model: see REPORT_zh.md and validation.json.', ';',"_chemical_name_common 'alanine'","_chemical_formula_sum 'C3 H7 N O2'",'_chemical_formula_weight 89.093',"_chemical_absolute_configuration 'unk'", "_space_group_name_H-M_alt 'P 21 21 21'","_space_group_name_Hall 'P 2ac 2ab'",'_space_group_IT_number 19',"_symmetry_cell_setting 'orthorhombic'",'loop_','_space_group_symop_id','_space_group_symop_operation_xyz']
 lines += [f"{i+1} '{op}'" for i,op in enumerate(OPS)]
 for name,v,s in zip('abc',cell,cell_su):lines.append(f'_cell_length_{name} {with_su(v,s)}')
 lines += ['_cell_angle_alpha 90','_cell_angle_beta 90','_cell_angle_gamma 90',f'_cell_volume {with_su(vol,vol_su)}','_cell_formula_units_Z 4','_cell_measurement_temperature 150.0','_diffrn_ambient_temperature 150.0',"_diffrn_radiation_type 'Mo Kalpha'",'_diffrn_radiation_wavelength 0.71073',"_diffrn_measurement_device_type 'HyPix Arc 100 (supplied acquisition metadata)'","_diffrn_measurement_method 'omega scans, 0.5 degree frames, 2 s exposure'","_diffrn_detector 'hybrid pixel area detector'",'_exptl_crystal_colour ?', '_exptl_crystal_description ?', '_exptl_crystal_size_max ?', '_exptl_crystal_size_mid ?', '_exptl_crystal_size_min ?', f'_exptl_crystal_density_diffrn {density:.5f}', '_exptl_crystal_F_000 192', '_exptl_absorpt_coefficient_mu ?', "_exptl_absorpt_correction_type 'none'",'_exptl_absorpt_process_details',';','No physical absorption correction. Relative run/module scale factors', 'derived from symmetry equivalents only; blocked/invalid pixels flagged.', ';',f"_diffrn_reflns_number {merge['n_used']-merge['n_rejected']}",f"_diffrn_reflns_av_R_equivalents {merge['Rint']:.6f}",f'_diffrn_reflns_theta_min {theta.min():.5f}',f'_diffrn_reflns_theta_max {theta.max():.5f}',f'_diffrn_measured_fraction_theta_max {actual_completeness:.6f}',f'_diffrn_reflns_theta_full {np.rad2deg(np.arcsin(.71073/1.6)):.5f}',f"_diffrn_measured_fraction_theta_full {merge['completeness_d0.8']['fraction']:.6f}",f"_reflns_number_total {stats['n_reflections']}",f"_reflns_number_gt {stats['n_gt2sigma']}","_reflns_threshold_expression 'I > 2 sigma(I)'","_refine_ls_structure_factor_coef 'Fsqd'","_refine_ls_matrix_type 'full'",f"_refine_ls_number_reflns {stats['n_reflections']}",f"_refine_ls_number_parameters {stats['n_parameters']}",'_refine_ls_number_restraints 0',f"_refine_ls_R_factor_gt {stats['R1_gt2sigma']:.6f}",f"_refine_ls_R_factor_all {stats['R1_all']:.6f}",f"_refine_ls_wR_factor_ref {stats['wR2']:.6f}",f"_refine_ls_goodness_of_fit_ref {stats['goodness_of_fit']:.6f}","_refine_ls_weighting_scheme 'calc'",f"_refine_ls_weighting_details '{m['weight_model']}'","_refine_ls_hydrogen_treatment 'constr'",'_refine_ls_abs_structure_Flack ?', '_refine_ls_abs_structure_details',';','Absolute configuration undetermined. Friedel equivalents merged.', 'Template hand chosen arbitrarily; no anomalous-scattering refinement.', ';',f"_refine_diff_density_max {diff['max_e_A3']:.6f}",f"_refine_diff_density_min {diff['min_e_A3']:.6f}",f"_refine_diff_density_rms {diff['rms_e_A3']:.6f}",'_refine_special_details',';','All non-H coordinates and anisotropic U freely refined on F squared.', 'H atoms ride; CA-H 0.98 A, NH3-H 0.91 A, CH3-H 0.98 A.', 'NH3 and CH3 torsions refined. H Uiso multipliers 1.2 (CA) and 1.5.', 'Neutral-atom four-Gaussian form factors; no dispersion terms.', 'Covariance scaled by residual variance; uncertainties remain conditional', 'on the intensity-error, geometry and riding-H models.', 'No independent standard crystallographic validation was available.', ';','loop_','_atom_site_label','_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z','_atom_site_U_iso_or_equiv','_atom_site_adp_type','_atom_site_occupancy','_atom_site_calc_flag']
 for i,(lab,typ) in enumerate(zip(m['labels'],m['types'])):
  cs=m['coordinate_su_nonH'][i] if i<6 else [np.nan]*3
  coords=' '.join(with_su(v,s) for v,s in zip(x[i],cs))
  lines.append(f"{lab} {typ} {coords} {with_su(ue[i],us[i,0]) if i<6 else f'{ue[i]:.6f}'} {'Uani' if i<6 else 'Uiso'} 1 {'d' if i<6 else 'calc'}")
 lines += ['loop_','_atom_site_aniso_label','_atom_site_aniso_U_11','_atom_site_aniso_U_22','_atom_site_aniso_U_33','_atom_site_aniso_U_12','_atom_site_aniso_U_13','_atom_site_aniso_U_23']
 for i in range(6):
  vals=[u[i,0,0],u[i,1,1],u[i,2,2],u[i,0,1],u[i,0,2],u[i,1,2]]
  lines.append(m['labels'][i]+' '+' '.join(with_su(v,s) for v,s in zip(vals,us[i,1:])))
 bonds=json.loads((OUT/'bond_geometry.json').read_text())
 lines += ['loop_','_geom_bond_atom_site_label_1','_geom_bond_atom_site_label_2','_geom_bond_distance','_geom_bond_site_symmetry_2']
 lines += [f"{b['atom1']} {b['atom2']} {with_su(b['distance'],b['su'])} ." for b in bonds]
 (OUT/'alanine.cif').write_text('\n'.join(lines)+'\n',encoding='utf8')
 fcf=['data_alanine_reflections',"_audit_creation_method 'Independent Python structure-factor calculation'",'_reflns_special_details',';','Fo2 and its standard uncertainty are divided by the refined intensity', 'scale to match Fc2. All allowed merged reflections included, negatives', 'retained. Phases in degrees; no unobserved intensities are invented.', ';','loop_','_refln_index_h','_refln_index_k','_refln_index_l','_refln_F_squared_meas','_refln_F_squared_sigma','_refln_F_squared_calc','_refln_phase_calc']
 for row in r:
  a,b,c,I,s,F,fr,fi=row[:8]
  fcf.append(f'{int(a)} {int(b)} {int(c)} {I/m["scale"]:.7f} {s/m["scale"]:.7f} {F/m["scale"]:.7f} {np.rad2deg(np.arctan2(fi,fr)):.6f}')
 (OUT/'alanine.fcf').write_text('\n'.join(fcf)+'\n',encoding='utf8')
 cart=x*cell
 (OUT/'alanine_molecule.xyz').write_text('13\nOne unwrapped asymmetric-unit molecule; Cartesian A; arbitrary hand\n'+'\n'.join(f'{t} {v[0]:.7f} {v[1]:.7f} {v[2]:.7f}' for t,v in zip(m['types'],cart))+'\n')
 (OUT/'export_statistics.json').write_text(json.dumps(dict(volume_A3=vol,volume_su_A3=vol_su,density_g_cm3=density,F000=192,theta_min=float(theta.min()),theta_max=float(theta.max()),actual_d_min=float(d.min()),completeness_at_actual_d_min=actual_completeness,n_possible_actual_dmin=int(possible.sum()),cif='alanine.cif',fcf='alanine.fcf'),indent=2))
 print('Exported alanine.cif, alanine.fcf and alanine_molecule.xyz',flush=True)
if __name__=='__main__':main()
