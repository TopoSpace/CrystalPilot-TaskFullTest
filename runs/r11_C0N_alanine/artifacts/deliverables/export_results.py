"""Export transparent CIF/FCF/XYZ records and compact diagnostic figures."""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from refine_aniso import LABELS,TYPES,H_LABELS,H_PARENTS,H_FACTORS,u_matrices
from solve import SIGNS,TRANS
from geometry import WAVELENGTH

OUT=Path(__file__).resolve().parent


def su(value,sigma,figures=2):
    if not np.isfinite(sigma) or sigma<=0:
        return f"{value:.7f}"
    digits=max(0,figures-1-int(math.floor(math.log10(sigma))))
    return f"{value:.{digits}f}({int(round(sigma*10**digits))})"


def loop(lines,names,rows):
    lines.append("loop_")
    lines.extend(names)
    lines.extend(" ".join(map(str,row)) for row in rows)


def write_cif(model,validation,absorption):
    z=np.array(model["parameter_vector"])
    esd=np.array(model["parameter_esd"])
    cell=np.array(model["cell"])
    cell_esd=np.array(validation["geometry_stability"]["adopted_internal_cell_esd"])
    volume=np.prod(cell)
    vs=volume*np.linalg.norm(cell_esd/cell)
    xyz=np.array(model["xyz"])
    hx=np.array(model["hydrogen_xyz"])
    u=np.array(model["uij"])
    cov=np.load(OUT/"final_parameter_covariance.npy")
    ueq=u[:,:3].mean(axis=1)
    us=[]
    for i in range(6):
        ids=np.arange(18+6*i,21+6*i)
        us.append(np.sqrt(np.sum(cov[np.ix_(ids,ids)])/9))
    stats=model["statistics"]
    refs=pd.read_csv(OUT/"final_refinement_reflections.csv")
    obs=pd.read_csv(OUT/"merged_corrected_observations.csv")
    obs=obs[obs.d>=.75]
    ri=np.abs(obs.I_scaled-obs.merge_mean).sum()/np.maximum(obs.I_scaled,0).sum()
    completeness=next(r for r in validation["completeness"] if r["d_min"]==.75)
    theta_min=np.degrees(np.arcsin(WAVELENGTH/(2*refs.d.max())))
    theta_max=np.degrees(np.arcsin(WAVELENGTH/(2*refs.d.min())))
    fulltheta=np.degrees(np.arcsin(WAVELENGTH/(2*.8)))
    cent=np.loadtxt(OUT/"geometry_flexible_residuals.csv",delimiter=",")
    norm=np.linalg.norm(cent[:,10:13]/cell,axis=1)
    celltheta=np.degrees(np.arcsin(WAVELENGTH*norm/2))
    l=["data_alanine",
       "_audit_creation_method 'Self-written Python 3.12 / NumPy / SciPy algorithms'",
       "_chemical_name_common 'alanine; absolute configuration undetermined'",
       "_chemical_formula_moiety 'C3 H7 N O2'",
       "_chemical_formula_sum 'C3 H7 N O2'",
       "_chemical_formula_weight 89.093",
       f"_cell_length_a {su(cell[0],cell_esd[0])}",
       f"_cell_length_b {su(cell[1],cell_esd[1])}",
       f"_cell_length_c {su(cell[2],cell_esd[2])}",
       "_cell_angle_alpha 90","_cell_angle_beta 90","_cell_angle_gamma 90",
       f"_cell_volume {su(volume,vs)}","_cell_formula_units_Z 4",
       f"_cell_measurement_reflns_used {len(cent)}",
       f"_cell_measurement_theta_min {celltheta.min():.3f}",
       f"_cell_measurement_theta_max {celltheta.max():.3f}",
       "_cell_measurement_temperature 150",
       "_space_group_crystal_system orthorhombic",
       "_space_group_name_H-M_alt 'P 21 21 21'","_space_group_IT_number 19",
       "_exptl_crystal_description ?","_exptl_crystal_colour ?",
       f"_exptl_crystal_density_diffrn {4*89.093/(.602214076*volume):.5f}",
       "_exptl_crystal_F_000 192","_exptl_absorpt_coefficient_mu 0.11620",
       "_exptl_absorpt_correction_type numerical",
       f"_exptl_absorpt_correction_T_min {absorption['Tmin']:.6f}",
       f"_exptl_absorpt_correction_T_max {absorption['Tmax']:.6f}",
       "_exptl_absorpt_process_details",
       ";",
       "Self-written numerical ray integration through nine supplied crystal faces.",
       "mu=0.11620 mm^-1 is supplied metadata, not independently recalculated.",
       "Uniform illuminated volume assumed. The listed transmission extrema refer",
       "only to this numerical correction. Additional scan-dependent empirical",
       "quadratic scale functions were fitted to equivalent intensities, without Fc.",
       "Scale-function complexity was selected using held-out reflection groups.",
       ";",
       "_diffrn_ambient_temperature 150",
       "_diffrn_radiation_type 'Mo Kalpha'","_diffrn_radiation_wavelength 0.71073",
       "_diffrn_measurement_device_type 'KM4 Xcalibur2 / HyPix-Arc 100 (input metadata)'",
       "_diffrn_measurement_method 'omega rotation; 0.5 degree frames; 2 s exposures'",
       f"_diffrn_reflns_number {len(obs)}",
       f"_diffrn_reflns_theta_min {theta_min:.3f}",
       f"_diffrn_reflns_theta_max {theta_max:.3f}",
       f"_diffrn_reflns_theta_full {fulltheta:.3f}",
       f"_diffrn_measured_fraction_theta_max {completeness['fraction']:.5f}",
       "_diffrn_measured_fraction_theta_full 1.00000",
       f"_diffrn_reflns_av_R_equivalents {ri:.6f}",
       f"_diffrn_reflns_limit_h_min {int(obs.h.min())}",
       f"_diffrn_reflns_limit_h_max {int(obs.h.max())}",
       f"_diffrn_reflns_limit_k_min {int(obs.k.min())}",
       f"_diffrn_reflns_limit_k_max {int(obs.k.max())}",
       f"_diffrn_reflns_limit_l_min {int(obs.l.min())}",
       f"_diffrn_reflns_limit_l_max {int(obs.l.max())}",
       f"_reflns_number_total {len(refs)}",
       f"_reflns_number_gt {stats['n_observed']}",
       "_reflns_threshold_expression 'I > 2 sigma(I)'",
       "_reflns_special_details",";",
       "Friedel mates were merged under orthorhombic Laue symmetry.",
       "All signed intensities within d >= 0.75 A were retained for refinement,",
       "apart from screw-axis systematic absences. No Fo/Fc outlier omission.",
       "Raw shadow masks were derived from diffuse background, not calculated F.",
       "Robust equivalent-reflection weights are retained in the observation CSV.",
       "The R-equivalent reported above includes all non-shadowed observations;",
       "the smaller robust-filtered diagnostic is reported separately in the summary.",
       "Measurements beyond 0.75 A are retained in merged_corrected.csv.",
       ";",
       "_refine_ls_structure_factor_coef Fsqd","_refine_ls_matrix_type full",
       f"_refine_ls_number_reflns {len(refs)}","_refine_ls_number_parameters 57",
       "_refine_ls_number_restraints 0",
       f"_refine_ls_R_factor_gt {stats['R1_gt2sigma']:.6f}",
       f"_refine_ls_R_factor_all {stats['R1_all']:.6f}",
       f"_refine_ls_wR_factor_ref {stats['wR2']:.6f}",
       f"_refine_ls_goodness_of_fit_ref {stats['gof']:.6f}",
       "_refine_ls_shift/su_max ?",
       "_refine_ls_weighting_scheme calc",
       "_refine_ls_weighting_details 'w=1/[sigma(I)^2+(0.03P)^2]; P=(max(I,0)+2Fc^2)/3'",
       "_refine_ls_hydrogen_treatment constr","_refine_ls_extinction_method none",
       "_refine_ls_abs_structure_Flack ?",
       "_refine_ls_abs_structure_details",
       ";",
       "Absolute configuration has not been determined. Friedel pairs are merged;",
       "anomalous-dispersion terms are neglected. The displayed hand is arbitrary.",
       ";",
       f"_refine_diff_density_max {model['difference_map']['max_e_A3']:.5f}",
       f"_refine_diff_density_min {model['difference_map']['min_e_A3']:.5f}",
       f"_refine_diff_density_rms {model['difference_map']['rms_e_A3']:.5f}",
       "_refine_special_details",";",
       "Six non-H atoms have free coordinates and anisotropic displacement tensors.",
       "Seven H atoms ride on idealized geometry: N-H 0.91, methine C-H 1.00,",
       "methyl C-H 0.98 A. NH3 and CH3 torsions were optimized independently.",
       "Uiso(H)=1.5 Ueq(parent), except methine H: 1.2 Ueq(C).",
       "No heavy-atom distance or angle restraints were used.",
       "Coordinate and ADP standard uncertainties use the full numerical normal",
       "matrix multiplied by goodness-of-fit squared. Cell uncertainties include",
       "leave-one-scan-out variation, but not external detector calibration error.",
       "Scattering coefficients are explicitly recorded in refine.py.",
       "No external crystallographic validation software was used. The custom",
       "format and numerical checks are not a substitute for independent checkCIF.",
       ";",
       "_computing_data_collection 'Recorded by input experiment; no collection software invoked'",
       "_computing_cell_refinement 'Self-written detector-plane / Ewald least squares'",
       "_computing_data_reduction 'Self-written TY6 decoder, summation integration and scaling'",
       "_computing_structure_solution 'Self-written symmetry-constrained charge flipping'",
       "_computing_structure_refinement 'Self-written anisotropic F^2 least squares'",
       "_computing_publication_material 'Self-written CIF exporter and numerical checks'"]
    loop(l,["_space_group_symop_id","_space_group_symop_operation_xyz"],
         [[1,"'x,y,z'"],[2,"'x+1/2,-y+1/2,-z'"],
          [3,"'-x,y+1/2,-z+1/2'"],[4,"'-x+1/2,-y,z+1/2'"]])
    atoms=[]
    for i,(label,typ,x) in enumerate(zip(LABELS,TYPES,xyz)):
        atoms.append([label,typ,*[su(a,b) for a,b in zip(x,esd[3*i:3*i+3])],
                      su(ueq[i],us[i]),"Uani",1,"d"])
    for i,(label,x) in enumerate(zip(H_LABELS,hx)):
        atoms.append([label,"H",*[f"{a:.6f}" for a in x],
                      f"{ueq[H_PARENTS[i]]*H_FACTORS[i]:.6f}","Uiso",1,"calc"])
    loop(l,["_atom_site_label","_atom_site_type_symbol","_atom_site_fract_x",
            "_atom_site_fract_y","_atom_site_fract_z","_atom_site_U_iso_or_equiv",
            "_atom_site_adp_type","_atom_site_occupancy","_atom_site_calc_flag"],atoms)
    adps=[]
    for i,label in enumerate(LABELS):
        order=[0,1,2,5,4,3]
        adps.append([label,*[su(u[i,j],esd[18+6*i+j]) for j in order]])
    loop(l,["_atom_site_aniso_label","_atom_site_aniso_U_11","_atom_site_aniso_U_22",
            "_atom_site_aniso_U_33","_atom_site_aniso_U_12","_atom_site_aniso_U_13",
            "_atom_site_aniso_U_23"],adps)
    geom=validation["molecular_geometry"]
    loop(l,["_geom_bond_atom_site_label_1","_geom_bond_atom_site_label_2",
            "_geom_bond_distance","_geom_bond_site_symmetry_2"],
         [[b["atom1"],b["atom2"],su(b["distance"],b["esd"]),"."] for b in geom["bonds"]])
    loop(l,["_geom_angle_atom_site_label_1","_geom_angle_atom_site_label_2",
            "_geom_angle_atom_site_label_3","_geom_angle","_geom_angle_site_symmetry_1",
            "_geom_angle_site_symmetry_3"],
         [[a["atom1"],a["atom2"],a["atom3"],su(a["angle"],a["esd"]),".","."]
          for a in geom["angles"]])
    hb=[]
    for h in geom["hydrogen_bonds"]:
        sym=f"{h['symop']}_"+''.join(str(t+5) for t in h["translation"])
        hb.append([h["donor"],h["hydrogen"],h["acceptor"],"0.91",
                   f"{h['H_A']:.4f}",f"{h['D_A']:.4f}",f"{h['D_H_A']:.2f}",sym])
    loop(l,["_geom_hbond_atom_site_label_D","_geom_hbond_atom_site_label_H",
            "_geom_hbond_atom_site_label_A","_geom_hbond_distance_DH",
            "_geom_hbond_distance_HA","_geom_hbond_distance_DA","_geom_hbond_angle_DHA",
            "_geom_hbond_site_symmetry_A"],hb)
    (OUT/"alanine.cif").write_text("\n".join(l)+"\n",encoding="ascii")
    fcf=["data_alanine_reflections","_audit_creation_method 'Self-written F^2 calculation'",
         "_reflns_Friedel_coverage ?"]
    loop(fcf,["_refln_index_h","_refln_index_k","_refln_index_l","_refln_F_squared_meas",
              "_refln_F_squared_sigma","_refln_F_squared_calc","_refln_phase_calc"],
         [[int(r.h),int(r.k),int(r.l),f"{r.I:.8f}",f"{r.sigma:.8f}",
           f"{r.Fc2:.8f}",f"{np.degrees(r.phase_rad):.7f}"] for r in refs.itertuples()])
    (OUT/"alanine_reflections.fcf").write_text("\n".join(fcf)+"\n")
    (OUT/"alanine_merged.hkl").write_text("".join(
        f"{int(r.h):4d}{int(r.k):4d}{int(r.l):4d}{r.I:8.2f}{r.sigma:8.2f}\n"
        for r in refs.itertuples())+"   0   0   0    0.00    0.00\n")
    allxyz=np.concatenate([xyz,hx])
    alltypes=TYPES+["H"]*7
    mol=["13","One alanine molecular unit; crystallographic hand is not assigned."]
    mol.extend(f"{t} {x[0]:.7f} {x[1]:.7f} {x[2]:.7f}" for t,x in zip(alltypes,allxyz*cell))
    (OUT/"alanine_molecule.xyz").write_text("\n".join(mol)+"\n")
    uc=["52","P212121 cell; a,b,c = "+", ".join(map(str,cell))]
    for sg,tr in zip(SIGNS,TRANS):
        uc.extend(f"{t} {x[0]:.7f} {x[1]:.7f} {x[2]:.7f}"
                  for t,x in zip(alltypes,((allxyz*sg+tr)%1)*cell))
    (OUT/"alanine_unitcell.xyz").write_text("\n".join(uc)+"\n")
    return dict(volume=float(volume),density=float(4*89.093/(.602214076*volume)),
                final_nonshadowed_observations=len(obs),Rint_all_final_resolution=float(ri),
                theta_max=float(theta_max),theta_full=float(fulltheta))


def plots(model,validation):
    refs=pd.read_csv(OUT/"final_refinement_reflections.csv")
    obs=pd.read_csv(OUT/"merged_corrected_observations.csv")
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    ax=axes[0,0]
    ax.scatter(np.sqrt(np.maximum(refs.I,0)),np.sqrt(refs.Fc2),s=10,alpha=.7,color="#136c70")
    lim=np.sqrt(max(refs.I.max(),refs.Fc2.max()))*1.02
    ax.plot([0,lim],[0,lim],color="#6b6b6b",lw=1)
    ax.set(xlabel="Observed amplitude (relative)",ylabel="Calculated amplitude (relative)",
           title="Final refinement")
    ax=axes[0,1]
    residual=(refs.I-refs.Fc2)*np.sqrt(refs.weight)
    ax.scatter(1/refs.d,residual,s=10,color="#aa4352")
    ax.axhline(0,color="#777777",lw=1)
    ax.set(xlabel="1/d (1/angstrom)",ylabel="Weighted intensity residual",title="Resolution dependence")
    ax=axes[1,0]
    c=validation["completeness"]
    ax.plot([r["d_min"] for r in c[:6]],[r["fraction"]*100 for r in c[:6]],"o-",color="#426297")
    ax.set(xlabel="Minimum d (angstrom)",ylabel="Cumulative completeness (%)",ylim=(75,101))
    ax=axes[1,1]
    for run,g in obs.groupby("run"):
        g=g.sort_values("omega")
        ax.plot(g.omega,g.scale_factor,label=f"Scan {int(run)}",lw=1)
    ax.set(xlabel="Omega (degree)",ylabel="Empirical scale",title="Equivalent-reflection scaling")
    ax.legend(fontsize=8,ncol=2)
    fig.tight_layout()
    fig.savefig(OUT/"quality_diagnostics.png",dpi=180)
    plt.close(fig)
    xyz=np.array(model["xyz"])*model["cell"]
    hx=np.array(model["hydrogen_xyz"])*model["cell"]
    points=np.r_[xyz,hx]
    center=points.mean(axis=0)
    _,_,v=np.linalg.svd(points-center)
    screen=(points-center)@v.T
    fig,ax=plt.subplots(figsize=(7,5))
    for i,j in [(0,3),(1,3),(3,4),(2,4),(4,5)]:
        ax.plot(screen[[i,j],0],screen[[i,j],1],color="#777777",lw=3,zorder=1)
    for i,parent in enumerate(H_PARENTS):
        ax.plot(screen[[parent,6+i],0],screen[[parent,6+i],1],color="#aaaaaa",lw=1,zorder=1)
    colors={"O":"#be4452","N":"#487eae","C":"#414747","H":"#c5c9c9"}
    for i,(point,t) in enumerate(zip(screen,TYPES+["H"]*7)):
        ax.scatter(point[0],point[1],s=160 if t!="H" else 45,color=colors[t],zorder=3)
        if i<6:
            ax.annotate(LABELS[i],point[:2],xytext=(7,7),textcoords="offset points",fontsize=11)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Alanine: refined molecular geometry\nAbsolute configuration undetermined")
    fig.tight_layout()
    fig.savefig(OUT/"molecule.png",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    model=json.loads((OUT/"model_final.json").read_text())
    validation=json.loads((OUT/"validation.json").read_text())
    absorption=json.loads((OUT/"absorption.json").read_text())
    summary=write_cif(model,validation,absorption)
    (OUT/"export_statistics.json").write_text(json.dumps(summary,indent=2))
    plots(model,validation)
    print(json.dumps(summary))
