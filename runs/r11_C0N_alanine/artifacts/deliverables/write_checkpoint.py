from pathlib import Path
import json

OUT=Path(__file__).resolve().parent
m=json.loads((OUT/"model_isotropic.json").read_text())
lines=["data_alanine_intermediate",
       "_audit_creation_method 'Self-written NumPy/SciPy pipeline'",
       "_chemical_formula_sum 'C3 H7 N O2'",
       "_cell_length_a %.7f"%m["cell"][0],
       "_cell_length_b %.7f"%m["cell"][1],
       "_cell_length_c %.7f"%m["cell"][2],
       "_cell_angle_alpha 90",
       "_cell_angle_beta 90",
       "_cell_angle_gamma 90",
       "_cell_formula_units_Z 4",
       "_space_group_name_H-M_alt 'P 21 21 21'",
       "_space_group_IT_number 19",
       "_diffrn_radiation_wavelength 0.71073",
       "_diffrn_ambient_temperature 150",
       "_refine_ls_R_factor_all %.6f"%m["R1"],
       "_refine_special_details",
       ";",
       "INTERMEDIATE MODEL: non-H atoms only; no claim of publication quality.",
       "Hydrogens, anisotropic displacement and final validation are pending.",
       "Absolute structure is unknown. Cell is an approximate self-fitted metric.",
       ";",
       "loop_","_space_group_symop_id","_space_group_symop_operation_xyz",
       "1 'x,y,z'","2 'x+1/2,-y+1/2,-z'",
       "3 '-x,y+1/2,-z+1/2'","4 '-x+1/2,-y,z+1/2'",
       "loop_","_atom_site_label","_atom_site_type_symbol","_atom_site_fract_x",
       "_atom_site_fract_y","_atom_site_fract_z","_atom_site_U_iso_or_equiv",
       "_atom_site_occupancy"]
count={}
for typ,xyz,u in zip(m["types"],m["xyz"],m["uiso"]):
    count[typ]=count.get(typ,0)+1
    lines.append("%s%d %s %.7f %.7f %.7f %.6f 1"%(typ,count[typ],typ,*xyz,u))
(OUT/"alanine_intermediate.cif").write_text("\n".join(lines)+"\n")
