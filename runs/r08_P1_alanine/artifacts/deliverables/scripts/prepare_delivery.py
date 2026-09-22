"""Create the deliverable structure and enrich only supported CIF metadata.
Calculated coordinates, displacement parameters, refinement statistics and
embedded SHELXL RES/HKL text are not changed. Original CIF is retained.
"""
import json
import shutil
import math
import gemmi
import numpy as np
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from run_command import ROOT

out=ROOT/'deliverables'
source=ROOT/'work/refine/10_publication'
backup=ROOT/'work/validation/pre_LIST4_delivery'
assert out.is_dir()
def check_target(target):
    if target.exists():
        previous=backup/target.name
        assert previous.exists() and previous.read_bytes()==target.read_bytes(), target
for extension in ['ins','res','hkl','fcf','lst']:
    target=out/f'alanine.{extension}'
    check_target(target)
    shutil.copy2(source/f'alanine.{extension}',target)
provenance=out/'provenance'
provenance.mkdir(exist_ok=True)
shutil.copy2(source/'alanine.cif',provenance/'alanine_shelxl_original.cif')
experiments=ExperimentListFactory.from_json_file(str(ROOT/'work/dials/refined_joint_esd.expt'),check_format=False)
crystal=experiments[0].crystal
reflections=flex.reflection_table.from_file(str(ROOT/'work/dials/refined_joint_esd.refl'))
used=reflections.select(reflections.get_flags(reflections.flags.used_in_refinement))
spacings=np.array(list(crystal.get_unit_cell().d(used['miller_index'])))
theta=np.degrees(np.arcsin(0.71073/(2*spacings)))
cell_measurement={'reflections':len(used),'theta_min':float(theta.min()),'theta_max':float(theta.max()),'cell':crystal.get_unit_cell().parameters(),'cell_esd':crystal.get_cell_parameter_sd(),'volume':crystal.get_unit_cell().volume(),'volume_esd_full_covariance':crystal.get_cell_volume_sd()}
(provenance/'cell_measurement.json').write_text(json.dumps(cell_measurement,indent=2),encoding='utf-8')
doc=gemmi.cif.read_file(str(source/'alanine.cif'))
b=doc.sole_block()
changes={}
def setvalue(tag,value):
    changes[tag]={'previous':b.find_value(tag),'new':value}
    b.set_pair(tag,value)
def text(tag,value):
    setvalue(tag,gemmi.cif.quote(value))
text('_audit_creation_method','SHELXL-2019/3; metadata added by prepare_delivery.py')
setvalue('_audit_creation_date','2026-09-22')
text('_chemical_name_systematic','2-azaniumylpropanoate')
text('_chemical_name_common','alanine')
text('_chemical_formula_moiety','C3 H7 N O2')
setvalue('_chemical_absolute_configuration','unk')
setvalue('_cell_volume','424.89(15)')
setvalue('_cell_measurement_temperature','150')
setvalue('_diffrn_ambient_temperature','150')
setvalue('_cell_measurement_reflns_used',str(len(used)))
setvalue('_cell_measurement_theta_min',f'{theta.min():.3f}')
setvalue('_cell_measurement_theta_max',f'{theta.max():.3f}')
setvalue('_exptl_crystal_density_method',gemmi.cif.quote('not measured'))
setvalue('_exptl_absorpt_correction_type','multi-scan')
text('_exptl_absorpt_process_details',
     'DIALS 3.30 physical model: scan-dependent scale and decay, low-order spherical-harmonic absorption correction; scale_interval=30 deg, decay_interval=60 deg, anomalous=True. No analytical face-based absorption correction was applied.')
text('_exptl_absorpt_special_details',
     'Absolute transmission factors are not available from this relative empirical correction. Crystal dimensions were not established. Instrument shadow masks were derived from persistent low raw-pixel background, not Fo-Fc selection; their definitions and unmasked controls are retained.')
text('_diffrn_source','microfocus X-ray tube (input generator record)')
setvalue('_diffrn_source_voltage','50')
setvalue('_diffrn_source_current','1')
text('_diffrn_measurement_device_type','Rigaku/Oxford Diffraction instrument with HyPix-Arc 100 detector')
text('_diffrn_detector','hybrid pixel area detector, two panels')
text('_diffrn_measurement_method','omega scans at multiple kappa/phi settings')
setvalue('_diffrn_detector_area_resol_mean','10')
text('_diffrn_special_details',
     '616 production frames in six scans, 0.5 deg per frame, 2 s per frame; 30 pre-experiment frames excluded. Collection 2024-02-09. Temperature monitor range 149.93-150.05 K; a separate parameter file lists 149 K. The monitor record is used, without inventing a temperature standard uncertainty. Wavelength 0.71073 A is from frame headers. Nominal detector distance 47 mm; panel geometry was refined independently. Legacy detector/monochromator description strings in the parameter file are not treated as reliable hardware identification.')
text('_computing_data_collection','CrysAlisPro 42.90a (supplied input metadata)')
text('_computing_cell_refinement','DIALS 3.30; scan-varying orientations and constrained common orthorhombic cell')
text('_computing_data_reduction','DIALS 3.30; local raw-pixel shadow masks')
text('_computing_structure_solution','SHELXT 2018/2 (intrinsic phasing)')
text('_computing_structure_refinement','SHELXL-2019/3 (Sheldrick, 2015)')
text('_computing_publication_material','SHELXL-2019/3; Gemmi and local metadata/validation scripts')
setvalue('_atom_sites_solution_primary','dual')
setvalue('_atom_sites_solution_secondary','difmap')
text('_refine_ls_abs_structure_details',
     'Parsons quotient estimate from 339 selected quotients; see SHELXL listing. The standard uncertainty is too large to determine absolute structure. The deposited hand is arbitrary; no L/D or R/S assignment and no inversion-twin fraction are established from these data.')
text('_refine_special_details',
     'Non-H atoms were refined anisotropically. Three ammonium H atoms are supported by an N-H omit difference map and refined with independent coordinates and isotropic displacement parameters without restraints. C-H atoms ride on C; the methyl group rotates, with Uiso(H)=1.5Ueq(C), and the methine H has Uiso(H)=1.2Ueq(C). The refinement uses d>=0.77 A because the Laue completeness is 99.83% to this limit; the sparse measured extension to about 0.708 A is retained in HKL and tested separately. No individual reflection OMIT instructions, extinction correction, solvent mask or twin model are used. Standard integration/partiality/scaling filters are recorded in the reduction logs. DIALS exact common-cell constraints suppress its default ESD output; a local covariance extension expanded the reduced normal-matrix covariance through the exact constraint Jacobian. This did not change fitted cell parameters or installed software. The cell-volume uncertainty 0.15 A^3 includes full cell covariance, whereas the original SHELXL value 0.10 A^3 neglects general cell correlations. Full covariance and code are supplied. Cell/atomic uncertainties do not include all possible instrument systematics.')
check_target(out/'alanine.cif')
doc.write_file(str(out/'alanine.cif'))
(provenance/'cif_metadata_changes.json').write_text(json.dumps(changes,indent=2),encoding='utf-8')
# Directly usable final reduction outputs; all earlier controls are archived later.
processing=out/'processing'
processing.mkdir(exist_ok=True)
for name in ['masked_input.expt','integrated.expt','integrated.refl','scaled.expt','scaled.refl','scaled_unmerged.mtz','scaled_merged.mtz','mask_definition.json','intensity_audit.json','intensity_audit.txt','geometry_0.mask','geometry_1.mask','geometry_2.mask','geometry_3.mask','dials.integrate.log','dials.symmetry.log','dials.scale.log']:
    path=ROOT/'work/masked_data'/name
    if path.exists():
        shutil.copy2(path,processing/name)
for name in ['refined_joint_esd.expt','refined_joint_esd.refl','refined_joint_esd.log','common_cell_covariance.npz','common_cell_covariance.json','joint_cell.phil','scale.phil']:
    shutil.copy2(ROOT/'work/dials'/name,processing/name)
for name in ['raw_frames_manifest.json','environment.json','commands.jsonl']:
    shutil.copy2(ROOT/name,provenance/name)
print(json.dumps(cell_measurement,indent=2))
print('Prepared final SHELXL files, enriched CIF, and final reduction outputs.')
