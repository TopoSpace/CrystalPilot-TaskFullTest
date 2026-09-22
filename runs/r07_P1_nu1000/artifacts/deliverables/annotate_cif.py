"""Remove invented/default experimental metadata after a SHELXL reproduction.
Only reporting metadata are changed; atom sites, ADPs, reflections and refinement
statistics are not changed. Requires gemmi. Usage: python annotate_cif.py file.cif
"""
from pathlib import Path
import sys,datetime
import gemmi
DETAILS='''Framework-only average structure solved from the supplied start.hkl and
start.ins, with element identity and linker chemistry inferred only from the
supplied synthesis image. No reference structure was used.

P6/mmm; wavelength 0.68883 A; refinement against F-squared to d_min=1.10 A.
All observations within this resolution range, including weak and negative
intensities, were retained; no individual reflection OMIT instructions were used.
Non-H atoms were refined anisotropically. Aromatic C-H atoms were placed in
riding positions (C-H 0.93 A, Uiso(H)=1.2Ueq(C)). Node-bound hydrogen atoms were
not located or assigned. The distinction between node oxide, hydroxide and
coordinated water, and the full framework proton count, remain unresolved.
The reported formula C88 H44 O32 Zr6, formula weight, F(000), calculated density
and absorption coefficient refer ONLY to the explicit coordinate model.
They exclude unassigned node hydrogen atoms and all pore contents and must not
be interpreted as an experimentally established bulk chemical composition.

Weak ligand scattering required soft aromatic bond-length/1,3-distance and
planarity restraints, carboxylate distance similarity, RIGU and SIMU restraints.
The exact numerical instructions are preserved in the embedded RES and the
accompanying INS. Some large/anisotropic displacements remain; no unsupported
split-site disorder model was introduced.

Pore scattering was included as additive A/B structure-factor contributions
from PLATON SQUEEZE (supplied FAB), while the original measured HKL was retained.
No solvent or guest molecule was assigned. The solvent-accessible volume is
approximately 79 percent of the cell; four low-order independent reflections
are missing. Both this very large void fraction and the limited effective
resolution are outside the routine ideal conditions described by the SQUEEZE
manual. Mask-dependent residuals and recovered electron counts are therefore
not independent validation of a unique pore model or solvent stoichiometry.
Unmasked models and mask/resolution/scattering-factor sensitivity trials are
retained in the calculation records.

Anomalous dispersion was calculated at the input wavelength by the
Cromer-Liberman method in Gemmi 0.7.5; attenuation coefficients came from cctbx.
This wavelength lies near the Zr K edge. The actual experimental energy
calibration and sample-specific anomalous scattering are unknown. Alternative
Henke and Sasaki tabulations were examined and archived, not treated as
additional experimental information.

Experimental temperature, crystal dimensions, instrument, absorption-correction
history, cell estimated standard deviations and collection protocol were not
provided. In particular, SHELXL's automatically inserted 293(2) K values have
been replaced with unknown values. Synthesis temperatures are not diffraction
measurement temperatures. This is a reproducible framework model with explicit
limitations, not an unqualified publication-ready experimental report.'''

def annotate(path):
    path=Path(path); doc=gemmi.cif.read_file(str(path)); block=doc.sole_block()
    for tag in ['_cell_measurement_temperature','_diffrn_ambient_temperature']:
        block.set_pair(tag,'?')
    block.set_pair('_chemical_name_common',gemmi.cif.quote('NU-1000 (framework-only average model)'))
    block.set_pair('_computing_structure_solution',gemmi.cif.quote('SHELXT 2018/2'))
    block.set_pair('_computing_structure_refinement',gemmi.cif.quote('SHELXL 2019/3'))
    block.set_pair('_computing_publication_material',gemmi.cif.quote('Gemmi 0.7.5; Python metadata annotation'))
    block.set_pair('_audit_creation_method',gemmi.cif.quote('SHELXL output; unknown experimental metadata explicitly annotated'))
    block.set_pair('_audit_creation_date',datetime.date.today().isoformat())
    block.set_pair('_refine_special_details',';\n'+DETAILS+'\n;')
    doc.write_file(str(path))
    # Round-trip syntax validation is intentionally part of the reproduction.
    gemmi.cif.read_file(str(path))
    return path
if __name__=='__main__':
    print('Annotated:',annotate(sys.argv[1]))
