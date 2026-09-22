"""Remove unsupported experimental defaults without altering refined numbers."""
import argparse
from pathlib import Path
import gemmi

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()
doc = gemmi.cif.read_file(args.input)
block = doc.sole_block()
for tag in ("_cell_measurement_temperature", "_diffrn_ambient_temperature"):
    block.set_pair(tag, "?")
block.set_pair("_chemical_name_common",
               gemmi.cif.quote("NU-1000 framework-only average structure"))
block.set_pair("_chemical_formula_moiety", "?")
block.set_pair("_computing_structure_solution", gemmi.cif.quote("SHELXT 2018/2"))
block.set_pair("_computing_structure_refinement", gemmi.cif.quote("SHELXL 2019/3"))
block.set_pair("_computing_publication_material",
               gemmi.cif.quote("cctbx, smtbx, gemmi; supplied scripts"))
block.set_pair("_atom_sites_solution_primary", "dual")
block.set_pair("_atom_sites_solution_secondary", "difmap")
details = """\
Independent solution from the supplied HKLF4 intensities, unit cell and
synthesis diagram. No reference structure or sample-specific external
information was used. SHELXT phasing used d >= 1.2 A; the final SHELXL
full-matrix refinement used d >= 1.0 A. No individual reflections were
manually omitted. All original observations and negative intensities are
retained in the separate HKL and embedded archive.

The P6/mmm average framework has 20 independent non-H sites (2 Zr, 6 O,
12 C); all are anisotropic with full physical occupancy, including the
appropriate SHELX special-position occupancy factors. Six aromatic H
sites were placed geometrically and ride on their parent carbon atoms,
C-H = 0.93 A, Uiso(H) = 1.2 Ueq(C). The 0.93 A distance is a modeling
convention, not evidence for a measurement temperature. No node-bound
hydrogens were located or added. O3 and O6 are triply bridging node O;
O4 and O5 are terminal node-coordinated O. Their protonation and possible
terminal ligand identities are not established.

The formula C88 H44 O32 Zr6 and derived formula weight, density and F(000)
describe coordinate-model atoms only. They exclude unknown node-bound
hydrogens and the unassigned pore content, and must not be used as the
complete chemical or analytical formula of the specimen.

DFIX/DANG/FLAT restraints regularize the aromatic and carboxylate
geometry; RIGU/SIMU restrain ligand displacement parameters. All
restraints and uncertainties are explicitly listed in the embedded
RES and separate INS. A no-geometrical-restraint control retained the
same connectivity with a maximum non-H shift of about 0.035 A.

No explicit pore solvent/guest atoms were assigned. Unassigned pore
scattering was estimated with the installed smtbx.masks.mask algorithm
and supplied to SHELXL as fixed partial structure factors via ABIN.
Probe radius 1.3 A, shrink radius 1.3 A, grid step 0.25 A, grid
160 x 160 x 72, set completion enabled, maximum 100 cycles.
The input source map and the exact FAB are archived. The mask volume
is approximately 17381 A^3 (78.67 percent of the cell). Its electron
count is processing-dependent and is not assigned a solvent formula.
The mask is not treated as a collection of independently observed atoms;
formal coordinate standard uncertainties are conditional on the chosen
mask and restraints. No independent Rfree is claimed.

Wavelength-specific Sasaki dispersion factors from the installed cctbx
tables were used (Zr f-prime = -9.0413504, f-double-prime = 2.7716699).
A Henke-table control changed non-H coordinates by up to about 0.030 A;
the given wavelength is sensitive to the Zr-edge dispersion model.
The tabulated f-double-prime was used to calculate atomic absorption
cross sections, not to invent an experimental absorption correction.

Limitations: high intensity-merging residuals, weak high-resolution
data, several large/elongated positive-definite displacement ellipsoids,
and residual low-angle mismatches remain. No split-site model is asserted.
Temperature, crystal dimensions, cell uncertainties, raw images and
data-reduction/absorption-correction history were not supplied.
Default 293 K temperature fields have deliberately been replaced by ?.
This is a reproducible framework model, not a fully validated
publication-ready determination of the complete solvated specimen.
"""
block.set_pair("_refine_special_details", gemmi.cif.quote(details))
block.set_pair("_exptl_special_details", gemmi.cif.quote(
    "Only integrated HKLF4 intensities, an initial INS and a synthesis "
    "diagram were supplied. Measurement temperature and data-reduction "
    "history are unknown. No new experiment was performed."))
doc.write_file(args.output)
gemmi.cif.read_file(args.output)
print("Curated and syntax-checked:", args.output)
