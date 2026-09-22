# Methods and Reproducibility

## Scope and Provenance

This is a self-written reduction, solution, refinement and checking pipeline.
Only the Python interpreter specified by `environment.json` was invoked.
Numerical dependencies are NumPy, SciPy, pandas and matplotlib; all other imports
are Python standard-library modules or source files written in this directory.
No crystallographic program, library, reference coordinates, external structure
database, downloaded table, or package source code was used.

The specimen identity and C3 H7 N O2 formula come from the task and supplied
experiment metadata. Z=4, the nominal wavelength, detector geometry, a starting
orientation matrix, scan axes and recorded temperature are input metadata.
The orientation matrix is an indexing seed, not a structural coordinate model.
The initial matrix was checked against diffraction peaks and subsequently fitted.
Atomic phases and positions were solved here, without a molecular reference model.

`input_provenance.json` identifies metadata files and the 616 formal raw frames.
The 30 pre-experiment frames were not needed and were not included in refinement.
The supplied crystal information also contains earlier cell estimates. Those
estimates were not substituted for the independently fitted final cell.

## Raw TY6 Decoder

`rawio.py` was inferred from the supplied binary files. It does not import or
reuse a detector-format reader. The 800 x 775 images have 6576-byte headers.
After a compressed-byte-count word, each image row has its own indexed interval.
Each row encodes one initial ordinary difference, blocks of sixteen differences,
and fifteen remaining differences. A block control byte specifies two widths
for groups of eight, in low-nibble/high-nibble order. Values are packed least
significant bits first and are offset by 2^(width-1)-1. Eight-bit escape values
carry signed 16- or 32-bit differences following that block. Cumulative sums
recover the integer photon counts independently for each row.

Assertions check the exact consumed length against every row pointer.
All 616 decoded frame totals agree with the corresponding header numeric means
times 620000 pixels to floating-point roundoff, at most 7.3e-12 counts.
The total is 22,210,409 counts; the maximum pixel is 12,764, below the input
overflow threshold. Six selected frames were independently decoded again and
compared exactly to cached arrays. All formal-frame SHA-256 fingerprints were
rechecked. Unused trailing bytes are not interpreted or needed to recover counts.

## Peak Locations and Geometry

`decode_all.py` locates 3 x 3 local summed maxima and measures background-subtracted
5 x 5 centroids. The first set contains 6,220 peak observations. Motor steps are
converted using the supplied header conventions: 12800 steps/degree for omega
and detector angle, 6400 steps/degree for kappa and phi. The actual frame headers
are used, not an assumed correspondence to forward/reverse collection order.

`find_conventions.py` enumerated axis, image and angle sign choices against the
supplied orientation seed. In the selected convention the incident beam is +X;
detector image axes and the three motor rotations have the recorded negative
sign convention, and the kappa axis tilts toward -X from +Z.

`fit_geometry.py` groups compatible indexed peaks across frames to estimate
intensity-weighted centroids. The original simplified module model is retained
as `geometry.json`. `geometry_flexible.py` uses 880 stronger interior centroids
to fit a common orthorhombic reciprocal metric, its orientation, both independent
detector module planes, kappa geometry and five relative scan omega offsets.
The fixed pixel pitch is 0.100 mm. The residual reciprocal-vector RMS is about
0.000332 in dimensionless wavelength-scaled units.

Scattering vectors are q = s_out - s_in. Indexing uses q0 = R^-1 q and h = UB^-1 q0.
The reciprocal columns have lengths wavelength/a, wavelength/b, wavelength/c.
The Ewald condition is q_X = -|q|^2/2. Rotation intersections are solved analytically
and rays intersect the separately calibrated detector planes.

Internal cell uncertainties use the larger of conditional normal-matrix
uncertainties and leave-one-scan-out jackknife estimates. These do not encompass
unknown external pixel-pitch, spectral-wavelength or detector calibration errors.
The final cell differs by approximately 0.06-0.09% from a supplied earlier cell.
No assertion of external cell accuracy is made.

## Integration and Correction

`integrate_flexible.py` enumerates all lattice points to d=0.70 A and predicts
3,895 interior observations. Detector edges, the module gap, rotation endpoints,
and very small rotation Lorentz denominators are excluded by geometry.
Integration uses a four-pixel-radius disk, a six-to-eight-pixel background annulus,
and a rotation half-window

    clamp(0.35 + 0.18 |q| / |q_Y|, 0.65, 2.8) degrees.

Several complete detector frames cover each window. Measured strong-peak
centroids may recenter a window within fixed bounds; weak peaks use predictions.
Unrelated bright pixels are clipped from background estimates only. The signal
itself is summed, never thresholded to positive intensity.
The Poisson variance includes signal counts and uncertainty in estimated
background. Negative net observations are retained.

The rotation correction is proportional to |q_Y|. Unpolarized beam polarization
is P = (1 + cos(2 theta)^2)/2, consistent with the supplied 0.5 factor.
All scans have 0.5-degree steps and 2-second exposures, so the common exposure
and scan-speed factors are absorbed into the overall scale.

`shadow_mask.py` caps individual diffuse-count contributions at two counts,
averages four temporal blocks per scan, smooths spatially and compares to the
radial diffuse-background median. A conservative minimum transmission ratio
below 0.55 within the signal footprint flags a shadow. This excludes 68 predicted
observations before merging. It does not use calculated structure factors.
All flags and original observations are retained in CSV.

`absorption.py` normalizes the nine supplied face normals, intersects plane
triplets to find the convex crystal, and samples its volume with fixed-seed Sobol
points. For each reflection, the incident and exit path lengths through the
polyhedron are found from the bounding planes. The average exp[-mu (t_in+t_out)]
is the transmission. mu=0.11620 mm^-1 is the value in the supplied CAP_shape header.
There are 29,407 in-crystal quadrature points. The half-sample convergence check
changes sampled transmissions by at most about 4.5e-6.
This assumes uniform illumination; the supplied crystal shape is not independently
remeasured. Numerical transmission is 0.989255-0.995397.

## Symmetry and Merging

The observed odd h00, 0k0 and 00l classes are not significantly above noise.
The metric, equivalence agreement and these screw absences support P 21 21 21.
The operations used are:

    x,y,z
    x+1/2,-y+1/2,-z
    -x,y+1/2,-z+1/2
    -x+1/2,-y,z+1/2

Equivalent intensities and Friedel pairs are merged as |h|,|k|,|l|.
Absences are exported separately and not used as independent allowed reflections.
Anomalous terms are neglected; neither a Flack nor a Hooft estimate is claimed.

`merge.py` alternates equivalent means with intensity-scale estimation.
Mean-dependent variance floors prevent zero-count shadows from dominating a
strong group by having artificially tiny counting uncertainties.
Discrepant equivalent observations are continuously downweighted, not silently
removed from the saved observation file. A diagnostic threshold of robust weight
greater than 0.5 distinguishes the quoted filtered Rint from the unfiltered one.
There are 51 low-weight observations after the physical shadow mask.

Per-scan quadratic omega scale functions were selected only after improving
held-out equivalent-reflection groups: filtered validation Rint decreased from
about 0.05115 for constant scan scales to 0.04887 for quadratic scales.
No structural Fc values are involved in this scale selection.
The reported unfiltered Rint at the final resolution is about 0.06763; the
robust-filtered diagnostic over the full reduced set is about 0.04073.

All 623 independent allowed reflections with d >= 0.75 A, including negatives,
enter the final refinement. This is 623/642 = 97.04% coverage; coverage to 0.80 A
is 540/540. The reduced table to about 0.70 A is also preserved, but its coverage
is only 649/788 = 82.36%. Missing reflections are never fabricated.

## Ab Initio Solution

`solve.py` uses fixed-seed, symmetry-constrained charge flipping on a 32 x 32 x 64
grid. Shell intensity normalization removes approximate radial falloff.
Unmeasured Fourier coefficients are set to zero in this implementation.
Twenty-four random starts and up to 900 iterations per start are recorded.
Symmetry-independent density maxima provide six non-H atom candidates.

`refine.py` first optimizes an all-carbon six-site model, then tries the 60
assignments of one N and two O among those sites. The best assignment forms a
carboxylate carbon bonded to two oxygens and a tetrahedral alpha carbon with N
and methyl substituents. No reference alanine coordinates were used.
`canonical_molecule` determines atom roles from periodic connectivity, so the
workflow does not depend on a particular peak ordering or origin.

The final solving input is preserved as `phasing_input.csv`; it contains no
phases or calculated structure factors. `solution_initial.json`,
`solution_density.npz` and `phasing_history.csv` record the solution stage.

## Refinement

`refine_aniso.py` evaluates complex structure factors directly over all four
symmetry operations. Neutral-atom four-Gaussian scattering coefficients are
explicit numerical constants in `refine.py`. They were not read from local
crystallographic software or downloaded. Independent tabulation verification
remains a limitation of this restricted workflow.

For each operation, both coordinates and the anisotropic tensor are transformed.
The displacement factor is exp[-2 pi^2 (h_i/a_i) U_ij (h_j/a_j)].
The six independent heavy atoms have 18 coordinate and 36 ADP parameters.
There is one overall intensity-scale parameter and two terminal H-group torsions:
57 parameters in total. All heavy-atom distances and angles are unrestrained.

Seven idealized riding H atoms are present. N-H=0.91 A, methine C-H=1.00 A,
methyl C-H=0.98 A. Uiso(H) is 1.5 Ueq(parent), or 1.2 for the methine hydrogen.
NH3 and CH3 torsions are optimized. The H positions are not independent observed
electron-density maxima and should not be described as freely refined.

Refinement minimizes sum w (Iobs-Icalc)^2 with

    w = 1/[sigma(Iobs)^2 + (0.03 P)^2]
    P = [max(Iobs,0) + 2 Icalc]/3.

Weights are updated in successive cycles. This coefficient was fixed, not tuned
to force goodness-of-fit to one. The final goodness-of-fit remains about 1.31.
No independent reflection is discarded on the basis of Fo/Fc disagreement.
Negative measured intensities remain in F^2 least squares; their amplitudes are
set to zero only when reporting conventional R1.

The full numerical normal matrix supplies covariance, scaled by goodness-of-fit
squared. Geometry uncertainties propagate this covariance and the adopted
internal cell uncertainties. Cell-coordinate cross covariances are unavailable.
These uncertainties do not capture all data-reduction or calibration systematics.

## Checks and Limits

- Raw frame row boundaries, totals, selected exact redecodes, and SHA-256 checks.
- Systematic absence intensities and equivalent-reflection agreement.
- Positive eigenvalues for every non-H anisotropic displacement tensor.
- Correct covalent geometry and three plausible N-H...O contacts.
- Difference Fourier synthesis over measured reflections only, with no invented
  missing amplitudes. The final grid is 64 x 64 x 128.
- Five-fold parameter-refinement cross-validation over all 623 reflections:
  combined held-out R1 is about 0.04097, versus 0.03781 in the full fit.
  Initial phases, space-group selection and reduction used the full data, so this
  is a conditional refinement check, not a completely blind pipeline validation.
- Leave-one-scan-out cell stability.
- Analytic isotropic amplitude derivatives checked by central differences.
- Friedel conjugacy and exact screw absences checked numerically.
- A self-written CIF tokenizer verifies loop dimensions and atom/reflection counts.
- A separately coded operation-by-operation structure-factor calculation reads
  the rounded exported CIF and checks it against the FCF and reported R values.

These checks are not an external CIF dictionary validation or a substitute for
checkCIF/PLATON. Absolute configuration, external cell calibration, crystal habit
and color, and independent instrument validation remain unresolved. The final
structure is checkable and has near-publication conventional residuals, but it
should not be represented as already independently publication-validated.

## Reproduction Order

`reproduce.py --full` runs the saved stages in this order:

1. decode_all.py
2. fit_geometry.py
3. geometry_flexible.py
4. integrate.py
5. shadow_mask.py
6. integrate_flexible.py
7. absorption.py
8. final_merge.py
9. solve.py
10. refine.py
11. refine_aniso.py --final
12. validate.py
13. export_results.py
14. check_outputs.py

It enforces the interpreter path in `environment.json`, fixes random seeds where
used, captures separate logs, and never writes to `inputs`.
Only generic linear algebra, optimization, FFT, image filtering, spatial convex
hulls and quadrature tools from the authorized libraries are used.
