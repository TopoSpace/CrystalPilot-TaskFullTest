"""Local covariance extension for DIALS equal-shift constraints.
Does NOT modify the installed software. DIALS 3.30 skips ESDs whenever equal-shift
constraints are used (engine.py, AdaptLstbx.calculate_esds). Here the same final
constrained normal matrix is inverted, then expanded with V_full = T V_reduced T^T.
The scale is DIALS' own final reduced chi-square. Refinement steps are unchanged.
The linear map T is exact for the equal-shift constraints, not a numerical fit.
"""
import json
from pathlib import Path
import numpy as np
from scitbx.array_family import flex
from dials.algorithms.refinement import engine
from dials.command_line.refine import run

original = engine.AdaptLstbx.calculate_esds

def constrained_esds(self):
    cm = self._constr_manager
    if cm is None:
        return original(self)
    if self.history.get_nrows() == 0 or self.cf is None:
        raise RuntimeError('No final normal matrix for uncertainty calculation')
    cf_inv = self.cf.matrix_packed_u_as_upper_triangle().matrix_inversion()
    inv_normal = cf_inv.matrix_multiply_transpose(cf_inv).as_numpy_array()
    chi2 = float(self.history['reduced_chi_squared'][-1])
    reduced_cov = chi2 * inv_normal
    T = np.zeros((cm._n_full_params, reduced_cov.shape[0]))
    for column, row in enumerate(cm._unconstrained_idx):
        T[row, column] = 1.0
    for column, rows in enumerate(cm._constrained_gps, start=cm._n_unconstrained_params):
        T[np.asarray(list(rows), dtype=int), column] = 1.0
    assert np.all(T.sum(axis=1) == 1)
    assert np.linalg.matrix_rank(T) == T.shape[1]
    full_cov = T @ reduced_cov @ T.T
    assert np.allclose(full_cov, full_cov.T, atol=1e-10)
    assert np.all(np.diag(full_cov) >= 0)
    np.savez('common_cell_covariance.npz', reduced_covariance=reduced_cov, expansion_matrix=T, full_covariance=full_cov)
    Path('common_cell_covariance.json').write_text(json.dumps({'method': 'V_full = T * (reduced_chi_square * inverse_normal_matrix) * T_transpose', 'reduced_chi_square': chi2, 'expanded_parameters': T.shape[0], 'independent_parameters': T.shape[1], 'constraint_groups': [list(g) for g in cm._constrained_gps], 'reduced_covariance_min_eigenvalue': float(np.linalg.eigvalsh(reduced_cov).min())}, indent=2), encoding='utf-8')
    self.parameter_var_cov = flex.double(full_cov)
    self._parameters.calculate_model_state_uncertainties(self.parameter_var_cov)
    self._parameters.set_param_esds(flex.sqrt(self.parameter_var_cov.matrix_diagonal()))
    print('Computed constrained covariance by exact linear propagation:', T.shape)

if __name__ == '__main__':
    engine.AdaptLstbx.calculate_esds = constrained_esds
    run(['refined.expt', 'refined.refl', 'joint_cell.phil', 'set_scan_varying_errors=True', 'output.experiments=refined_joint_esd.expt', 'output.reflections=refined_joint_esd.refl', 'output.log=refined_joint_esd.log'])
