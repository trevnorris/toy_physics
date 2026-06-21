# Path-A C0d Maxwell Gauge Identification

Verdict: **MIXED_GAUGE_PLUS_RESIDUAL**

## Scope

Diagnosis only: the saved C0b matrix is loaded and its SVD modes are recomputed. No gauge fix, deflation, recrawl, or changed-xi Jacobian assembly is performed.

## Operator Sources

```yaml
field_layout: {'source': 'stage1_solver.coupled_branch.pack_closed_coupled_fields/unpack_closed_coupled_fields', 'order': ['psi_real', 'psi_imag', 'a0', 'ar', 'aw', 'r0', 'mu'], 'spatial_a_lanes': ['ar', 'aw'], 'scalar_potential_lane': 'a0', 'gradient_embedding': 'zeros in psi_real, psi_imag, a0, r0, mu; d_r(lambda) in ar; d_w(lambda) in aw'}
discrete_gradient: {'source': 'stage1_solver.operators.tensor_center_gradient_r/w', 'boundary_closure': "the operators' one-sided centered-gradient closures; no extra lambda or A boundary condition imposed"}
raw_divergence: {'source': 'stage1_solver.operators.axisymmetric_vector_divergence(ar, aw)', 'definition': 'cell-average r^-2 d_r(r^2 A_r) + d_w A_w'}
weighted_gauge_residual: {'source': 'stage1_solver.operators.localized_maxwell_operator ar/aw gauge-control blocks', 'object': 'grad(Z(w) * axisymmetric_vector_divergence(ar, aw))', 'weight_source': 'stage1_solver.coupled_branch.localization_weight_torch'}
```

## Gauge Subspace

```yaml
status: MEASURED
construction: Full nodal scalar basis lambda_k, one scalar field per (r,w) cell. Each column is embedded as delta a0=0, delta ar=tensor_center_gradient_r(lambda_k), delta aw=tensor_center_gradient_w(lambda_k), with all non-A closed-state lanes zero. G is the SVD-orthonormalized column space. G_harm is formed from right singular directions of Z*D_A restricted to G whose weighted-divergence singular value is below the configured relative threshold.
grid: [16, 16]
scalar_basis_dim: 256
full_state_dim: 1297
dim_G: 255
dim_G_harm: 2
gradient_matrix_shape: [1297, 256]
gradient_rank_rtol: 1e-12
gradient_rank_threshold: 3.424787275790426e-11
weighted_divergence_operator: Z(w) * axisymmetric_vector_divergence restricted to G
harmonic_weighted_divergence_rtol: 1e-06
weighted_divergence_harmonic_threshold: 2.3590123904374058e-05
weighted_divergence_harmonic_indices: [253, 254]
```

## Controls

```yaml
status: PASS
positive_control: {'type': 'held_out_smooth_discrete_gradient', 'p_g_fraction': 1.0, 'non_a_remainder': 2.220446049250313e-16, 'pass': True}
negative_control: {'type': 'transverse_orthogonalized_random_A', 'p_g_fraction': 5.876660946583589e-30, 'non_a_remainder': 0.0, 'pass': True}
thresholds: {'positive_p_g_min': 0.9999999999, 'positive_non_a_max': 1e-12, 'negative_p_g_max': 1e-10, 'negative_non_a_max': 1e-12}
```

## Maxwell-Lane Mode Gate

A mode is `MAXWELL_GAUGE` iff `||P_G v||^2 >= 0.9` and `||grad(Z*D_A*A[v])||/||A[v]|| <= 0.1`. Raw `||D_A*A[v]||/||A[v]||` is reported for context but is not the gate.

| mode | sigma | P_G | P_G_harm | residual | weighted | raw_div | ar | aw | non_A | class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.700699e-11 | 9.995140e-01 | 9.971544e-01 | 4.860324e-04 | 2.507593e-01 | 1.487354e-01 | 7.729896e-03 | 9.922686e-01 | 1.499120e-06 | GENUINE_MAXWELL_STIFFNESS |
| 2 | 2.622571e-10 | 8.298041e-01 | 2.916985e-03 | 1.701959e-01 | 4.688636e+00 | 2.785336e+00 | 9.562146e-01 | 4.377191e-02 | 1.351818e-05 | GENUINE_MAXWELL_STIFFNESS |
| 3 | 5.926281e-10 | 9.998625e-01 | 2.631188e-01 | 1.375026e-04 | 1.219922e-01 | 3.282671e+00 | 9.864547e-01 | 1.352430e-02 | 2.101846e-05 | GENUINE_MAXWELL_STIFFNESS |
| 4 | 2.811268e-08 | 9.989308e-01 | 1.940663e-01 | 1.069190e-03 | 4.016166e-02 | 1.284280e-01 | 2.223673e-01 | 7.776291e-01 | 3.625923e-06 | MAXWELL_GAUGE |

## Full Cluster Context

| mode | sigma | spatial_A | P_G | weighted | raw_div | class |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 7.434694e-14 | 6.141704e-09 | 6.141699e-09 | 1.072884e-02 | 4.928965e-01 | NOT_A_MAXWELL_LANE_MODE |
| 1 | 4.700699e-11 | 9.999985e-01 | 9.995140e-01 | 2.507593e-01 | 1.487354e-01 | GENUINE_MAXWELL_STIFFNESS |
| 2 | 2.622571e-10 | 9.999865e-01 | 8.298041e-01 | 4.688636e+00 | 2.785336e+00 | GENUINE_MAXWELL_STIFFNESS |
| 3 | 5.926281e-10 | 9.999790e-01 | 9.998625e-01 | 1.219922e-01 | 3.282671e+00 | GENUINE_MAXWELL_STIFFNESS |
| 4 | 2.811268e-08 | 9.999964e-01 | 9.989308e-01 | 4.016166e-02 | 1.284280e-01 | MAXWELL_GAUGE |

## Gauge vs Stiffness Discriminator

Saved-matrix evidence distinguishes penalized residual gauge from genuine Maxwell stiffness as follows: high projection into the proper multi-lambda gradient range plus small weighted gauge residual is gauge; low `P_G` or large weighted residual is classified as Maxwell stiffness under the C0d gate.

## Verdict Support

```yaml
thresholds: {'p_g_min': 0.9, 'weighted_gauge_residual_max': 0.1, 'maxwell_lane_fraction_min': 0.9}
controls_status: PASS
mode_count: 4
required_mode_count: 4
mode_classifications:
  - {'mode_index': 1, 'sigma': 4.7006989691328044e-11, 'p_g_fraction': 0.9995139675622529, 'p_g_harm_fraction': 0.9971544035997102, 'weighted_gauge_residual': 0.2507592599846196, 'raw_divergence': 0.14873541668580761, 'spatial_a_energy_fraction': 0.9999985008801926, 'non_a_remainder': 1.4991198074021383e-06, 'classification': 'GENUINE_MAXWELL_STIFFNESS'}
  - {'mode_index': 2, 'sigma': 2.6225709826789087e-10, 'p_g_fraction': 0.829804081276686, 'p_g_harm_fraction': 0.0029169854810530875, 'weighted_gauge_residual': 4.68863599534441, 'raw_divergence': 2.785335556924372, 'spatial_a_energy_fraction': 0.999986481819913, 'non_a_remainder': 1.3518180086991016e-05, 'classification': 'GENUINE_MAXWELL_STIFFNESS'}
  - {'mode_index': 3, 'sigma': 5.926280559311367e-10, 'p_g_fraction': 0.9998624973765375, 'p_g_harm_fraction': 0.2631187529404603, 'weighted_gauge_residual': 0.12199215262929851, 'raw_divergence': 3.2826711927454153, 'spatial_a_energy_fraction': 0.9999789815431896, 'non_a_remainder': 2.1018456810351083e-05, 'classification': 'GENUINE_MAXWELL_STIFFNESS'}
  - {'mode_index': 4, 'sigma': 2.8112679473453323e-08, 'p_g_fraction': 0.9989308103957472, 'p_g_harm_fraction': 0.19406629929577604, 'weighted_gauge_residual': 0.04016165859397893, 'raw_divergence': 0.1284279956506222, 'spatial_a_energy_fraction': 0.9999963740765402, 'non_a_remainder': 3.6259234598157164e-06, 'classification': 'MAXWELL_GAUGE'}
```

## Recommended Next Step

Next step: do not apply a whole-wall gauge conclusion. Pin/deflate only the gated gauge modes and investigate the remaining Maxwell-lane residual modes with their lane split and weighted residual evidence.

Combined-fix design if `WALL_IS_ALL_GAUGE`:

Pin/deflate the global U(1) phase generator already confirmed by C0c, and deflate the A-sector discrete-gradient subspace G or add an equivalent compatible weighted gauge constraint for grad(Z*D_A*A); then rerun the C0 crawl. This C0d run does not implement that fix.

## Guard Confirmation

```yaml
faithful_operator_boundary: {'physical_residual_entry_point': 'stage1_solver.coupled_branch.patha_closed_branch_residual', 'operators_modified_by_c0': False, 'frozen_physics_modified_by_c0': False, 'physical_export_guard_modified_by_c0': False}
scope_guard: {'diagnosis_only': True, 'gauge_fix_or_deflation_implemented': False, 'recrawl_implemented': False, 'xi_reassembly_implemented': False, 'faithful_operators_touched_by_c0d': False, 'frozen_physics_touched_by_c0d': False, 'physical_export_guard_touched_by_c0d': False}
```

Machine artifact: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0d_maxwell_gauge/pathA_C0d_maxwell_gauge_identification.json`.
