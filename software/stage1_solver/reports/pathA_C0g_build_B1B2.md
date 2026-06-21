# Path-A C0g Build B-1/B-2/B-4 Staged Report

> **SUPERSEDED-IN-PART (2026-06-21):** this `FOLD_DISSOLVED` correctly refutes the shallow "fold at τ≈0.0291233" (the crawl
> crosses it), but the subsequent deeper crawl (`reports/pathA_C0g_deepcrawl.md`, verdict `STILL_GOING`) found a NEW
> genuinely-real near-singularity ~1e-5 deeper (~τ=0.0291132): σ_min collapses 3.95e-4→7.03e-6 with NO tangent reversal.
> The fold-vs-conditioning question is NOT settled — it moved deeper; pseudo-arclength (B-3) is the earned tool to resolve
> it. Read with the deepcrawl report. Resume: `decisions/13` §0 (4).

B-2 verdict: **FOLD_DISSOLVED**
B-4 match gate: **PASS**

## Scope

- Built B-1 gauge-fix only in solver/path coordinates.
- Built B-4 analytic/sparse Jacobian assembly and verified it before B-2.
- Ran B-2 re-confirmation on the gauge-fixed crawl.
- B-3 pseudo-arclength was not built.
- The original `patha_closed_branch_residual` remains the convergence/progress arbiter.

## B-4 Sparse Jacobian Match Gate

```yaml
status: PASS
reference_tau: 0.029125
tolerances: {'max_abs': 1e-10, 'max_rel': 1e-08}
coverage: {'random_jv_probe_count': 5, 'random_jv_blocks': ['psi_real_rows', 'psi_imag_rows', 'a0_rows', 'ar_rows', 'aw_rows', 'wall_rows', 'mass_row'], 'explicit_dense_checks': ['wall_rows_full', 'mass_row_full', 'mu_column_full'], 'random_jv_reference': 'torch.func.jvp(original_residual)', 'dense_row_reference': 'torch.func.jacrev(original_residual)'}
assembly_speedup_x: 74.32453792381257
colored_seconds: 40.03063162486069
fast_seconds: 0.5385924049187452
```

Random-Jv per-block mismatch:

```yaml
{'psi_real_rows': {'max_abs': 7.105427357601002e-15, 'max_rel': 4.1929963340803407e-16, 'reference_linf': 20.745154802468363}, 'psi_imag_rows': {'max_abs': 8.881784197001252e-15, 'max_rel': 5.001307728287879e-16, 'reference_linf': 22.971977262556276}, 'a0_rows': {'max_abs': 1.0658141036401503e-14, 'max_rel': 3.284430376426575e-16, 'reference_linf': 35.94438175496331}, 'ar_rows': {'max_abs': 1.0658141036401503e-14, 'max_rel': 4.4385238068499245e-16, 'reference_linf': 33.51443033826936}, 'aw_rows': {'max_abs': 1.0658141036401503e-14, 'max_rel': 3.96629014285476e-16, 'reference_linf': 30.599733149393906}, 'wall_rows': {'max_abs': 1.1102230246251565e-16, 'max_rel': 1.1102230246251565e-16, 'reference_linf': 0.5947599476200837}, 'mass_row': {'max_abs': 3.469446951953614e-18, 'max_rel': 3.469446951953614e-18, 'reference_linf': 0.006615563341309325}}
```

Explicit dense wall/mass/mu mismatch:

```yaml
{'wall_rows_full': {'max_abs': 1.7763568394002505e-15, 'max_rel': 1.7154445442896305e-16, 'reference_linf': 10.355081691876224}, 'mass_row_full': {'max_abs': 0.0, 'max_rel': 0.0, 'reference_linf': 0.014812047727710342}, 'mu_column_full': {'max_abs': 0.0, 'max_rel': 0.0, 'reference_linf': 0.1100958579737913}}
```

## B-1 Gauge Fix

Method: scaled analytic gauge-complement Newton step. The linear solve acts on `row_scale * J_original * col_scale` and constrains the step to `Q_perp` of the analytic gauge generators; no residual/operator row is changed.

Gauge-sector removal evidence:

```yaml
status: PASS
gauge_rank: 511
before_gauge_sector_sigma_min: 2.5900271760812645e-13
scaled_gauge_singular_min_retained: 0.00013067938180631383
after_bordered_gauge_sector_sigma_min_not_a_lift_gate: 0.9999999999999996
after_bordered_identity_metric_is_by_construction: True
source_helpers: ["_c0c_generators_for_state(name='phase')", '_c0d_scalar_gradient_matrix', '_c0d_build_gauge_subspace', '_c0e_coupled_gauge_matrix']
```

Genuine gauge-fixed Newton proof solve:

```yaml
status: PASS
tau: 0.02925
nonconverged_seed: True
seed_description: C0f2 proof-tau state plus global phase and non-gauge r0 perturbation; not a pure global-phase orbit seed
newton_iters: 2
newton_iters_total: 2
zero_epsilon_newton_iters: 2
initial_original_residual_max: 7.163841370421038e-05
final_original_residual_linf: 3.085029523508531e-10
```

Gauge-aligned path-only proof:

```yaml
status: PASS
tau: 0.02925
tolerance: 1e-06
invariant_deltas: {'density_abs_linf': 2.471383497154589e-10, 'F_rw_abs_linf': 1.5606244915332687e-13, 'a0_gradient_abs_linf': 1.8390756148858878e-11, 'r0_abs_linf': 5.582230944156663e-08, 'mu_abs': 1.1822935519489874e-08, 'max': 5.582230944156663e-08}
original_residual: {'without_gauge_fix_linf': 6.577858049916507e-12, 'with_gauge_fix_linf': 3.085029523508531e-10, 'l2_norm_abs_delta': 2.219078686105478e-09, 'arbiter': 'stage1_solver.coupled_branch.patha_closed_branch_residual'}
```

`r0` mode preservation:

```yaml
status: PASS
mode1_lane: r0
mode1_lane_fraction: 0.9925962163519322
mode1_projection_into_gauge_basis: 1.0732491662781598e-26
mode1_sigma: 0.0008301685974640145
```

## B-2 Re-Confirm Gate

Verdict: **FOLD_DISSOLVED** via `FOLD_DISSOLVED precedence branch`.

Timeout finalization:

```yaml
{'status': 'not_used'}
```

Full-budget fold-vs-dissolve convergence ladder:

| tau | nominal_tau | relative_to_tau_fold | converged | residual_linf | newton_iters | max_newton_iters | backtrack_index | wall_seconds | measurement_status | init_source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 3.000000e-02 | above | true | 3.388320e-10 | - | 20 | 0 | 0.000000e+00 | MEASURED | c0f2_genuine_warm_start_converged_state |
| 2.950000e-02 | 2.950000e-02 | above | true | 1.108803e-12 | - | 20 | 0 | 0.000000e+00 | MEASURED | c0f2_genuine_warm_start_converged_state |
| 2.925000e-02 | 2.925000e-02 | above | true | 3.085030e-10 | 2 | 20 | 0 | 4.391008e+00 | MEASURED | c0f2_converged_state_non_gauge_perturbed_for_gauge_fixed_proof |
| 2.912500e-02 | 2.912500e-02 | above | true | 8.829016e-10 | 7 | 20 | 0 | 2.006818e+01 | MEASURED | previous_c0g_gauge_fixed_converged_state |
| 2.912400e-02 | 2.912400e-02 | above | true | 8.478202e-10 | 4 | 20 | 0 | 1.139922e+01 | MEASURED | previous_c0g_gauge_fixed_converged_state |
| 2.912250e-02 | 2.912250e-02 | below | true | 6.615192e-10 | 5 | 20 | 0 | 1.412413e+01 | MEASURED | previous_c0g_gauge_fixed_converged_state |
| 2.912200e-02 | 2.912200e-02 | below | true | 5.966944e-10 | 4 | 20 | 0 | 9.256509e+00 | MEASURED | previous_c0g_gauge_fixed_converged_state |
| 2.912000e-02 | 2.912000e-02 | below | true | 1.591087e-09 | 5 | 20 | 0 | 1.407935e+01 | MEASURED | previous_c0g_gauge_fixed_converged_state |

Full-budget dissolve contest:

```yaml
{'status': 'FAIL', 'tau_fold_reference': 0.0291233, 'below_fold_margin': 5e-07, 'above_fold_sanity_tau': 0.029124, 'above_fold_sanity_converged': True, 'below_fold_attempt_count': 3, 'below_fold_converged_count': 3, 'persistent_failure_guard': {'crawl_persistent_failure_evidence': False, 'failed_attempt_count': 0, 'attempted_backtracking': False, 'full_newton_budget': True, 'max_tau_backtracks': 5}, 'failed_below_fold_rows': []}
```

Closer sampling:

```yaml
{'status': 'MEASURED', 'tau_fold_reference': 0.0291233, 'prior_last_converged_tau': 0.029125, 'prior_last_converged_distance': 1.7000000000003124e-06, 'closest_attempt_tau': 0.029124, 'closest_attempt_distance': 6.999999999993123e-07, 'closer_than_prior_last_converged': True, 'attempts': [{'tau': 0.03, 'converged': True, 'residual_linf': 3.3883204575424486e-10, 'newton_iters': None, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'c0f2_genuine_warm_start_converged_state'}, {'tau': 0.0295, 'converged': True, 'residual_linf': 1.1088031534600695e-12, 'newton_iters': None, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'c0f2_genuine_warm_start_converged_state'}, {'tau': 0.02925, 'converged': True, 'residual_linf': 3.085029523508531e-10, 'newton_iters': 2, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'c0f2_converged_state_non_gauge_perturbed_for_gauge_fixed_proof'}, {'tau': 0.029125, 'converged': True, 'residual_linf': 8.829015643585514e-10, 'newton_iters': 7, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'previous_c0g_gauge_fixed_converged_state'}, {'tau': 0.029124, 'converged': True, 'residual_linf': 8.478202168094029e-10, 'newton_iters': 4, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'previous_c0g_gauge_fixed_converged_state'}, {'tau': 0.0291225, 'converged': True, 'residual_linf': 6.615192282843907e-10, 'newton_iters': 5, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'previous_c0g_gauge_fixed_converged_state'}, {'tau': 0.029122, 'converged': True, 'residual_linf': 5.966943883439768e-10, 'newton_iters': 4, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'previous_c0g_gauge_fixed_converged_state'}, {'tau': 0.02912, 'converged': True, 'residual_linf': 1.5910866912188792e-09, 'newton_iters': 5, 'backtrack_index': 0, 'max_newton_iters': 20, 'init_source': 'previous_c0g_gauge_fixed_converged_state'}]}
```

Cosθ support:

| tau | status | cos_theta | call | stability |
| --- | --- | --- | --- | --- |
| 3.000000e-02 | MEASURED | 5.138510e-01 | FOLD_SUPPORT | STABLE |
| 2.950000e-02 | MEASURED | 5.572308e-01 | FOLD_SUPPORT | STABLE |
| 2.925000e-02 | MEASURED | 9.068289e-01 | FOLD_SUPPORT | STABLE |
| 2.912500e-02 | MEASURED | 9.155648e-01 | FOLD_SUPPORT | STABLE |
| 2.912400e-02 | MEASURED | 9.157093e-01 | FOLD_SUPPORT | STABLE |
| 2.912250e-02 | MEASURED | 9.159419e-01 | FOLD_SUPPORT | STABLE |
| 2.912200e-02 | MEASURED | 9.160246e-01 | FOLD_SUPPORT | STABLE |
| 2.912000e-02 | MEASURED | 9.163912e-01 | FOLD_SUPPORT | STABLE |

σ_min² fit support:

```yaml
status: FAIL
linear_r2: 0.7893207736661636
linear_slope: 0.005400955978226895
tau_fold_zero_crossing: 0.029030374662435004
monotone: False
threshold_r2: 0.95
```

Bordered-conditioning trend:

```yaml
status: FAIL
cond_Jb_max_over_min: 1.0000029841608764
ratio_growth_factor: 7.254629027649889
ratio_increasing_toward_fold: False
absolute_cond_bar_used: False
```

Step-6 measured rows:

| tau | status | cond_JQ_perp | cond_Jb | sigma_min_JQ_perp | sigma_min_Jb | call |
| --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | MEASURED | 4.373305e+03 | 5.699423e+02 | 2.245655e-03 | 1.723145e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.950000e-02 | MEASURED | 6.537164e+03 | 5.699420e+02 | 1.502323e-03 | 1.723146e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.925000e-02 | MEASURED | 5.366692e+03 | 5.699406e+02 | 1.829979e-03 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912500e-02 | MEASURED | 2.111712e+04 | 5.699406e+02 | 4.650695e-04 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912400e-02 | MEASURED | 2.240638e+04 | 5.699406e+02 | 4.383096e-04 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912250e-02 | MEASURED | 2.487272e+04 | 5.699406e+02 | 3.948475e-04 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912200e-02 | MEASURED | 2.589380e+04 | 5.699406e+02 | 3.792773e-04 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912000e-02 | MEASURED | 3.172661e+04 | 5.699406e+02 | 3.095487e-04 | 1.723150e-02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |

FOLD_DISSOLVED precedence check:

```yaml
{'status': 'PASS', 'below_fold_margin': 5e-07, 'required_count': 2, 'below_fold_converged_count': 3, 'branch_continuous_count': 3, 'no_r0_turning_signature': True, 'rows': [{'tau': 0.0291225, 'residual_linf': 6.615192282843907e-10, 'newton_iters': 5, 'backtrack_index': 0, 'init_source': 'previous_c0g_gauge_fixed_converged_state', 'used_existing_b2c': False}, {'tau': 0.029122, 'residual_linf': 5.966943883439768e-10, 'newton_iters': 4, 'backtrack_index': 0, 'init_source': 'previous_c0g_gauge_fixed_converged_state', 'used_existing_b2c': False}, {'tau': 0.02912, 'residual_linf': 1.5910866912188792e-09, 'newton_iters': 5, 'backtrack_index': 0, 'init_source': 'previous_c0g_gauge_fixed_converged_state', 'used_existing_b2c': False}]}
```

## Scope Guard

```yaml
staged_items_built: ['B-4', 'B-1', 'B-2']
B3_pseudoarclength_built: False
B4_analytic_sparse_assembly_built: True
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
patha_closed_branch_residual_touched: False
faithful_operators_touched: False
xi_or_grad_div_penalty_touched: False
physical_export_permitted_touched: False
prefer_existing_b2c_background_predictor: False
absolute_cond_gt_1e10_bar_used_for_B2: False
```

Machine JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_build_B1B2/pathA_C0g_build_B1B2.json`
