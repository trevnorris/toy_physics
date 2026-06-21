# Path-A C0f Globalization Probe

C0f-3 verdict: **DIAGNOSTIC_INCOMPLETE**

## Default Config

```yaml
grid: [16, 16]
depth_sequence: [0.03, 0.0295, 0.02925, 0.029125, 0.0290625, 0.029, 0.0285, 0.028]
max_newton_iters: 20
max_newton_iters_override: None
max_tau_backtracks: 5
line_search: armijo
max_line_search_iters: 20
line_search_shrink: 0.5
epsilon_schedule: [{'core_epsilon': 0.08, 'k1_radius_epsilon': 0.8, 'preconditioner_diagonal_shift': 1e-12, 'scale_floor': 1e-14, 'scale_max': 100000000.0, 'scale_min': 1e-08, 'use_jacobi_scaling': True}, {'core_epsilon': 0.02, 'k1_radius_epsilon': 0.4, 'preconditioner_diagonal_shift': 1e-12, 'scale_floor': 1e-14, 'scale_max': 100000000.0, 'scale_min': 1e-08, 'use_jacobi_scaling': True}, {'core_epsilon': 0.0, 'k1_radius_epsilon': 0.0, 'preconditioner_diagonal_shift': 1e-12, 'scale_floor': 1e-14, 'scale_max': 100000000.0, 'scale_min': 1e-08, 'use_jacobi_scaling': True}]
use_wall_predictor: True
eos_final_only_after_first: True
continuation_K_values_at_tau_0p03: [0.08, 0.18]
path_only_aid_confirmation:
  eos_continuation_scope: path continuation; final accepted solve uses final EOS K
  epsilon_aids_scope: preconditioner residual path only
  final_accepted_solve: zero core_epsilon and zero k1_radius_epsilon
  jacobi_scaling_scope: linear system row/column scaling only; recovered step checked by original residual
  jvp_residual: original physical residual
  newton_residual_merit_line_search_convergence: original physical residual
  physical_residual_entry_point: stage1_solver.coupled_branch.patha_closed_branch_residual
  preconditioner_residual_entry_point: stage1_solver.patha_c0_conditioning_spike.c0_preconditioner_residual
  wall_predictor_scope: initial guess for tau continuation only
```

## C0f-0 Default Crawl

| tau | solver_converged | b2c_linf_pass | Linf | L2 | iters | backtracks | alpha | tau_bt |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | true | true | 3.388320e-10 | 1.469884e-09 | 0 | 0 | - | 0 |
| 2.950000e-02 | true | true | 1.108803e-12 | 1.936762e-12 | 0 | 0 | - | 0 |
| 2.925000e-02 | true | true | 6.577858e-12 | 2.107759e-11 | 0 | 0 | - | 0 |
| 2.912500e-02 | true | true | 3.393808e-11 | 1.191665e-10 | 0 | 0 | - | 0 |
| 2.906250e-02 | false | false | nan | nan | - | 0 | - | 0 |

## C0f-0a Numeric Fold Detector

```yaml
status: MEASURED
call: FOLD_RISK
reason: growth_condition_false
growth_factor: 3.2017184946023476
growth_threshold: 10.0
monotone_growth_window: True
growth_condition: False
smaller_tau_failure_condition: {'attempted_target_taus': [0.0290625], 'candidate_nominal_tau': 0.0290625, 'failed_all_backtracks': False, 'max_backtrack_index': 0, 'max_tau_backtracks': 5, 'min_abs_delta_tau': 6.250000000000006e-05, 'min_tau_backtrack_delta': 5e-05, 'reason': 'smaller_tau_failed_but_backtrack_exhaustion_not_observed'}
```
| from | to | full | psi | R0 | mu | A |
| --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 2.950000e-02 | 2.352990e+01 | 3.830275e+00 | 2.841289e+01 | 3.935157e+00 | 3.421919e+00 |
| 2.950000e-02 | 2.925000e-02 | 3.629589e+01 | 5.365849e+00 | 4.371010e+01 | 5.690603e+00 | 4.728489e+00 |
| 2.925000e-02 | 2.912500e-02 | 7.533612e+01 | 1.024663e+01 | 9.045974e+01 | 1.117796e+01 | 8.914167e+00 |

## C0f-1 Merit Sweep

```yaml
status: NOT_MEASURED
reason: first_stalled_or_incomplete_tau_has_no_saved_state
tau: 0.0290625
initial_l2: None
initial_linf: None
linear_rel_resid: None
any_alpha_reduces_true_l2: None
smallest_reducing_alpha: None
diagnostic_interpretation: None
```

## C0f-2 1-P_G Gauge Re-Confirm

```yaml
status: MEASURED
reason: None
tau: 0.02899
classifier: dimensionless 1-P_G / 1-P_cpl; raw curl reported only
thresholds: {'gauge_fraction_max': 0.02, 'transverse_candidate_fraction_min': 0.05}
old_c0e_labels_ignored: True
raw_curl_used_as_classifier: False
state_artifact: /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/states/attempt_001_tau_0p02899.npz
matrix_path: /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/matrices/attempt_tau_0p02899_bt_0.npz
phase_mode: {'classification': 'PHASE_GAUGE_BY_PHASE_AND_COUPLED_CAPTURE', 'mode_index': 0, 'one_minus_p_cpl': 1.3485569312798873e-05, 'one_minus_p_g': None, 'one_minus_phase_capture': 0.00581814140092074, 'raw_curl_reported_only': {'cell_volume_weighted_norm': 0.005948691935595193, 'unweighted_norm': 0.012603332949842101}, 'sigma': 1.43587972095615e-12}
```
| mode | sigma | 1-P_G | 1-P_G_A | 1-P_cpl | boundary_curl | class |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.720320e-11 | 6.117609e-03 | 7.101717e-04 | 8.404271e-04 | 6.057128e-01 | GAUGE_BY_1_MINUS_P |
| 2 | 2.569643e-10 | 1.700267e-01 | 1.698799e-01 | 1.698549e-01 | 6.054464e-01 | TRANSVERSE_CANDIDATE_BY_1_MINUS_P |
| 3 | 5.917734e-10 | 4.445522e-04 | 2.302309e-04 | 2.308251e-04 | 6.054124e-01 | GAUGE_BY_1_MINUS_P |
| 4 | 2.808634e-08 | 1.078823e-03 | 1.063877e-03 | 1.068411e-03 | 6.060483e-01 | GAUGE_BY_1_MINUS_P |

Mixed controls (exact gradient + epsilon transverse):

| k | eps | expected | 1-P_G | 1-P_cpl | curl_w | boundary_curl |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.000000e-01 | 9.900990e-03 | 9.900990e-03 | 9.900990e-03 | 1.232799e-01 | 4.891918e-01 |
| 2 | 1.000000e-01 | 9.900990e-03 | 9.900990e-03 | 9.900990e-03 | 1.331442e-01 | 5.308010e-01 |
| 4 | 1.000000e-01 | 9.900990e-03 | 9.900990e-03 | 9.900990e-03 | 1.696431e-01 | 6.206627e-01 |
| 7 | 1.000000e-01 | 9.900990e-03 | 9.900990e-03 | 9.900990e-03 | 2.224611e-01 | 7.423880e-01 |

## Verdict Support

```yaml
accepted_default_crawl_through_target: False
deepest_accepted_tau: 0.029125
deepest_linf_pass_tau: 0.029125
fold_call: FOLD_RISK
fold_status: MEASURED
merit_status: NOT_MEASURED
reason: required fold or merit evidence not measured for a substantive verdict
recommended_next_step: Rerun the missing bounded diagnostic evidence before selecting C0g.
```

## Scope Guard

```yaml
depth_continuation: tau_only
existing_closed_newton_crawl_used: True
frozen_physics_touched_by_c0f: False
gauge_deflation_implemented: False
operators_touched_by_c0f: False
physical_export_guard_touched_by_c0f: False
pseudo_arclength_implemented: False
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
solver_logic_changed_by_c0f: False
trust_region_dogleg_lm_implemented: False
```

## Verification

```yaml
chunk_gates: {'status': 'PASS', 'command': 'PYTHONPATH=src timeout 600 pytest tests/test_patha_static_balance.py tests/test_patha_closed_newton.py tests/test_patha_closed_validation.py -q', 'summary': '11 passed, 18 warnings in 18.16s'}
focused_c0f_tests: {'status': 'PASS', 'command': 'PYTHONPATH=src pytest tests/test_patha_c0f_globalization_probe.py -q', 'summary': '4 passed'}
full_python_suite: {'status': 'PASS', 'command': 'PYTHONPATH=software/stage1_solver/src timeout 600 pytest software/stage1_solver/tests -q', 'workdir': '/var/projects/toy_physics', 'summary': '111 passed, 19 warnings in 165.54s'}
```

Machine artifact: `runs/pathA_C0f_globalization_probe/pathA_C0f_globalization_probe.json`

