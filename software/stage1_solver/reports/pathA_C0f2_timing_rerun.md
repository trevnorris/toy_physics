# Path-A C0f2 Timing Rerun

Phase status: **crawl_stalled_backtrack_exhausted**
Slow-vs-fold call: **SLOW_BUT_DESCENT_DIRECTION_EXISTS**

## Per-Tau Timing Table

| tau | wall_s | iters | bt | Linf | L2 | J_asm_s | lin_solve_s | resid_s | JVP_asm | JVP_lin | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 1.701442e+02 | 0 | 0 | 3.388320e-10 | 1.469884e-09 | 1.645622e+02 | 5.497807e+00 | 7.570402e-02 | 1012 | 34 | converged_initial |
| 2.950000e-02 | 1.230486e+02 | 0 | 0 | 1.108803e-12 | 1.936762e-12 | 1.194701e+02 | 3.527486e+00 | 4.654647e-02 | 759 | 23 | converged_initial |
| 2.925000e-02 | 1.259141e+02 | 0 | 0 | 6.577858e-12 | 2.107759e-11 | 1.217491e+02 | 4.113256e+00 | 4.731946e-02 | 759 | 26 | converged_initial |
| 2.912500e-02 | 1.672723e+02 | 0 | 0 | 3.393808e-11 | 1.191665e-10 | 1.618822e+02 | 5.326958e+00 | 5.676752e-02 | 1012 | 34 | converged_initial |
| 2.906250e-02 | 2.969393e+02 | 6 | 29 | 2.423837e-05 | 8.719489e-05 | 2.861694e+02 | 1.035383e+01 | 4.048594e-01 | 1771 | 65 | line_search_exhausted |
| 2.909375e-02 | 4.250044e+02 | 9 | 56 | 9.988005e-06 | 3.645407e-05 | 4.084976e+02 | 1.589598e+01 | 5.935040e-01 | 2530 | 100 | line_search_exhausted |
| 2.911250e-02 | 3.705365e+02 | 8 | 63 | 1.561491e-06 | 5.790302e-06 | 3.558303e+02 | 1.408762e+01 | 6.015772e-01 | 2277 | 92 | line_search_exhausted |

## Depth And Extrapolation

```yaml
deepest_linf_pass_tau: 0.029125
deepest_accepted_tau: 0.029125
reached_tau_0p028: False
extrapolation_status: MEASURED
extrapolation_assumed_target_tau: 0.02
extrapolation_model: recent_median_per_fine_tau_step
extrapolation_model_detail: Uses the median wall-clock over the deepest measured converged steps and the median accepted tau spacing; no claim of asymptotic validity.
extrapolation_deepest_measured_tau: 0.029125
extrapolation_median_recent_step_seconds: 146.5932092849398
extrapolation_median_tau_step: 0.00024999999999999675
extrapolation_remaining_steps: 37
extrapolation_measured_elapsed_seconds: 586.3791984578129
extrapolation_remaining_seconds: 5423.948743542773
extrapolation_estimated_total_seconds_from_tau_0p03: 6010.327942000586
extrapolation_blocked_by_measured_stall: True
extrapolation_validity: conditional_trend_only_not_a_completion_prediction
```

## Slow Vs Fold Evidence

```yaml
call: SLOW_BUT_DESCENT_DIRECTION_EXISTS
reason: merit sweep found true residual decrease for at least one alpha
linear_rel_resid: 1.0509779632646936e-12
best_alpha: 0.0078125
best_actual_l2: 8.715214325356727e-05
fold_call: FOLD_RISK
fold_growth_factor: 3.2017184946023476
fold:
  status: MEASURED
  call: FOLD_RISK
  reason: growth_condition_false
  growth_factor: 3.2017184946023476
  growth_threshold: 10.0
  growth_condition: False
merit:
  status: MEASURED
  reason: None
  tau: 0.0290625
  initial_l2: 8.719489474327922e-05
  initial_linf: 2.423836883522386e-05
  linear_rel_resid: 1.0509779632646936e-12
  any_alpha_reduces_true_l2: True
  best_alpha: 0.0078125
  best_actual_l2: 8.715214325356727e-05
  diagnostic_interpretation: LOCALIZED_GLOBALIZATION_EVIDENCE
```

| alpha | actual_l2 | pred_l2 | actual_linf | pred_linf | reduces |
| --- | --- | --- | --- | --- | --- |
| 1.000000e+00 | 3.734802e-04 | 9.163991e-17 | 1.373951e-04 | 4.071179e-17 | false |
| 5.000000e-01 | 1.817822e-04 | 4.359745e-05 | 6.703762e-05 | 1.211918e-05 | false |
| 2.500000e-01 | 1.155941e-04 | 6.539617e-05 | 4.215896e-05 | 1.817878e-05 | false |
| 1.250000e-01 | 9.442557e-05 | 7.629553e-05 | 3.234524e-05 | 2.120857e-05 | false |
| 6.250000e-02 | 8.875218e-05 | 8.174521e-05 | 2.808854e-05 | 2.272347e-05 | false |
| 3.125000e-02 | 8.743364e-05 | 8.447005e-05 | 2.611287e-05 | 2.348092e-05 | false |
| 1.562500e-02 | 8.717778e-05 | 8.583247e-05 | 2.516300e-05 | 2.385964e-05 | true |
| 7.812500e-03 | 8.715214e-05 | 8.651368e-05 | 2.469753e-05 | 2.404901e-05 | true |
| 3.906250e-03 | 8.716497e-05 | 8.685429e-05 | 2.446716e-05 | 2.414369e-05 | true |
| 1.953125e-03 | 8.717780e-05 | 8.702459e-05 | 2.435257e-05 | 2.419103e-05 | true |
| 9.765625e-04 | 8.718581e-05 | 8.710974e-05 | 2.429542e-05 | 2.421470e-05 | true |
| 4.882812e-04 | 8.719022e-05 | 8.715232e-05 | 2.426688e-05 | 2.422653e-05 | true |
| 2.441406e-04 | 8.719252e-05 | 8.717361e-05 | 2.425262e-05 | 2.423245e-05 | true |
| 1.220703e-04 | 8.719370e-05 | 8.718425e-05 | 2.424549e-05 | 2.423541e-05 | true |
| 6.103516e-05 | 8.719430e-05 | 8.718957e-05 | 2.424193e-05 | 2.423689e-05 | true |
| 3.051758e-05 | 8.719459e-05 | 8.719223e-05 | 2.424015e-05 | 2.423763e-05 | true |
| 1.525879e-05 | 8.719474e-05 | 8.719356e-05 | 2.423926e-05 | 2.423800e-05 | true |
| 7.629395e-06 | 8.719482e-05 | 8.719423e-05 | 2.423881e-05 | 2.423818e-05 | true |
| 3.814697e-06 | 8.719486e-05 | 8.719456e-05 | 2.423859e-05 | 2.423828e-05 | true |
| 1.907349e-06 | 8.719488e-05 | 8.719473e-05 | 2.423848e-05 | 2.423832e-05 | true |
| 9.536743e-07 | 8.719489e-05 | 8.719481e-05 | 2.423842e-05 | 2.423835e-05 | true |

## Genuineness And Scope

```yaml
genuine_continuation: True
prefer_existing_b2c_background_predictor: False
used_existing_b2c_permitted: False
warm_start_rule: previous_c0_converged_state
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
max_newton_iters: 20
max_tau_backtracks: 5
line_search: armijo
max_line_search_iters: 20
operators_touched_by_c0f2: False
frozen_physics_touched_by_c0f2: False
physical_export_guard_touched_by_c0f2: False
solver_logic_changed_by_c0f2: False
trust_region_dogleg_lm_implemented: False
pseudo_arclength_implemented: False
gauge_deflation_implemented: False
depth_continuation: tau_only
```

Machine JSON: `runs/pathA_C0f2_timing_rerun/pathA_C0f2_timing_rerun.json`
Checkpoint JSON: `runs/pathA_C0f2_timing_rerun/pathA_C0f2_timing_checkpoint.json`

