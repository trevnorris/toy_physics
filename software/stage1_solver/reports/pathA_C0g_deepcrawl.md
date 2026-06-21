# Path-A C0g Deeper Crawl

Work-1 verdict: **STILL_GOING**

## Scope

- Continuation is gauge-fixed and path-only: `row_scale * J_original * col_scale`, step constrained to `Q_perp`.
- The original `patha_closed_branch_residual` is the sole convergence arbiter at `Linf <= 1e-6`.
- B-3 pseudo-arclength was not built.
- Frozen residual/operators/`(1/xi)` penalty/export/physics are not edited by this runner.

## B-4 Reconfirm

```yaml
status: PASS
reference_tau: 0.029125
assembly_speedup_x: 80.60344287391108
tolerances: {'max_abs': 1e-10, 'max_rel': 1e-08}
```

## Deeper Ladder

| tau | converged | residual_linf | newton_iters | backtracks | wall_seconds | min_R0 | min_rho | mu | sigma_min_JQ_perp | admissibility |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.912250e-02 | true | 6.615192e-10 | 5 | 0 | 1.412413e+01 | 7.938496e-01 | 7.215848e-06 | 2.035272e+00 | 3.948475e-04 | PASS |
| 2.912200e-02 | true | 5.966944e-10 | 4 | 0 | 9.256509e+00 | 7.940100e-01 | 7.216663e-06 | 2.035241e+00 | 3.792773e-04 | PASS |
| 2.912000e-02 | true | 1.591087e-09 | 5 | 0 | 1.407935e+01 | 7.947284e-01 | 7.220304e-06 | 2.035100e+00 | 3.095487e-04 | PASS |
| 2.911813e-02 | true | 8.858940e-10 | 7 | 6 | 1.610475e+01 | 7.955946e-01 | 7.224675e-06 | 2.034931e+00 | 2.254878e-04 | PASS |
| 2.911625e-02 | true | 2.223002e-07 | 20 | 5 | 3.996254e+01 | 7.964222e-01 | 7.228790e-06 | 2.034770e+00 | 1.455266e-04 | PASS |
| 2.911443e-02 | true | 7.526114e-07 | 2 | 6 | 4.914247e+00 | 7.967120e-01 | 7.230142e-06 | 2.034714e+00 | 1.182048e-04 | PASS |
| 2.911380e-02 | true | 8.603467e-07 | 3 | 0 | 6.664580e+00 | 7.971193e-01 | 7.232158e-06 | 2.034636e+00 | 7.885724e-05 | PASS |
| 2.911320e-02 | true | 9.798546e-07 | 20 | 0 | 3.018961e+01 | 7.978618e-01 | 7.235834e-06 | 2.034492e+00 | 7.026590e-06 | PASS |
| 2.911310e-02 | false | 1.040272e-06 | 20 | 2 | 3.016354e+01 | 7.983671e-01 | 7.238331e-06 | 2.034394e+00 | - | - |
| 2.911300e-02 | false | 1.095938e-06 | 20 | 1 | 3.186863e+01 | 7.985136e-01 | 7.239047e-06 | 2.034366e+00 | - | - |
| 2.911280e-02 | false | 1.227104e-06 | 20 | 0 | 3.239970e+01 | 7.988153e-01 | 7.240515e-06 | 2.034308e+00 | - | - |
| 2.911265e-02 | false | 1.208956e-06 | 20 | 6 | 3.264085e+01 | 7.974310e-01 | 7.233648e-06 | 2.034576e+00 | - | - |
| 2.911262e-02 | false | 1.321974e-06 | 3 | 5 | 6.661046e+00 | 7.969387e-01 | 7.231171e-06 | 2.034672e+00 | - | - |
| 2.911250e-02 | false | 1.235552e-06 | 2 | 4 | 5.281659e+00 | 7.976577e-01 | 7.234772e-06 | 2.034532e+00 | - | - |
| 2.911086e-02 | false | 1.800293e-06 | 20 | 5 | 3.301523e+01 | 7.979366e-01 | 7.236055e-06 | 2.034479e+00 | - | - |
| 2.910898e-02 | false | 2.497993e-06 | 20 | 4 | 3.182716e+01 | 7.974098e-01 | 7.233305e-06 | 2.034582e+00 | - | - |
| 2.910728e-02 | false | 3.215789e-06 | 20 | 4 | 4.206734e+01 | 7.990428e-01 | 7.241277e-06 | 2.034268e+00 | - | - |
| 2.910500e-02 | false | 3.950694e-06 | 20 | 3 | 3.513304e+01 | 7.988123e-01 | 7.240003e-06 | 2.034314e+00 | - | - |
| 2.910172e-02 | false | 5.031264e-06 | 20 | 3 | 3.207732e+01 | 7.984381e-01 | 7.237954e-06 | 2.034388e+00 | - | - |
| 2.910013e-02 | false | 5.620580e-06 | 20 | 3 | 3.835454e+01 | 7.973375e-01 | 7.232369e-06 | 2.034603e+00 | - | - |
| 2.909000e-02 | false | 9.125243e-06 | 20 | 2 | 3.284333e+01 | 7.983822e-01 | 7.236928e-06 | 2.034407e+00 | - | - |
| 2.908719e-02 | false | 1.021696e-05 | 20 | 2 | 3.767882e+01 | 7.971231e-01 | 7.230449e-06 | 2.034653e+00 | - | - |
| 2.908583e-02 | false | 1.067776e-05 | 20 | 2 | 3.774863e+01 | 7.989471e-01 | 7.239443e-06 | 2.034301e+00 | - | - |
| 2.906000e-02 | false | 1.971160e-05 | 20 | 1 | 3.852653e+01 | 7.973614e-01 | 7.229887e-06 | 2.034626e+00 | - | - |
| 2.905813e-02 | false | 2.035831e-05 | 20 | 1 | 3.834057e+01 | 7.990064e-01 | 7.237968e-06 | 2.034309e+00 | - | - |
| 2.905722e-02 | false | 2.059922e-05 | 7 | 1 | 1.611831e+01 | 7.981421e-01 | 7.233647e-06 | 2.034476e+00 | - | - |
| 2.900000e-02 | false | 4.064629e-05 | 7 | 0 | 1.779291e+01 | 7.990743e-01 | 7.234589e-06 | 2.034335e+00 | - | - |
| 2.900000e-02 | false | 4.073542e-05 | 8 | 0 | 2.065262e+01 | 7.975072e-01 | 7.226714e-06 | 2.034639e+00 | - | - |
| 2.900000e-02 | false | 4.072752e-05 | 5 | 0 | 1.301051e+01 | 7.975408e-01 | 7.226896e-06 | 2.034632e+00 | - | - |

## Sigma Trend

```yaml
status: MEASURED
call: LINEAR_MONOTONE_FOLD_SUPPORT
linear_r2: 0.9566269269579496
linear_slope: 0.016699204434613493
tau_fold_zero_crossing: 0.029113886909519333
monotone_toward_stall: True
branch_tangent: NO_TANGENT_REVERSAL
```

## Verdict Support

```yaml
{'verdict': 'STILL_GOING', 'why': 'deeper crawl did not meet the strict no-fold or real-fold criteria', 'tau_fold_reference': 0.0291233, 'deepest_converged_tau': 0.0291132, 'deepest_margin_past_tau_fold': 1.0100000000002468e-05, 'new_deep_converged_count': 5, 'required_new_deep_converged_count': 5, 'all_converged_below_admissible': True, 'sigma_min_finite': True, 'sigma_min_finite_floor': 1e-08, 'sigma_min_min': 7.026590471648914e-06, 'sigma_min_collapse_factor': 56.193325904739346, 'sigma_min_collapsing_toward_zero': True, 'sigma_fit_supports_fold': True, 'branch_tangent_reversal': False, 'r0_trend': {'status': 'FAIL', 'direction': 'nonmonotone', 'first': 0.7938496490533923, 'last': 0.7975408112469868, 'delta': 0.003691162193594555, 'tolerance': 1e-08}, 'mu_trend': {'status': 'FAIL', 'direction': 'nonmonotone', 'first': 2.0352723136396618, 'last': 2.0346319798279167, 'delta': -0.0006403338117451085, 'tolerance': 2.035272313639662e-08}, 'persistent_failure_guard': {'crawl_persistent_failure_evidence': True, 'failed_attempt_count': 21, 'attempted_backtracking': True, 'full_newton_budget': True, 'max_tau_backtracks': 6}, 'failed_attempts': [{'tau': 0.029, 'nominal_tau': 0.029, 'residual_linf': 4.064628690073069e-05, 'backtrack_index': 0, 'newton_iters': 7, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029060000000000002, 'nominal_tau': 0.029, 'residual_linf': 1.9711598458833218e-05, 'backtrack_index': 1, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.02909, 'nominal_tau': 0.029, 'residual_linf': 9.125243420286117e-06, 'backtrack_index': 2, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029105, 'nominal_tau': 0.029, 'residual_linf': 3.95069446401837e-06, 'backtrack_index': 3, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.0291125, 'nominal_tau': 0.029, 'residual_linf': 1.2355521524336222e-06, 'backtrack_index': 4, 'newton_iters': 2, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029, 'nominal_tau': 0.029, 'residual_linf': 4.073542111516712e-05, 'backtrack_index': 0, 'newton_iters': 8, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029058125, 'nominal_tau': 0.029, 'residual_linf': 2.0358314906841757e-05, 'backtrack_index': 1, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.0290871875, 'nominal_tau': 0.029, 'residual_linf': 1.0216961690595605e-05, 'backtrack_index': 2, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029101718749999998, 'nominal_tau': 0.029, 'residual_linf': 5.031264198427851e-06, 'backtrack_index': 3, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029108984375, 'nominal_tau': 0.029, 'residual_linf': 2.4979933947785055e-06, 'backtrack_index': 4, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.0291126171875, 'nominal_tau': 0.029, 'residual_linf': 1.3219738242923468e-06, 'backtrack_index': 5, 'newton_iters': 3, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029, 'nominal_tau': 0.029, 'residual_linf': 4.072752199659781e-05, 'backtrack_index': 0, 'newton_iters': 5, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029057216796875, 'nominal_tau': 0.029, 'residual_linf': 2.0599221343283597e-05, 'backtrack_index': 1, 'newton_iters': 7, 'message': 'line search failed to reduce original residual'}, {'tau': 0.029085825195312498, 'nominal_tau': 0.029, 'residual_linf': 1.0677758429314746e-05, 'backtrack_index': 2, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.02910012939453125, 'nominal_tau': 0.029, 'residual_linf': 5.62058031994099e-06, 'backtrack_index': 3, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029107281494140623, 'nominal_tau': 0.029, 'residual_linf': 3.21578907698343e-06, 'backtrack_index': 4, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029110857543945313, 'nominal_tau': 0.029, 'residual_linf': 1.8002932197086485e-06, 'backtrack_index': 5, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029112645568847656, 'nominal_tau': 0.029, 'residual_linf': 1.208956076467066e-06, 'backtrack_index': 6, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.0291128, 'nominal_tau': 0.0291128, 'residual_linf': 1.2271042550866806e-06, 'backtrack_index': 0, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.029113, 'nominal_tau': 0.0291128, 'residual_linf': 1.095938063055621e-06, 'backtrack_index': 1, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}, {'tau': 0.0291131, 'nominal_tau': 0.0291128, 'residual_linf': 1.040271644558305e-06, 'backtrack_index': 2, 'newton_iters': 20, 'message': 'maximum Newton iterations reached'}]}
```

## Admissibility

Overall: **PASS**

| tau | status | positive_density_status | min_rho | sane_min_R0_status | min_R0 | residual_linf | epsilon_independence_status | zero_aid_residual_equivalence_max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.912250e-02 | PASS | PASS | 7.215848e-06 | PASS | 7.938496e-01 | 6.615192e-10 | PASS | 0.000000e+00 |
| 2.912200e-02 | PASS | PASS | 7.216663e-06 | PASS | 7.940100e-01 | 5.966944e-10 | PASS | 0.000000e+00 |
| 2.912000e-02 | PASS | PASS | 7.220304e-06 | PASS | 7.947284e-01 | 1.591087e-09 | PASS | 0.000000e+00 |
| 2.911625e-02 | PASS | PASS | 7.228790e-06 | PASS | 7.964222e-01 | 2.223002e-07 | PASS | 0.000000e+00 |
| 2.911813e-02 | PASS | PASS | 7.224675e-06 | PASS | 7.955946e-01 | 8.858940e-10 | PASS | 0.000000e+00 |
| 2.911443e-02 | PASS | PASS | 7.230142e-06 | PASS | 7.967120e-01 | 7.526114e-07 | PASS | 0.000000e+00 |
| 2.911380e-02 | PASS | PASS | 7.232158e-06 | PASS | 7.971193e-01 | 8.603467e-07 | PASS | 0.000000e+00 |
| 2.911320e-02 | PASS | PASS | 7.235834e-06 | PASS | 7.978618e-01 | 9.798546e-07 | PASS | 0.000000e+00 |

## Empty-Core Target

```yaml
{'status': 'MEASURED', 'first_tau': 0.0291225, 'deepest_tau': 0.0291132, 'first_min_R0': 0.7938496490533923, 'deepest_min_R0': 0.7978618226749459, 'first_min_rho': 7.2158480974477856e-06, 'deepest_min_rho': 7.2358341366518185e-06, 'approached_small_core_regime': False, 'deepest_solvable': True, 'interpretation': 'min_R0/min_rho did not move toward the small-core regime over this ladder'}
```

## Scope Guard

```yaml
B3_pseudoarclength_built: False
B4_analytic_sparse_assembly_built: True
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
patha_closed_branch_residual_touched: False
faithful_operators_touched: False
xi_or_grad_div_penalty_touched: False
physical_export_permitted_touched: False
prefer_existing_b2c_background_predictor: False
frozen_physics_touched: False
```

Machine JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_deepcrawl/pathA_C0g_deepcrawl.json`
