# Path-A C0g B-3 Run Pseudo-Arclength

> **ORCHESTRATOR CORRECTION (post dual-review, 2026-06-21).** This report's Outcome-A criteria were met without gaming
> (FIDELITY-CLEAN), but the adversarial review rated it **EARNED-BUT-OVERSTATED**, with two corrections that override the prose below:
> (1) **Likely turned EARLY:** the continuation turned at τ=0.0291159 where σ_min(J·Q_perp)≈1.0e-5, but the fixed-τ deep-point run
> reached σ_min≈4.8e-7 at τ=0.0291150 (DEEPER + ~20× more singular). A genuine fold turns AT the σ_min minimum; this turned shallow
> of it ⇒ "rounded the fold" vs "wobbled shallow of the real singularity" is UNRESOLVED (needs a smaller-ds re-run).
> (2) **No deeper throat — "Fresh Deep Tight State" is MISLABELED:** min_R0 RISES 0.7975→0.7992 across the trace (throat OPENS) and
> min_rho is flat ~7.2e-6; the τ_fold state (residual 6.7e-8) is branch-accepted (≤1e-6), NOT "tight" (≤1e-11), and NOT deeper than
> fixed-τ Newton already reached. Rounding this fold delivered NO deeper/empty-core throat. min_R0 is stuck ~0.798 across the whole
> τ range ⇒ τ does not control throat depth. See `decisions/13` §0.

Outcome: **A - GENUINE_FOLD (numerically real but EARNED-BUT-OVERSTATED — see correction above)**

## Frozen Schedule

```yaml
initial_ds: 0.001
min_ds: 4e-05
ds_shrink: 0.5
max_retries_per_step: 30
tau_scale: 1000.0
residual_tolerance: 1e-06
hard_wall_no_progress_N: 30
hard_wall_epsilon_res: 0.01
hard_wall_epsilon_adv: 1e-09
reseed_order: tangent -> tangent_reseed -> secant_fallback -> step_shrink
single_arbiter: original patha_closed_branch_residual Linf <= BACKGROUND_RESIDUAL_TOL
```

## Seed

```yaml
anchor_tau: 0.02911625
anchor_state_path: /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_followup/C1/states/continuation_B3_0_accepted_009_tau_0p02911625.npz
anchor_provenance: fresh_tight_reconverge_from_B3_followup_C1
fresh_accepted_state: false
```

## Accepted Trace

| accepted_index | s | solved_tau | dtauds | finite_difference_delta_tau_over_delta_s | original_residual_linf | sigma_min_JQ_perp | min_R0 | min_rho | predictor_step | predictor_source | fallback_used | secant_fallback_available | wall_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.000000e-03 | 2.911603e-02 | -2.217673e-04 | -2.217673e-04 | 2.165928e-08 | 4.435093e-05 | 7.974621e-01 | 7.234027e-06 | 1.000000e-03 | tangent | false | false | 4.193578e+00 |
| 1 | 2.000000e-03 | 2.911590e-02 | -1.283556e-04 | -1.283556e-04 | 4.405426e-08 | 1.009752e-05 | 7.978157e-01 | 7.235788e-06 | 1.000000e-03 | tangent | false | true | 3.460871e+00 |
| 2 | 3.000000e-03 | 2.911587e-02 | -2.949540e-05 | -2.949540e-05 | 6.679543e-08 | 2.441494e-05 | 7.981719e-01 | 7.237559e-06 | 1.000000e-03 | tangent | false | true | 3.511525e+00 |
| 3 | 4.000000e-03 | 2.911594e-02 | 7.126806e-05 | 7.126806e-05 | 8.943214e-08 | 5.884328e-05 | 7.985273e-01 | 7.239322e-06 | 1.000000e-03 | tangent | false | true | 3.601710e+00 |
| 4 | 5.000000e-03 | 2.911611e-02 | 1.699403e-04 | 1.699403e-04 | 1.115194e-07 | 9.284644e-05 | 7.988784e-01 | 7.241060e-06 | 1.000000e-03 | tangent | false | true | 3.628581e+00 |
| 5 | 6.000000e-03 | 2.911637e-02 | 2.628994e-04 | 2.628994e-04 | 1.326845e-07 | 1.261275e-04 | 7.992220e-01 | 7.242758e-06 | 1.000000e-03 | tangent | false | true | 3.673803e+00 |

## Rejected Attempts

| step_index | retry_index | s | tau | dtauds | original_residual_linf | sigma_min_JQ_perp | ds | predictor_source | fallback_used | fallback_failed_reason | fallback_not_applicable_reason | reason | ds_schedule_exhausted | no_forward_progress_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Decisive Evidence

```yaml
{'outcome': 'A', 'label': 'GENUINE_FOLD', 'why': 'accepted tau(s) has a local minimum with original-residual accepted far-side points and sigma recovery', 'tau_fold': 0.029115870381663817, 'turning_accepted_index': 2, 'accepted_point_geometry': {'status': 'PASS', 'decreasing_points_before': 2, 'increasing_points_after': 3, 'taus': [0.02911602823266734, 0.029115899877063114, 0.029115870381663817, 0.02911594164971883, 0.029116111590049242, 0.029116374489480627]}, 'single_arbiter': {'status': 'PASS', 'residual_tolerance': 1e-06, 'max_original_residual_linf': 1.326845179960151e-07}, 'sigma_recovery': {'status': 'PASS', 'near_minimum_sigma_min_JQ_perp': 1.0097524498429972e-05, 'far_side_sigma_min_JQ_perp': [5.884327912801298e-05, 9.284644147883705e-05, 0.00012612747102714705], 'monotone_increase_over_required_window': True, 'recovery_factor': 12.490929935030923, 'required_factor': 10.0, 'required_monotone_points': 3}}
```

## Fresh Deep Tight State

```yaml
{'tau': 0.029115870381663817, 's': 0.0030000000000000603, 'state_artifact': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/B3_run/states/accepted_002_tau_0p0291158703817.npz', 'original_residual_linf': 6.679542915132008e-08, 'min_R0': 0.7981719169001027, 'min_rho': 7.237559050747055e-06, 'provenance': 'corrector_accepted_fresh'}
```

## Deeper-State Mode Comparison

{'status': 'MEASURED', 'comparison': 'declared_anchor_vs_fresh_deepest_corrector_accepted_state', 'self_comparison': False, 'anchor_tau': 0.02911625, 'fresh_tau': 0.029115870381663817, 'fresh_state_artifact': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/B3_run/states/accepted_002_tau_0p0291158703817.npz', 'abs_physical_mode_cosine': 0.9999995680544816, 'anchor': {'label': 'declared_anchor', 'tau': 0.02911625, 'status': 'MEASURED', 'smallest_singular_values': [0.17269727105561206, 0.15578351370951643, 0.11134898043597898, 0.10298097396918486, 0.08522684958138727, 0.05998703863387571, 0.05876238512610493, 0.02410583987745364, 0.01723150165432005, 7.804139703532174e-05], 'sigma_min': 7.804139703532174e-05, 'sigma_max': 9.820931089682327, 'sigma_gap_ratio_next_over_min': 220.79950268600427, 'cluster_call': 'ISOLATED_NEAR_NULL', 'sector_energy_fractions': {'psi': 0.00115689549092288, 'A': 3.8515697894105776e-08, 'r0': 0.9940375383262159, 'mu': 0.004805527667163325}, 'fold_transversality': {'wT_Ftau': 0.4859017832364471, 'normalized_abs_wT_Ftau': 0.9175629363168817, 'call': 'FOLD_TRANSVERSAL'}, 'gauge_tests': {'scaled_removed_gauge_projection_fraction': 2.1775715599688314e-33, 'independent_expanded_candidate_projection_fraction': 0.8109343699020704, 'curl_mode_l2': 3.173633255722171e-15, 'scaled_residual_response_l2': 7.804139703528706e-05, 'high_physical_overlap_alone_is_not_used': True}, 'localization': {'throat_energy_fraction_first_quarter_w': 0.02215507934700143, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}, 'mode_artifact': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/mode_comparison/C2/declared_anchor_mode.npz', 'jacobian_source': 'stage1_solver.preconditioners.assemble_closed_coupled_autodiff_sparse_jacobian', 'gauge_rank': 511, 'redundant_gauge_rank_check': 511}, 'fresh_deepest': {'label': 'fresh_deepest', 'tau': 0.029115870381663817, 'status': 'MEASURED', 'smallest_singular_values': [0.17275637285175596, 0.15576125672127766, 0.1113489803998403, 0.10297902407629507, 0.0852268495766634, 0.06005920497682463, 0.05876261405394003, 0.024090147345133546, 0.0172315016062139, 2.4414942654765485e-05], 'sigma_min': 2.4414942654765485e-05, 'sigma_max': 9.82093108969053, 'sigma_gap_ratio_next_over_min': 705.776861730639, 'cluster_call': 'ISOLATED_NEAR_NULL', 'sector_energy_fractions': {'psi': 0.0011124789981277934, 'A': 3.7447185994762925e-08, 'r0': 0.9941522601549254, 'mu': 0.004735223399760956}, 'fold_transversality': {'wT_Ftau': -0.48388353908764653, 'normalized_abs_wT_Ftau': 0.9180594809548304, 'call': 'FOLD_TRANSVERSAL'}, 'gauge_tests': {'scaled_removed_gauge_projection_fraction': 2.1311327712215682e-33, 'independent_expanded_candidate_projection_fraction': 0.8111390107092907, 'curl_mode_l2': 2.601638262057067e-15, 'scaled_residual_response_l2': 2.4414942654788537e-05, 'high_physical_overlap_alone_is_not_used': True}, 'localization': {'throat_energy_fraction_first_quarter_w': 0.022171425759245196, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}, 'mode_artifact': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/mode_comparison/C2/fresh_deepest_mode.npz', 'jacobian_source': 'stage1_solver.preconditioners.assemble_closed_coupled_autodiff_sparse_jacobian', 'gauge_rank': 511, 'redundant_gauge_rank_check': 511}}

## Resolution Evidence

- The 24x24 reproduction is physical-supporting evidence from `reports/pathA_C0g_B3_deeppoint.md` Stage-3 and `runs/pathA_C0g_B3_deeppoint/pathA_C0g_B3_deeppoint.json`; this report uses the evidence, not the stale prior insufficiency label.

## Predictor / Fallback Audit

```yaml
{'secant_fallback_available_count': 5, 'fallback_used_count': 0, 'fallback_failed_reasons': [], 'predictor_sources': ['tangent'], 'attempted_flag_without_actual_use': False}
```

## Scope Guard

```yaml
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
convergence_judged_only_on_original_residual: True
branch_acceptance_tolerance: 1e-06
bordered_arclength_residual_is_not_acceptance_arbiter: True
B4_direction_tangent_svd_only: True
gauge_Q_perp_rebuilt_inside_every_solve: True
patha_closed_branch_residual_touched: False
faithful_operators_touched: False
xi_or_grad_div_penalty_touched: False
physical_export_permitted_touched: False
frozen_physics_touched: False
no_LM_no_PTC_no_Sobolev_regularization: True
```

## Git Diff Summary

```
No diff in src/stage1_solver/coupled_branch.py, src/stage1_solver/operators.py, or src/stage1_solver/preconditioners.py.
```

Machine JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/pathA_C0g_B3_run_pseudoarclength.json`
Progress JSONL: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_run_pseudoarclength/progress.jsonl`
