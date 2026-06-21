# Path-A C0g B-3 Pseudo-Arclength

B-3.0 reproduction: **FAIL**
B-3.2 verdict: **NOT_MEASURED**

## Method

- Reduced unknowns are `(Q_perp coefficients, tau)`.
- Declared arclength metric: `||dx||_2^2 + (1000.0 dτ)^2`.
- Predictor tangent first solves `row_scale * J_original * col_scale * Q_perp * z = -row_scale * F_tau * dτ` in the gauge complement; the augmented null vector is the fallback.
- The bordered corrector drives the unchanged original physical residual; the arclength row is not a convergence arbiter.
- `Q_perp` is rebuilt inside every predictor and corrector linear solve.

## B-4 Reconfirm

```yaml
status: PASS
reference_tau: 0.029125
tolerances: {'max_abs': 1e-10, 'max_rel': 1e-08}
assembly_speedup_x: 71.88268385316061
```

## Target

```yaml
deepcrawl_json_path: /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_deepcrawl/pathA_C0g_deepcrawl.json
seed_tau: 0.0291225
deepest_converged_tau: 0.0291132
deepest_residual_linf: 9.798546447628564e-07
```

## B-3.0 Ladder

| accepted_index | predecessor_state_id | predictor_step | corrector_iters | solved_tau | dtauds | original_residual_linf | min_R0 | min_rho | mu | sigma_min_JQ_perp | bordered_condition | wall_seconds | not_initialized_from_comparison_artifact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | single_initial_seed_from_deepcrawl | 6.660503e-04 | 1 | 2.912201e-02 | -7.506941e-04 | 7.986823e-13 | 7.940089e-01 | 7.216658e-06 | 2.035241e+00 | 3.793844e-04 | 1.095861e+06 | 8.355959e+00 | true |
| 1 | B3_0_accepted_000 | 9.106771e-04 | 1 | 2.912134e-02 | -7.376296e-04 | 2.195315e-12 | 7.942320e-01 | 7.217790e-06 | 2.035197e+00 | 3.577292e-04 | 1.091970e+06 | 7.658131e+00 | true |
| 2 | B3_0_accepted_001 | 9.559511e-04 | 1 | 2.912067e-02 | -7.178037e-04 | 3.058640e-12 | 7.944738e-01 | 7.219015e-06 | 2.035150e+00 | 3.342613e-04 | 1.070191e+06 | 1.336931e+01 | true |
| 3 | B3_0_accepted_002 | 9.793336e-04 | 1 | 2.912000e-02 | -6.941084e-04 | 3.875412e-12 | 7.947301e-01 | 7.220313e-06 | 2.035100e+00 | 3.093830e-04 | 1.037819e+06 | 1.292904e+01 | true |
| 4 | B3_0_accepted_003 | 9.599692e-04 | 1 | 2.911938e-02 | -6.662075e-04 | 4.091911e-12 | 7.949904e-01 | 7.221629e-06 | 2.035049e+00 | 2.841198e-04 | 9.918708e+05 | 1.358610e+01 | true |
| 5 | B3_0_accepted_004 | 1.014270e-03 | 1 | 2.911875e-02 | -6.346291e-04 | 5.784874e-12 | 7.952758e-01 | 7.223069e-06 | 2.034993e+00 | 2.564235e-04 | 9.365842e+05 | 1.691769e+01 | true |
| 6 | B3_0_accepted_005 | 1.091228e-03 | 1 | 2.911813e-02 | -5.958546e-04 | 8.707592e-12 | 7.955954e-01 | 7.224679e-06 | 2.034930e+00 | 2.254085e-04 | 8.622349e+05 | 1.392642e+01 | true |
| 7 | B3_0_accepted_006 | 1.209754e-03 | 1 | 2.911750e-02 | -5.467603e-04 | 1.444410e-11 | 7.959657e-01 | 7.226541e-06 | 2.034858e+00 | 1.894852e-04 | 7.616152e+05 | 1.455224e+01 | true |
| 8 | B3_0_accepted_007 | 1.420430e-03 | 1 | 2.911688e-02 | -4.817362e-04 | 2.862671e-11 | 7.964221e-01 | 7.228830e-06 | 2.034769e+00 | 1.452064e-04 | 6.240437e+05 | 1.410423e+01 | true |
| 9 | B3_0_accepted_008 | 2.004835e-03 | 1 | 2.911626e-02 | -3.888116e-04 | 1.044314e-10 | 7.971042e-01 | 7.232241e-06 | 2.034637e+00 | 7.903553e-05 | 4.349184e+05 | 2.067382e+01 | true |

## B-3.0 Reproduction Check

```yaml
{'status': 'FAIL', 'checks': [{'accepted_id': 'B3_0_accepted_000', 'solved_tau': 0.02912200507574685, 'comparison_tau': 0.029122, 'tau_abs_error': 5.075746850602414e-09, 'mismatch': {'rho_linf': 3.5892371592211036e-09, 'curl_A_linf': 8.517760179563476e-14, 'r0_linf': 1.1126047048115595e-06, 'mu_abs': 2.1729174681794916e-07}, 'max_gauge_invariant_mismatch': 1.1126047048115595e-06, 'status': 'PASS', 'comparison_artifact_read_after_run': True}, {'accepted_id': 'B3_0_accepted_001', 'solved_tau': 0.029121342082368048, 'comparison_tau': 0.029122, 'tau_abs_error': 6.579176319507218e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_002', 'solved_tau': 0.029120666844220967, 'comparison_tau': 0.02912, 'tau_abs_error': 6.668442209672998e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_003', 'solved_tau': 0.02912000025932065, 'comparison_tau': 0.02912, 'tau_abs_error': 2.593206482881527e-10, 'mismatch': {'rho_linf': 5.730148758265052e-09, 'curl_A_linf': 1.2575927710444117e-12, 'r0_linf': 1.683469272872884e-06, 'mu_abs': 3.3278814193238304e-07}, 'max_gauge_invariant_mismatch': 1.683469272872884e-06, 'status': 'PASS', 'comparison_artifact_read_after_run': True}, {'accepted_id': 'B3_0_accepted_004', 'solved_tau': 0.029119375337953265, 'comparison_tau': 0.02912, 'tau_abs_error': 6.246620467345587e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_005', 'solved_tau': 0.02911875055962929, 'comparison_tau': 0.029118125, 'tau_abs_error': 6.255596292879806e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_006', 'solved_tau': 0.029118126013791743, 'comparison_tau': 0.029118125, 'tau_abs_error': 1.0137917415276032e-09, 'mismatch': {'rho_linf': 2.7392167423201386e-09, 'curl_A_linf': 1.4508222452511293e-13, 'r0_linf': 8.036941501199379e-07, 'mu_abs': 1.586101152639685e-07}, 'max_gauge_invariant_mismatch': 8.036941501199379e-07, 'status': 'PASS', 'comparison_artifact_read_after_run': True}, {'accepted_id': 'B3_0_accepted_007', 'solved_tau': 0.029117502102375726, 'comparison_tau': 0.029118125, 'tau_abs_error': 6.228976242757256e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_008', 'solved_tau': 0.029116880500629452, 'comparison_tau': 0.02911625, 'tau_abs_error': 6.305006294525517e-07, 'status': 'SKIPPED_INTERMEDIATE', 'reason': 'accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau', 'comparison_artifact_read_after_run': False}, {'accepted_id': 'B3_0_accepted_009', 'solved_tau': 0.029116256342628825, 'comparison_tau': 0.02911625, 'tau_abs_error': 6.342628824929086e-09, 'mismatch': {'rho_linf': 2.2422527396100123e-06, 'curl_A_linf': 5.072980484069208e-13, 'r0_linf': 0.0006820129158060517, 'mu_abs': 0.00013294364724325547}, 'max_gauge_invariant_mismatch': 0.0006820129158060517, 'status': 'FAIL', 'comparison_artifact_read_after_run': True}], 'state_continuity': [{'accepted_id': 'B3_0_accepted_001', 'state_l2_distance_from_previous': 0.0006244521449010925, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_002', 'state_l2_distance_from_previous': 0.0006768627061253064, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_003', 'state_l2_distance_from_previous': 0.0007177005022285015, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_004', 'state_l2_distance_from_previous': 0.0007289706870385514, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_005', 'state_l2_distance_from_previous': 0.0007993713242040526, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_006', 'state_l2_distance_from_previous': 0.0008954004428273273, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_007', 'state_l2_distance_from_previous': 0.0010374241128451396, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_008', 'state_l2_distance_from_previous': 0.001279198249192796, 'distinct_above_floor': True}, {'accepted_id': 'B3_0_accepted_009', 'state_l2_distance_from_previous': 0.0019126473180686842, 'distinct_above_floor': True}], 'not_initialized_from_comparison_artifact': True, 'comparison_artifacts_read_only_after_run': True, 'tolerance': 5e-05}
```

## B-3.1 Ladder

| accepted_index | predecessor_state_id | predictor_step | corrector_iters | solved_tau | dtauds | original_residual_linf | min_R0 | min_rho | mu | sigma_min_JQ_perp | bordered_condition | wall_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Verdict

```yaml
{'verdict': 'NOT_MEASURED', 'why': 'B-3.0 reproduction gate failed; hard-region continuation was not attempted', 'scientific_verdict_earned': False}
```

## Scope Guard

```yaml
B3_pseudoarclength_built: True
B4_analytic_sparse_assembly_built: True
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
convergence_judged_only_on_original_residual: True
patha_closed_branch_residual_touched: False
faithful_operators_touched: False
xi_or_grad_div_penalty_touched: False
physical_export_permitted_touched: False
prefer_existing_b2c_background_predictor: False
frozen_physics_touched: False
no_LM_no_PTC_no_Sobolev_regularization: True
```

## Git Diff Summary

```
No diff in src/stage1_solver/coupled_branch.py or src/stage1_solver/operators.py.
```

Machine JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_pseudoarclength/pathA_C0g_B3_pseudoarclength.json`
