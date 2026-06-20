# Path-A C0b Wall Diagnosis

Verdict: **DIAGNOSTIC_INCOMPLETE**
Deepest converged tau: `3.000000e-02` versus prior floor `2.900000e-02`; R0_min: `7.478362e-01`.

## Single Arbiter Boundary

Newton merit, line search, and convergence are evaluated with `stage1_solver.coupled_branch.patha_closed_branch_residual`. The C0 matter epsilon and k1 radius floor are used only while assembling the preconditioner. The final accepted state for every converged row has all epsilon aids inactive and is checked against the original residual.

## Tau Attempts

| tau | delta_tau | backtrack | init | used_b2c | eps_probes | R0_min | min_rho | residual | converged | sigma_status | sigma_method | sigma_min | cond_ratio | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | - | 0 | existing_b2c_background_target_predictor | true | 3 | 7.478362e-01 | 6.948664e-06 | 3.387974e-10 | true | MEASURED | dense_svd | 2.064617e-15 | 9.724165e-02 | C0 tau attempt completed |
| 2.899000e-02 | -1.010000e-03 | 0 | previous_c0_converged_state | false | 3 | 7.990414e-01 | 7.247820e-06 | 5.322040e-05 | false | MEASURED | dense_svd | 1.435880e-12 | 1.010646e+00 | tau=2.899000000000e-02, eos_K=1.800000000000e-01, epsilon=0 failed after full schedule: maximum Newton iterations reached |

Notes:

- The residual column is the original unscaled physical residual from `patha_closed_branch_residual`; no scaled or preconditioner residual is used for convergence.
- The depth sequence includes failed attempts below the prior `tau ~= 0.029` floor.
- Below-floor attempts are required to show `init.source=previous_c0_converged_state` and `used_existing_b2c=false` in the machine JSON.

## True Sigma Diagnostics

```yaml
status: MEASURED
shallow:
  tau: 0.03
  status: MEASURED
  sigma_method: dense_svd
  sigma_min: 2.064616665800032e-15
  triplet_residual_rel: 6.975080485356868e-16
  v_min_energy_fractions: {'field': 1.0, 'R0': 7.273224323143215e-25, 'mu': 4.583896994797658e-27}
  u_min_energy_fractions: {'pde_rows': 1.0, 'wall_rows': 1.3882453944442058e-25, 'mass_constraint_row': 1.8638502075678564e-24}
  condition_field_block: 3.890361604709551e+16
  condition_full_bordered: 4.0007155223076154e+17
  cond_ratio: 0.09724164547609694
  gauge_projection_status: not_available
deepest:
  tau: 0.02899
  status: MEASURED
  sigma_method: dense_svd
  sigma_min: 1.4358797209561496e-12
  triplet_residual_rel: 6.859463053113423e-16
  v_min_energy_fractions: {'field': 1.0, 'R0': 1.548388527534765e-21, 'mu': 7.582746566121585e-24}
  u_min_energy_fractions: {'pde_rows': 1.0, 'wall_rows': 4.534003152745185e-24, 'mass_constraint_row': 3.7898946678067944e-23}
  condition_field_block: 581377447366182.5
  condition_full_bordered: 575253192991118.4
  cond_ratio: 1.0106461892774903
  gauge_projection_status: not_available
```

## Tracked-Triplet Fold Test

```yaml
status: MEASURED
call: NO_TRACKED_COMPLEMENT_CROSSING
persistent_track_indices: []
fold_track_indices: []
determinant_sign_support: {'status': 'NOT_USED', 'reason': 'singular-triplet tracking is primary; det sign omitted unless separately stability-checked'}
tracked_modes:
  - track_index: 0
    sigma_ratio: 19.313232685788968
    persistent_border_or_mass_mode: False
    fold_like_complement_mode: False
    points:
      - tau: 0.03, sigma: 7.434693840833268e-14, overlap: None, v_group: field, v_fraction: 1.0
      - tau: 0.02899, sigma: 1.43587972095615e-12, overlap: 0.9970816561017359, v_group: field, v_fraction: 1.0
  - track_index: 1
    sigma_ratio: 1.0041741689405186
    persistent_border_or_mass_mode: False
    fold_like_complement_mode: False
    points:
      - tau: 0.03, sigma: 4.7006989691328044e-11, overlap: None, v_group: field, v_fraction: 1.0
      - tau: 0.02899, sigma: 4.720320480768486e-11, overlap: 0.9971921023511197, v_group: field, v_fraction: 1.0
  - track_index: 2
    sigma_ratio: 1.0205973759436986
    persistent_border_or_mass_mode: False
    fold_like_complement_mode: False
    points:
      - tau: 0.03, sigma: 2.6225709826789087e-10, overlap: None, v_group: field, v_fraction: 1.0
      - tau: 0.02899, sigma: 2.569643078156986e-10, overlap: 0.9998088664856075, v_group: field, v_fraction: 1.0
  - track_index: 3
    sigma_ratio: 1.0014441913976562
    persistent_border_or_mass_mode: False
    fold_like_complement_mode: False
    points:
      - tau: 0.03, sigma: 5.926280559311367e-10, overlap: None, v_group: field, v_fraction: 1.0
      - tau: 0.02899, sigma: 5.917734218459453e-10, overlap: 0.9998446650537844, v_group: field, v_fraction: 1.0
  - track_index: 4
    sigma_ratio: 1.0009378599241157
    persistent_border_or_mass_mode: False
    fold_like_complement_mode: False
    points:
      - tau: 0.03, sigma: 2.8112679473453323e-08, overlap: None, v_group: field, v_fraction: 1.0
      - tau: 0.02899, sigma: 2.8086338422231964e-08, overlap: 0.9999827340422915, v_group: field, v_fraction: 1.0
```

## R0-vs-Tau

```yaml
status: MEASURED
statement: tau_decrease_does_not_approach_empty_core_on_measured_crawl
tau_decrease_approaches_empty_core: False
candidate_knobs: UNKNOWN_NOT_MEASURED
trend:
  - tau: 0.03, R0_min: 0.7478361661225943, min_rho: 6.948664305210435e-06, converged: True
  - tau: 0.02899, R0_min: 0.7990414237187238, min_rho: 7.247820476152243e-06, converged: False
```

## Aid Admissibility

```yaml
residual_equality_status: PASS
residual_equality_max_abs: 0.0
epsilon_independence_status: NOT_MEASURED
residual_equality_table:
  - label: tame, tau: 0.03, max_abs: 0.0, status: PASS
  - label: near_floor, tau: 0.03, max_abs: 0.0, status: PASS
  - label: deepest_attempted, tau: 0.02899, max_abs: 0.0, status: PASS
epsilon_independence_table:
  - epsilon: 0.08, residual: 3.3879743455145217e-10, state_group_deltas: {'psi': 0.0, 'A': 0.0, 'R0': 0.0, 'mu': 0.0}, observable_deltas: {'D0_inputs': {'density_linf_delta': 0.0, 'R0_linf_delta': 0.0, 'status': 'MEASURED'}, 'BZN_overlaps': {'status': 'NOT_MEASURED', 'reason': 'BdG/Maxwell overlap recomputation is outside this bounded diagnostic'}}, status: NOT_MEASURED
  - epsilon: 0.02, residual: 3.3879743455145217e-10, state_group_deltas: {'psi': 0.0, 'A': 0.0, 'R0': 0.0, 'mu': 0.0}, observable_deltas: {'D0_inputs': {'density_linf_delta': 0.0, 'R0_linf_delta': 0.0, 'status': 'MEASURED'}, 'BZN_overlaps': {'status': 'NOT_MEASURED', 'reason': 'BdG/Maxwell overlap recomputation is outside this bounded diagnostic'}}, status: NOT_MEASURED
  - epsilon: 0.0, residual: 3.3879743455145217e-10, state_group_deltas: {'psi': 0.0, 'A': 0.0, 'R0': 0.0, 'mu': 0.0}, observable_deltas: {'D0_inputs': {'density_linf_delta': 0.0, 'R0_linf_delta': 0.0, 'status': 'MEASURED'}, 'BZN_overlaps': {'status': 'NOT_MEASURED', 'reason': 'BdG/Maxwell overlap recomputation is outside this bounded diagnostic'}}, status: NOT_MEASURED
```

C0-1 and C0-2 are implemented as preconditioner-only epsilon aids. The conditioned residual used for the approximate inverse equals the physical residual at `core_epsilon=0` and `k1_radius_epsilon=0`; accepted rows are all evaluated with those epsilons inactive. No `log rho`, `sqrt rho`, density-only variable rewrite, faithful operator edit, or final-state k1 clamp is used.

Epsilon schedule used on each target:

| core_epsilon | k1_radius_epsilon | status |
| --- | --- | --- |
| 8.000000e-02 | 8.000000e-01 | path-only preconditioner aid |
| 2.000000e-02 | 4.000000e-01 | path-only preconditioner aid |
| 0.000000e+00 | 0.000000e+00 | final physical-limit polish |

Verdict gates use the exact C0b thresholds: `SPIKE_SUFFICIENT` requires a genuine below-floor convergence; structural near-null requires a border/gauge dominant fraction `>=0.7` and `cond_ratio<0.1`; production-solver required needs field dominance `>=0.7`, `cond_ratio>=0.1`, measured sigma, and no tracked fold.

## Tau Attempt Details

### tau=3.000000000000e-02

```yaml
target_tau: 3.000000000000e-02
nominal_target_tau: 0.03
delta_tau: None
backtrack_index: 0
start_from_tau: None
init: {'path': 'software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p03.json', 'source': 'existing_b2c_background_target_predictor'}
used_existing_b2c: True
final_original_residual_linf: 3.387974345515e-10
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: True
min_rho: 6.948664305210e-06
min_R0: 7.478361661226e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | start | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 3.387974e-10 | - | - | - | - | clean_tau_attempt_state | initial residual met tolerance |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 0 | 3.387974e-10 | - | - | - | - | previous_converged_epsilon_state | initial residual met tolerance |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 0 | 3.387974e-10 | - | - | - | - | previous_converged_epsilon_state | initial residual met tolerance |

### tau=2.899000000000e-02

```yaml
target_tau: 2.899000000000e-02
nominal_target_tau: 0.02899
delta_tau: -0.0010100000000000005
backtrack_index: 0
start_from_tau: 0.03
init: {'previous_tau': 0.03, 'source': 'previous_c0_converged_state', 'wall_predictor': {'applied': True, 'initial_wall_linf': 0.0004058400535512773, 'message': 'The solution converged.', 'method': 'scipy.root(hybr)', 'predicted_wall_linf': 6.852157730108388e-16, 'r0_max': 0.9860582941036897, 'r0_min': 0.7788412355446299}}
used_existing_b2c: False
final_original_residual_linf: 5.322040363521e-05
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: False
min_rho: 7.247820476152e-06
min_R0: 7.990414237187e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | start | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 2 | 5.322041e-05 | 7.812500e-03 | 4.477552e-01 | 10 | 1.000000e+00 | clean_tau_attempt_state | maximum Newton iterations reached |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 2 | 5.322040e-05 | 7.812500e-03 | 4.477552e-01 | 6 | 1.000000e+00 | clean_tau_attempt_state | maximum Newton iterations reached |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 2 | 5.322040e-05 | 7.812500e-03 | 4.477552e-01 | 7 | 1.000000e+00 | clean_tau_attempt_state | maximum Newton iterations reached |

## Verdict Support

```yaml
reason: measured_evidence_did_not_satisfy_exact_substantive_verdict_gates
crawl_persistent_failure_evidence: False
below_floor_failed_attempt_count: 1
attempted_backtracking: False
full_newton_budget: False
fold_call: NO_TRACKED_COMPLEMENT_CROSSING
dominant_v_group: field
dominant_v_fraction: 1.0
dominant_u_group: pde_rows
dominant_u_fraction: 1.0
cond_ratio: 1.0106461892774903
```

## Recommended Next Step

Review the measured diagnostics and tighten the next bounded diagnostic around the unresolved gate.

Machine artifact: `software/stage1_solver/runs/pathA_C0b_wall_diagnosis/pathA_C0b_diagnostic.json`.
