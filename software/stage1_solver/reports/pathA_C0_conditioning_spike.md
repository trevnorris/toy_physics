# Path-A C0 Conditioning Spike

Verdict: **PRODUCTION_SOLVER_REQUIRED**
Deepest converged tau: `2.912500e-02` versus prior floor `2.900000e-02`; R0_min: `7.931268e-01`.

## Single Arbiter Boundary

Newton merit, line search, and convergence are evaluated with `stage1_solver.coupled_branch.patha_closed_branch_residual`. The C0 matter epsilon and k1 radius floor are used only while assembling the preconditioner. The final accepted state for every converged row has all epsilon aids inactive and is checked against the original residual.

## Diagnostic Table

| tau | R0_min | min_rho | residual | converged | D0 | sigma_min_unscaled | sigma_min_scaled | cond_unscaled | cond_scaled | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 7.478362e-01 | 6.948664e-06 | 3.387974e-10 | true | 2.312082e-01 | 2.070636e-18 | 1.694917e-20 | 4.643876e+20 | 8.639955e+20 | 0.000000e+00 | C0 depth target completed |
| 2.950000e-02 | 7.655102e-01 | 7.059581e-06 | 3.631852e-13 | true | - | - | - | - | - | 0.000000e+00 | C0 depth target completed |
| 2.925000e-02 | 7.790815e-01 | 7.137544e-06 | 3.596429e-14 | true | 2.255325e-01 | - | - | - | - | 0.000000e+00 | C0 depth target completed |
| 2.912500e-02 | 7.931268e-01 | 7.212166e-06 | 3.950424e-10 | true | 2.246077e-01 | 4.591642e-18 | 2.153514e-19 | 2.094191e+20 | 6.800051e+19 | 0.000000e+00 | C0 depth target completed |
| 2.906250e-02 | 7.982904e-01 | 7.240140e-06 | 2.423841e-05 | false | - | - | - | - | - | 1.000000e+00 | tau=2.906250000000e-02, eos_K=1.800000000000e-01, aid_index=0 failed: line search failed to reduce original residual |
| 2.900000e-02 | 7.987182e-01 | 7.244681e-06 | 5.355903e-05 | false | - | 7.198995e-13 | 2.091340e-14 | 1.335711e+15 | 7.002211e+14 | 1.000000e+00 | tau=2.900000000000e-02, eos_K=1.800000000000e-01, aid_index=0 failed: line search failed to reduce original residual |
| 2.850000e-02 | 8.016771e-01 | 7.257061e-06 | 1.987476e-04 | false | - | - | - | - | - | 1.000000e+00 | tau=2.850000000000e-02, eos_K=1.800000000000e-01, aid_index=0 failed: line search failed to reduce original residual |
| 2.800000e-02 | 8.057971e-01 | 7.317173e-06 | 4.918900e-04 | false | - | 4.692438e-12 | 1.491783e-13 | 2.049207e+14 | 9.816444e+13 | 1.000000e+00 | tau=2.800000000000e-02, eos_K=1.800000000000e-01, aid_index=0 failed: line search failed to reduce original residual |

Notes:

- The residual column is the original unscaled physical residual from `patha_closed_branch_residual`; no scaled or preconditioner residual is used for convergence.
- Large-grid `sigma_min_*` values are recorded with method `one_norm_lu_estimate` in the JSON artifact. They are conditioning diagnostics, not physical evidence of stability. Rows without sigma values were outside the representative linear-diagnostic tau set.
- `D0` is reported only where an existing B2c extraction/evaluation bundle exists. Failed below-floor backgrounds are not extrapolated into D0/B/Z/N values.
- The depth sequence includes failed attempts below the prior `tau ~= 0.029` floor: `0.0290625`, `0.029`, `0.0285`, and `0.028`.

## Aid Admissibility

```yaml
no_variable_transform_used: True
residual_equality_status: PASS
residual_equality_max_abs: 0.0
residual_equality_checked_state_count: 3
epsilon_independence_status: PASS
epsilon_path_final_aids_inactive: True
original_residual_gating_status: PASS
jacobi_scaling_solution_neutral: True
jacobi_scaling_note: Scaling is applied as R*J*C dz = R*(-F), step=C*dz. Line search and convergence use the unscaled original residual.
```

C0-1 and C0-2 were implemented as preconditioner-only epsilon aids. The conditioned residual used for the approximate inverse equals the physical residual at `core_epsilon=0` and `k1_radius_epsilon=0`; the final accepted rows are all evaluated with those epsilons inactive. No `log rho`, `sqrt rho`, density-only variable rewrite, faithful operator edit, or final-state k1 clamp is used.

Epsilon schedule used on each target:

| core_epsilon | k1_radius_epsilon | status |
| --- | --- | --- |
| 8.000000e-02 | 8.000000e-01 | path-only preconditioner aid |
| 2.000000e-02 | 4.000000e-01 | path-only preconditioner aid |
| 0.000000e+00 | 0.000000e+00 | final physical-limit polish |

Verdict gate: `PRODUCTION_SOLVER_REQUIRED` is issued only after C0-1/C0-2/C0-3 were active in the path and C0-4 attempted the frozen-family tau sequence through and below `0.029`; every below-floor attempt failed the original residual line search and remained above the B2c `1e-6` background tolerance.

## Depth Attempts

### tau=3.000000000000e-02

```yaml
final_original_residual_linf: 3.387974345515e-10
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: True
min_rho: 6.948664305210e-06
min_R0: 7.478361661226e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 3.387974e-10 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 0 | 3.387974e-10 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 0 | 3.387974e-10 | - | - | - | - | initial residual met tolerance |

### tau=2.950000000000e-02

```yaml
final_original_residual_linf: 3.631851763775e-13
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: True
min_rho: 7.059580505090e-06
min_R0: 7.655102446834e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 3.631852e-13 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 0 | 3.631852e-13 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 0 | 3.631852e-13 | - | - | - | - | initial residual met tolerance |

### tau=2.925000000000e-02

```yaml
final_original_residual_linf: 3.596428710395e-14
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: True
min_rho: 7.137544065798e-06
min_R0: 7.790815134569e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 3.596429e-14 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 0 | 3.596429e-14 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 0 | 3.596429e-14 | - | - | - | - | initial residual met tolerance |

### tau=2.912500000000e-02

```yaml
final_original_residual_linf: 3.950424206506e-10
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: True
min_rho: 7.212165974964e-06
min_R0: 7.931267729874e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 3.950424e-10 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 2.000000e-02 | 4.000000e-01 | 0 | 3.950424e-10 | - | - | - | - | initial residual met tolerance |
| 1.800000e-01 | 0.000000e+00 | 0.000000e+00 | 0 | 3.950424e-10 | - | - | - | - | initial residual met tolerance |

### tau=2.906250000000e-02

```yaml
final_original_residual_linf: 2.423840866221e-05
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: False
min_rho: 7.240139876803e-06
min_R0: 7.982903785831e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 2.423841e-05 | - | 5.100233e+02 | 8 | 1.000000e+00 | line search failed to reduce original residual |

### tau=2.900000000000e-02

```yaml
final_original_residual_linf: 5.355903474474e-05
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: False
min_rho: 7.244680992989e-06
min_R0: 7.987182092757e-01
k1_clamp_active_in_path: True
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 5.355903e-05 | - | 3.697922e+01 | 11 | 1.000000e+00 | line search failed to reduce original residual |

### tau=2.850000000000e-02

```yaml
final_original_residual_linf: 1.987476224519e-04
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: False
min_rho: 7.257061241202e-06
min_R0: 8.016771247464e-01
k1_clamp_active_in_path: False
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 3 | 1.987476e-04 | - | 2.737109e+02 | 11 | 1.000000e+00 | line search failed to reduce original residual |

### tau=2.800000000000e-02

```yaml
final_original_residual_linf: 4.918899520096e-04
b2c_background_tolerance: 1.000000000000e-06
final_physical_converged: False
min_rho: 7.317173456532e-06
min_R0: 8.057970723908e-01
k1_clamp_active_in_path: False
```

| eos_K | core_eps | k1_eps | iters | residual | alpha | step_norm | gmres_iters | gmres_growth | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.800000e-01 | 8.000000e-02 | 8.000000e-01 | 0 | 4.918900e-04 | - | 2.972087e+02 | 12 | 1.000000e+00 | line search failed to reduce original residual |

## Production Solver Recommendation

If this verdict is PRODUCTION_SOLVER_REQUIRED, the supported next step is a production linear solver with PETSc-style block preconditioning, multigrid or mesh grading near the wall. This C0 run does not implement that rebuild.

Machine artifact: `software/stage1_solver/runs/pathA_C0_conditioning_spike/pathA_C0_diagnostic_table.json`.
