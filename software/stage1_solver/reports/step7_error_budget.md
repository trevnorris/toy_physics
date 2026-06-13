# Step 7 Numerical Error Budget

Overall engineering gate: PASS
Config hash: `69b62c686ca71ff2`
Diagnostics digest: `aef3d65ff750d973`

**Scope framing:** this is a numerical-validity budget on the Step 4 target-blind surrogate observables only. It is not a physical section-G observable extraction, uses no extraction map or section-H target, emits no physical packet, and keeps `physical_export_permitted = False`.

**Honest limitations:** the budget covers only discretization, Newton-GMRES solver floor, boundary truncation, and conservation/Gauss closure. It excludes free_choice / physical-parameter uncertainty (GATE A: none frozen, and this will dominate eventual real observables) and model-form / parent-action uncertainty. No precision is claimed below the solver floor; null spatial sectors remain null diagnostics.

## Component Noise Floors

| component | value | units | source_step | source_commit | source_field |
| --- | --- | --- | --- | --- | --- |
| solver | 1.068892e-08 | absolute | Step 4 | c9b8b2c | noise_floor.preliminary_numerical_floor |
| discretization | 3.477501e-04 | relative raw-field self-difference; per-observable u_disc uses Step 4 finest_error_estimate | Step 4 | c9b8b2c | noise_floor.last_raw_field_relative_l2_change and observable_summary[*].finest_error_estimate |
| boundary | 6.850290e-05 | relative interior raw-field L2 | Step 5 | 63cd885 | boundary_error_metric.max_relative_l2 |
| conservation | 8.237735e-03 | relative Gauss closure | Step 6 | 4a03797 | max finest-level sponge_on independent Gauss-closure relative_residual |

Governing floor by observable class:

| observable_class | governing_component | governing_observable | governing_u_total | null_class |
| --- | --- | --- | --- | --- |
| branch scalar | discretization | chemical_potential | 2.871536e-04 | false |
| charge/mass integral | conservation | density_mass | 8.238020e-03 | false |
| energy coupled | conservation | field_energy_like_integral | 1.684863e-02 | false |
| gauge/Maxwell coupled | conservation | scalar_gauge_l2 | 1.926349e-04 | false |
| mass/density pointwise | discretization | peak_density | 5.548176e-05 | false |
| null spatial current | null | spatial_current_l2 | 0.000000e+00 | true |
| null spatial gauge | null | spatial_gauge_l2 | 0.000000e+00 | true |
| raw-field aggregate | boundary | raw_field_l2_norm | 6.852170e-05 | false |
| solver-floor diagnostic | solver | final_residual_linf | 1.068892e-08 | false |

## Combination Rule

RSS of independent numerical axes (discretization, boundary, conservation where applicable), then floored at the Step 4 solver/Newton-GMRES floor.
Grid self-convergence, boundary truncation, and Gauss/conservation closure are treated as independent numerical axes. The solver floor is a hard reporting floor, not an independent axis to RSS twice.
For sensitivity, each row also reports max(u_solver, u_disc, u_boundary, u_conservation).
Boundary relative floors are converted as `u_boundary(O) = boundary_rel * |O|` for non-null solution observables except the residual diagnostic. Conservation relative floors are converted as `u_conservation(O) = conservation_rel * |O|` only for the charge/mass integral, scalar-gauge/Maxwell, and energy-coupled surrogates.
For `density_mass`, the conservation floor is scoped through the localized Gauss-law proportionality `surface_integral Z F^{i0} dA = mu0*q_star*volume_integral rho dV`: the charge density `J^0 = q_star*rho` uses the same `rho` field as the mass integral `volume_integral rho dV`, so the Gauss-closure relative residual genuinely bounds that integral rather than switching categories.
`floor_limited = true` means the Step 4 finite-grid estimate is below the solver floor; the composed total can still be larger when boundary or conservation components also apply.
For `density_mass`, `verdict = solver-floor diagnostic` and `floor_limited = true` describe the Step 4 convergence/discretization behavior only; the composed Step 7 budget is conservation-dominated (`u_total = 8.238020e-03`), so the row is not a claim of ~1e-8 precision.
Near-zero pointwise/residual surrogates such as `min_density` and `final_residual_linf` can carry relative uncertainty near or above 1; this is the budget correctly declining to claim precision on a near-zero quantity, not a defect.

## Per-Observable Budget

| label | observable_class | finest_value | u_solver | u_disc | u_boundary | u_conservation | u_total | relative_uncertainty | dominant_component | verdict | floor_limited | u_max_alternative | rss_over_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| density mass integral | charge/mass integral | 1.000000e+00 | 1.068892e-08 | 2.305045e-12 | 6.850290e-05 | 8.237735e-03 | 8.238020e-03 | 8.238020e-03 | conservation | solver-floor diagnostic | true | 8.237735e-03 | 1.000035e+00 |
| peak density | mass/density pointwise | 2.800545e-01 | 1.068892e-08 | 5.205938e-05 | 1.918455e-05 | 0.000000e+00 | 5.548176e-05 | 1.981106e-04 | discretization | expected-order convergence | false | 5.205938e-05 | 1.065740e+00 |
| minimum density | mass/density pointwise | 2.309396e-07 | 1.068892e-08 | 2.284870e-07 | 1.582003e-11 | 0.000000e+00 | 2.284870e-07 | 9.893798e-01 | discretization | expected-order convergence | false | 2.284870e-07 | 1.000000e+00 |
| raw field L2 norm | raw-field aggregate | 1.000273e+00 | 1.068892e-08 | 1.142831e-07 | 6.852161e-05 | 0.000000e+00 | 6.852170e-05 | 6.850299e-05 | boundary | expected-order convergence | false | 6.852161e-05 | 1.000001e+00 |
| A0 L2 norm | gauge/Maxwell coupled | 2.337608e-02 | 1.068892e-08 | 4.898945e-06 | 1.601329e-06 | 1.925659e-04 | 1.926349e-04 | 8.240685e-03 | conservation | expected-order convergence | false | 1.925659e-04 | 1.000358e+00 |
| spatial gauge L2 norm | null spatial gauge | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - | null | null diagnostic | false | 0.000000e+00 | 1.000000e+00 |
| spatial current L2 norm | null spatial current | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - | null | null diagnostic | false | 0.000000e+00 | 1.000000e+00 |
| field-energy-like integral | energy coupled | 2.044752e+00 | 1.068892e-08 | 3.632239e-04 | 1.400715e-04 | 1.684413e-02 | 1.684863e-02 | 8.239935e-03 | conservation | expected-order convergence | false | 1.684413e-02 | 1.000267e+00 |
| chemical potential | branch scalar | 2.093577e+00 | 1.068892e-08 | 2.487751e-04 | 1.434161e-04 | 0.000000e+00 | 2.871536e-04 | 1.371594e-04 | discretization | expected-order convergence | false | 2.487751e-04 | 1.154270e+00 |
| final residual Linf | solver-floor diagnostic | 9.284885e-09 | 1.068892e-08 | 8.704916e-12 | 0.000000e+00 | 0.000000e+00 | 1.068892e-08 | 1.151217e+00 | solver | solver-floor diagnostic | true | 1.068892e-08 | 1.000000e+00 |

## RSS-vs-Max Sensitivity

Rows where RSS exceeds the conservative max alternative by the configured 1.050000e+00 ratio are surfaced:

| label | u_total | u_max_alternative | rss_over_max | dominant_component |
| --- | --- | --- | --- | --- |
| peak density | 5.548176e-05 | 5.205938e-05 | 1.065740e+00 | discretization |
| chemical potential | 2.871536e-04 | 2.487751e-04 | 1.154270e+00 | discretization |

## Pass Checks

Counted gates only; `passed` is the conjunction of these checks.

- non_null_uncertainties_floor_at_solver: PASS
- solver_floor_limited_rows_flagged: PASS
- null_sectors_remain_null: PASS
- conservation_floor_scoped: PASS
- boundary_floor_scoped: PASS
- observable_set_matches_step4: PASS

Asserted-by-construction checks, reported but not counted as physics gates:
- target_blind_surrogate_budget_only_not_a_physics_gate: PASS - Step 7 composes only target-blind surrogate observables and imports no target/reference/export map; asserted by module construction, not an independent physics measurement.
- physical_export_permitted_is_false_not_a_physics_gate: PASS - Step 7 emits no physical packet; the external export guard lives in the firewalled physical model and is untouched - asserted by construction, not read here.
- prior_step4_passed_not_a_physics_gate: PASS - Default pinned provenance records Step 4 PASS at commit c9b8b2c; this is not a live Step 7 rerun.
- prior_step5_passed_not_a_physics_gate: PASS - Default pinned provenance records Step 5 PASS at commit 63cd885; this is not a live Step 7 rerun.
- prior_step6_passed_not_a_physics_gate: PASS - Default pinned provenance records Step 6 PASS at commit 4a03797; this is not a live Step 7 rerun.

## Provenance

- step4: `stage1_solver.convergence.run_step4` at `c9b8b2c`, report `software/stage1_solver/reports/step4_convergence_study.md`, recorded input `code-pinned run_step4 output snapshot`.
- step5: `stage1_solver.boundary_characterization.run_step5` at `63cd885`, report `software/stage1_solver/reports/step5_boundary_characterization.md`, recorded input `code-pinned run_step5 output scalar`.
- step6: `stage1_solver.conservation_diagnostics.run_step6` at `4a03797`, report `software/stage1_solver/reports/step6_conservation_diagnostics.md`, recorded input `code-pinned run_step6 output scalar`.

## Reproduction

```bash
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.error_budget_harness
timeout 600 env PYTHONPATH=software/stage1_solver/src pytest software/stage1_solver/tests/test_stage1_solver.py
```

Machine-readable table: `software/stage1_solver/runs/step7_error_budget/step7_error_budget/error_budget_table.json`.
Target-blindness: Step 7 imports no benchmark targets, no references, no extraction map, and no physical export path; relative scales are the observable's own magnitude or prior target-blind aggregate norms.
Export guard: `physical_export_permitted` remains false; no physical packet is emitted.

