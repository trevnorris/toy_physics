# Path-A Chunk 1c Self-Consistent Closed Validation

Overall engineering gate: PASS
Config hash: `6ca201495f2b28c2`
S_Sigma digest: `63501899d5957528`

Scope: target-blind engineering validation of the closed placeholder solve. No calibration, physical packet, coefficient extraction, or export is performed.

## Method

Closed placeholder label: Path-A chunk-1b closed engineering placeholder; no physical packet; target-blind
Refinement ladder: `((4, 4), (8, 8), (16, 16))` with ratio `2`.
Preconditioner rebuild policy for validation ladder: `once_per_continuation`. The true matrix-free JVP residual remains the solve operator.

## Self-Consistent Balance Terms

The counted balance check here is term non-degeneracy. Residual-floor rows are solve diagnostics because the Newton residual already contains the wall residual. The source and conservative flux recomputes are parallel reconstructions of the same kernels: useful edit-drift pins, but not physics gates.

| term | linf | l2 | relative_to_dominant |
| --- | --- | --- | --- |
| flux_divergence | 6.391311e-02 | 6.915119e-02 | 1.000000e+00 |
| gradient_square | 3.283017e-05 | 1.684117e-05 | 5.136688e-04 |
| potential_gradient | 5.794288e-02 | 6.286577e-02 | 9.065883e-01 |
| source | 5.998099e-03 | 6.279460e-03 | 9.384771e-02 |

Dominant term Linf `6.391311e-02`; balance residual Linf `6.294184e-14`; residual relative to dominant `9.848032e-13`; solver-floor reference `2.000000e-09`. The nontrivial-term floor `1.000000e-08` is a degeneracy tripwire, not a calibrated non-degeneracy scale.
Edit-drift pins: source max_abs_diff `0.000000e+00`; flux-divergence max_abs_diff `0.000000e+00`.
For the source, no alternate discretization is meaningful because it is an exact reduced source. Correctness is supported by reciprocity with the forward wall-to-matter force, the chunk-1a dual-engine MMS, and the fidelity audit; the row above only pins wiring/edit drift.

Counted balance gates:
| gate | status |
| --- | --- |
| balance_terms_nontrivial | PASS |

Balance identity checks, reported but not counted as physics gates:
| gate | status |
| --- | --- |
| independent_source_recompute_matches_not_a_physics_gate | PASS |
| independent_flux_recompute_matches_not_a_physics_gate | PASS |
| balance_residual_relative_below_solver_floor_not_a_physics_gate | PASS |

## Closed Grid-Convergence Ladder

| level | grid | dof | wall_clock_seconds | converged | final_residual_linf | final_tolerance | r0_min | r0_max | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | closed_l0_nr_4_nw_4 | 85 | 1.822957e+01 | true | 1.186493e-09 | 2.000000e-09 | 9.878671e-01 | 9.970928e-01 | closed validation continuation completed |
| 1 | closed_l1_nr_8_nw_8 | 329 | 4.102487e+01 | true | 1.456064e-09 | 2.000000e-09 | 9.876428e-01 | 9.984848e-01 | closed validation continuation completed |
| 2 | closed_l2_nr_16_nw_16 | 1297 | 4.099620e+01 | true | 1.532886e-09 | 2.000000e-09 | 9.875850e-01 | 9.992279e-01 | closed validation continuation completed |

Raw-field self-difference, including R0:

| coarse_grid | fine_grid | raw_field_relative_l2_change | r0_relative_l2_change | raw_field_linf_change | r0_linf_change | observed_order_from_raw_l2_change | observed_order_from_r0_l2_change |
| --- | --- | --- | --- | --- | --- | --- | --- |
| closed_l0_nr_4_nw_4 | closed_l1_nr_8_nw_8 | 5.095828e-03 | 2.475113e-05 | 3.238039e-03 | 3.213911e-05 | - | - |
| closed_l1_nr_8_nw_8 | closed_l2_nr_16_nw_16 | 1.258458e-03 | 6.942618e-06 | 9.422037e-04 | 9.553670e-06 | 2.017669e+00 | 1.833953e+00 |

Independent flux discretization check:

The non-conservative reconstruction uses center gradients only: `g = wall_center_gradient(R0)`, `q = T_w(R0, w_center) * g`, `-wall_center_gradient(q)`. It is compared with the conservative face-flux operator on interior wall cells only; the mouth one-sided stencil and zero-traction exit cell are excluded because those closures intentionally differ.

| grid | spacing | interior_cells_compared | interior_max_abs_diff | observed_order | h2_scaled_tolerance | decreases_or_at_floor | finest_within_tolerance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| nr_4_nw_4 | 3.000000e-01 | 2 | 1.134133e-03 | - | 4.530294e-02 | true | - |
| nr_8_nw_8 | 1.500000e-01 | 6 | 6.006613e-04 | 9.169667e-01 | 1.149867e-02 | true | - |
| nr_16_nw_16 | 7.500000e-02 | 14 | 2.063890e-04 | 1.541186e+00 | 2.879399e-03 | true | true |

R0-aware surrogate observables:

| observable | finest_grid | finest_value | last_observed_order | richardson_estimate | finest_error_estimate |
| --- | --- | --- | --- | --- | --- |
| density_mass | closed_l2_nr_16_nw_16 | 5.000000e-02 | -6.078231e-01 | - | - |
| raw_matter_gauge_l2_norm | closed_l2_nr_16_nw_16 | 2.236107e-01 | 2.102767e+00 | 2.236107e-01 | 1.880148e-08 |
| r0_l2_norm | closed_l2_nr_16_nw_16 | 1.086359e+00 | 1.340812e+00 | 1.086357e+00 | 1.950569e-06 |
| r0_min | closed_l2_nr_16_nw_16 | 9.875850e-01 | 1.956428e+00 | 9.875650e-01 | 2.005921e-05 |
| r0_max | closed_l2_nr_16_nw_16 | 9.992279e-01 | 9.054881e-01 | 1.000079e+00 | 8.510334e-04 |
| r0_range | closed_l2_nr_16_nw_16 | 1.164287e-02 | 1.012967e+00 | 1.242955e-02 | 7.866886e-04 |
| chemical_potential | closed_l2_nr_16_nw_16 | 2.037131e+00 | 1.938572e+00 | 2.041208e+00 | 4.076902e-03 |
| final_residual_linf | closed_l2_nr_16_nw_16 | 1.532886e-09 | 1.811077e+00 | 1.563504e-09 | 3.061800e-11 |
| balance_residual_relative_to_dominant | closed_l2_nr_16_nw_16 | 9.848032e-13 | 4.051379e-01 | -2.101321e-11 | 2.199801e-11 |
| mass_request | closed_l2_nr_16_nw_16 | 5.000000e-02 | - | - | - |

Laptop-limit stop reason: completed configured ladder.

Convergence gates:
| gate | status |
| --- | --- |
| minimum_three_levels | PASS |
| all_completed_levels_converged | PASS |
| fixed_ratio_ladder | PASS |
| raw_field_self_difference_decreases_or_at_floor | PASS |
| independent_flux_discretization_converges | PASS |

Convergence identity checks, reported but not counted as physics gates:
| gate | status |
| --- | --- |
| balance_residual_relative_stable_under_refinement_not_a_physics_gate | PASS |

## Closed Conservation Diagnostics

Number/charge and the independent Gauss closure are evaluated on the solved matter/gauge fields. The energy-density confinement term uses frozen geometry because `conservation_diagnostics._energy_density` calls `confinement_potential_torch` without `radius=`; the closed-vs-frozen comparison is apples-to-apples for that channel because both branches use the same frozen confinement term.

The 4x frozen-baseline envelope is a sanity bound for the closed solve. The live conservation signal is the independent center-gradient Gauss reconstruction; structurally-zero channels are disclosed separately.
Closed/frozen metric threshold factor: `4.0`.
null_floor_label keys on transport-signal norm, not residual magnitude; that is correct for this isotropic stationary branch and should be revisited for a current-carrying branch.

Frozen baseline solve: grid `frozen_baseline_nr_16_nw_16`, converged=True, final_residual_linf=1.560697e-09, wall_clock=3.410430e+01s.

| metric | activity | closed_r0_solved | frozen_geometry_baseline | threshold | basis | status |
| --- | --- | --- | --- | --- | --- | --- |
| local_l2_relative_max | STRUCTURALLY ZERO | 0.000000e+00 | 0.000000e+00 | 1.000000e-10 | stationary current sectors have zero transport signal on this branch | PASS |
| local_linf_relative_max | STRUCTURALLY ZERO | 0.000000e+00 | 0.000000e+00 | 1.000000e-10 | stationary current sectors have zero transport signal on this branch | PASS |
| independent_gauss_relative_max | LIVE | 2.839914e-03 | 2.846931e-03 | 1.138772e-02 | center-gradient Gauss reconstruction independent of operator fluxes | PASS |
| operator_maxwell_absolute_max | SOLVER-FLOOR IDENTITY | 1.118897e-16 | 7.118438e-15 | 2.941845e-09 | operator Gauss residual is bounded by the converged Newton residual | PASS |
| fv_budget_closure_absolute_max | STRUCTURALLY ZERO | 0.000000e+00 | 0.000000e+00 | 1.000000e-10 | finite-volume budget identity closes algebraically for zero flux/sink | PASS |

Independent center-gradient Gauss closure rows:

| surface | reconstruction | relative_residual | absolute_residual | surface_flux | enclosed_mu0_charge |
| --- | --- | --- | --- | --- | --- |
| nested_surface_0 | independent_center_gradient | 2.839914e-03 | -7.243665e-06 | 2.543420e-03 | 2.550664e-03 |
| nested_surface_1 | independent_center_gradient | 3.547494e-04 | 1.759497e-06 | 4.961590e-03 | 4.959831e-03 |

Conservation gates:
| gate | status |
| --- | --- |
| closed_conservation_consistent_with_frozen_baseline | PASS |
| stationary_current_sectors_labelled_null_or_floor | PASS |

Identity checks, reported but not counted as physics gates:
| gate | status |
| --- | --- |
| fv_budget_identity_roundoff_not_a_physics_gate | PASS |
| closed_operator_gauss_residual_at_solver_floor_not_a_physics_gate | PASS |

## Counted Gates

| gate | status |
| --- | --- |
| balance_terms_nontrivial | PASS |
| minimum_three_levels | PASS |
| all_completed_levels_converged | PASS |
| fixed_ratio_ladder | PASS |
| raw_field_self_difference_decreases_or_at_floor | PASS |
| independent_flux_discretization_converges | PASS |
| closed_conservation_consistent_with_frozen_baseline | PASS |
| stationary_current_sectors_labelled_null_or_floor | PASS |
| target_token_scan_clean | PASS |

## Identity / Not Physics Gates

| gate | status |
| --- | --- |
| independent_source_recompute_matches_not_a_physics_gate | PASS |
| independent_flux_recompute_matches_not_a_physics_gate | PASS |
| balance_residual_relative_below_solver_floor_not_a_physics_gate | PASS |
| balance_residual_relative_stable_under_refinement_not_a_physics_gate | PASS |
| fv_budget_identity_roundoff_not_a_physics_gate | PASS |
| closed_operator_gauss_residual_at_solver_floor_not_a_physics_gate | PASS |

## Stop Conditions

| condition | status |
| --- | --- |
| edit_drift_pin_mismatch_not_a_physics_gate | not tripped |
| independent_flux_discretization_not_converging | not tripped |
| balance_terms_degenerate | not tripped |
| fewer_than_three_levels_completed | not tripped |
| r0_nonpositive | not tripped |
| closed_conservation_materially_worse_than_frozen | not tripped |

## Reproduction

```bash
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_closed_validation
timeout 600 env PYTHONPATH=software/stage1_solver/src pytest software/stage1_solver/tests -q
```

Machine-readable table: `software/stage1_solver/runs/patha_chunk1c_self_consistent_validation/patha_chunk1c_self_consistent_validation/closed_validation_table.json`.
Target-token scan: PASS (scanned 6 Path-A 1c files).

