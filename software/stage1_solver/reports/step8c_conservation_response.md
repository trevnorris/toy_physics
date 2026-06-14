# Step 8c Conservation Balance And Surrogate Response

Overall engineering gate: PASS_WITH_RESOLUTION_LIMIT
Config hash: `10e452493d025c62`
Diagnostics digest: `ac42ec665807195c`

**Scope:** engineering-smoke, target-blind, forced/CAP driven balance plus surrogate low-frequency methodology. No physical packet is exported, and the open matter/gauge-to-wall source is not supplied.

## Sources

- linearized_skeleton: notes/moving_throat_pde_program_compact.md lines 1383-1455
- wall_modal_split: notes/moving_throat_pde_program_compact.md lines 1298-1308
- p2_degeneracy: notes/moving_throat_pde_program_compact.md lines 1340-1364
- open_backreaction: notes/moving_throat_pde_program_compact.md lines 1377-1381
- wall_to_matter: notes/moving_throat_pde_program_compact.md lines 1075-1089
- bulk_angular_reduction: Step-8a directive bulk angular reduction paragraph
- maxwell_scalarization: Step-8a engineering-smoke retained-lane scalarization; the full vector-harmonic Maxwell l=2 reduction is deferred.
- time_derivative_terms: notes/moving_throat_pde_program_compact.md lines 1383-1455: matter i*hbar*d_t doublet, boxed Maxwell equation, boxed wall equation.
- absorber_reference: docs/boundary_and_noise_methods.md sections 0, 3, 5: Step-8 CAP decision and a generic CAP in the spirit of Manolopoulos 2002; this implementation does not use the wavelength-sized pole profile and claims no analytic reflection-free guarantee.
- maxwell_temporal_truncation: The driven Maxwell temporal block is an engineering-smoke truncation: only the diagonal -Z(w)*omega^2*A_N self-term is retained.  The omega-linear A0<->Ai temporal cross-couplings and gauge-term temporal pieces are dropped; the A0 coefficient is a placeholder, not the correct N=0 gauge-term reduction.
- compact_linearized_hierarchy: notes/moving_throat_pde_program_compact.md lines 1383-1455: linearized matter, Maxwell, and wall hierarchy with external f_ext.
- compact_open_reduced_lanes: notes/moving_throat_pde_program_compact.md lines 1377-1381: full coupled matter/gauge renormalization remains open.
- compact_confinement_coupling: notes/moving_throat_pde_program_compact.md lines 1063-1089: wall-to-matter confinement perturbation delta V_conf.
- compact_reachability_note: notes/moving_throat_pde_program_compact.md lines 2677-2710: conservative response bundle, cited only for reachability, not extracted.
- compact_wp3_card: notes/moving_throat_pde_program_compact.md lines 6846-6856: actual-branch card, deferred here.
- scope_decision: software/stage1_solver/decisions/02_step8c_s_eta_scope.md: surrogate-only 8c scope and open return-source decision.
- parent_status: docs/branch_realization_parent_status_decision.md: effective_closure keeps Path-A derivations deferred.

## Phasor Current

Frechet derivative of the Step-3 gauge-covariant number current at the WP1 background, applied to complex Step-8b perturbation phasors.
Algebra: delta j_i=(hbar/m)(delta psi_R grad_i psi_I0 + psi_R0 grad_i delta psi_I - delta psi_I grad_i psi_R0 - psi_I0 grad_i delta psi_R) -(q/m)(delta A_i rho0 + A_i0 delta rho), with delta rho=2(psi_R0 delta psi_R + psi_I0 delta psi_I). Charge current is q times this number current; Noether energy flux is S=mu*j.

| level | grid | sector | phasor_current_l2 | static_branch_current_l2 | non_null_vs_static_ratio |
| --- | --- | --- | --- | --- | --- |
| 0 | nr_4_nw_4 | number | 1.844666e-03 | 0.000000e+00 | static_exactly_null |
| 0 | nr_4_nw_4 | charge | 6.456332e-04 | 0.000000e+00 | static_exactly_null |
| 0 | nr_4_nw_4 | energy | 3.602605e-03 | 0.000000e+00 | static_exactly_null |
| 1 | nr_8_nw_8 | number | 1.986992e-03 | 0.000000e+00 | static_exactly_null |
| 1 | nr_8_nw_8 | charge | 6.954471e-04 | 0.000000e+00 | static_exactly_null |
| 1 | nr_8_nw_8 | energy | 4.088657e-03 | 0.000000e+00 | static_exactly_null |
| 2 | nr_16_nw_16 | number | 2.071903e-03 | 0.000000e+00 | static_exactly_null |
| 2 | nr_16_nw_16 | charge | 7.251662e-04 | 0.000000e+00 | static_exactly_null |
| 2 | nr_16_nw_16 | energy | 4.319045e-03 | 0.000000e+00 | static_exactly_null |
Static branch current is symmetry-null in this control; `static_exactly_null` replaces the old divide-by-floor ratio artifact. The soft control is the static null, while the load-bearing transport tooth is the driven non-null phasor_current_l2.

## Forced/CAP Balance

This table is a diagnostic decomposition, not a closed-conservation gate. Because `operator_exchange = -(number_divergence + static_projection)`, the `number_divergence` term cancels algebraically in `balance_number = number_divergence - source_number`. The residual that closes to the solve floor is the continuity projection of the converged driven LU residual, not an independent source-vs-divergence measurement. The per-term columns remain for inspection, and the FV closure row is an identity check.
| level | grid | sector | interior_source_l2 | interior_balance_l2 | interior_balance_l2_relative | observed_order | integrated_injected_abs | integrated_cap_absorbed_abs | integrated_harmonic_storage_abs | integrated_operator_exchange_abs | integrated_balance_residual_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | nr_16_nw_16 | number | 3.531355e-03 | 6.259553e-17 | 1.772564e-14 | -1.585002e+00 | 9.528088e-03 | 4.069381e-03 | 2.223726e-03 | 9.730392e-03 | 3.457785e-18 |
| 2 | nr_16_nw_16 | charge | 1.235974e-03 | 2.190844e-17 | 1.772564e-14 | -1.585002e+00 | 3.334831e-03 | 1.424283e-03 | 7.783041e-04 | 3.405637e-03 | 1.210225e-18 |
| 2 | nr_16_nw_16 | energy | 7.361386e-03 | 1.304853e-16 | 1.772564e-14 | -1.585002e+00 | 1.986205e-02 | 8.482943e-03 | 4.635531e-03 | 2.028376e-02 | 7.208023e-18 |

## Independent Gauss Residual Decrease

Independent center-gradient E=-grad(delta A0) reconstruction, reused from Step 6. At this CPU-bounded Step-8c resolution the relative residual decreases under refinement but does not close.
The relative residual is O(1)-O(8) on the finest grid because the enclosed phasor charge is near zero and the driven A0 lane is the Step-8b temporal-truncation engineering-smoke operator (`-Z omega^2 A_N` diagonal only). This is a decrease-rate diagnostic, not a closure claim.
| level | grid | surface | surface_flux_abs | enclosed_mu0_charge_abs | residual_abs | relative_residual | observed_order |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | nr_16_nw_16 | nested_surface_0 | 4.697278e-04 | 6.560826e-05 | 5.337324e-04 | 8.135140e+00 | 3.603464e+00 |
| 2 | nr_16_nw_16 | nested_surface_1 | 6.044218e-04 | 4.304361e-04 | 7.788215e-04 | 1.809378e+00 | 4.462891e+00 |

## Low-Frequency Surrogate Methodology

target-blind first-order Taylor fit in omega over a low-frequency band.
The band is bounded by omega<=0.5 to stay away from the Step-8b omega=6 resolution-limited regime; the choice is numerical, not target-driven.
The functionals were fixed before the runs and are not forcing overlaps: an interior linear-density mean, a CAP-free scalar-gauge norm, and a mid-wall eta mean. They are CAP/operator-methodology coefficients only.
Coefficient verdict `order_converges` means the extraction order-converges under refinement. The live 8c Gauss floor is kept honestly in `u_total`; when relative_uncertainty exceeds one, these rows are not precision measurements.
| functional_label | coefficient_label | finest_value | observed_order | finest_error_estimate | u_total | relative_uncertainty | precision_label | dominant_component | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| interior linear density mean, real phasor component | low-frequency Taylor constant term | 6.599336e-04 | 2.092378e+00 | 7.096733e-06 | 5.368657e-03 | 8.135147e+00 | not_precision_measured_resolution_limited | conservation | order_converges |
| interior linear density mean, real phasor component | low-frequency Taylor linear term | 1.701685e-05 | 3.057120e+00 | 9.123289e-08 | 1.384345e-04 | 8.135141e+00 | not_precision_measured_resolution_limited | conservation | order_converges |
| scalar gauge response L2 outside the CAP layer | low-frequency Taylor constant term | 1.430281e-03 | 1.400733e+00 | 1.233252e-05 | 1.163554e-02 | 8.135144e+00 | not_precision_measured_resolution_limited | conservation | order_converges |
| scalar gauge response L2 outside the CAP layer | low-frequency Taylor linear term | 6.714504e-05 | 1.804174e+00 | 4.994078e-07 | 5.462345e-04 | 8.135143e+00 | not_precision_measured_resolution_limited | conservation | order_converges |
| mid-wall eta mean, real phasor component | low-frequency Taylor constant term | 1.518587e-03 | 3.741275e+00 | 1.174082e-06 | 1.235391e-02 | 8.135140e+00 | not_precision_measured_resolution_limited | conservation | order_converges |
| mid-wall eta mean, real phasor component | low-frequency Taylor linear term | -2.076996e-05 | 6.669239e-01 | 2.263760e-05 | 1.704762e-04 | 8.207827e+00 | not_precision_measured_resolution_limited | conservation | order_converges |

Omega sampling stability on the finest level:
Count only the predeclared primary omega subset index 0 and require relative_change <= 0.35; absolute_change and u_total do not rescue the gate.
| functional | coefficient | omega_set_index | counts_for_stability | absolute_change | relative_change | stability_verdict |
| --- | --- | --- | --- | --- | --- | --- |
| interior_linear_density_mean_real | taylor_c0 | 0 | true | 7.541069e-07 | 1.142701e-03 | counted_stable |
| interior_linear_density_mean_real | taylor_c1 | 0 | true | 7.937967e-07 | 4.664768e-02 | counted_stable |
| interior_linear_density_mean_real | taylor_c0 | 1 | false | 1.309765e-06 | 1.984691e-03 | stress_diagnostic_within_relative_tol |
| interior_linear_density_mean_real | taylor_c1 | 1 | false | 1.190695e-05 | 6.997153e-01 | stress_diagnostic_exceeds_relative_tol |
| scalar_gauge_cap_free_l2 | taylor_c0 | 0 | true | 2.437950e-06 | 1.704525e-03 | counted_stable |
| scalar_gauge_cap_free_l2 | taylor_c1 | 0 | true | 2.566263e-06 | 3.821970e-02 | counted_stable |
| scalar_gauge_cap_free_l2 | taylor_c0 | 1 | false | 4.234334e-06 | 2.960491e-03 | stress_diagnostic_within_relative_tol |
| scalar_gauge_cap_free_l2 | taylor_c1 | 1 | false | 3.849395e-05 | 5.732955e-01 | stress_diagnostic_exceeds_relative_tol |
| wall_midband_eta_mean_real | taylor_c0 | 0 | true | 1.703978e-06 | 1.122081e-03 | counted_stable |
| wall_midband_eta_mean_real | taylor_c1 | 0 | true | 1.793661e-06 | 8.635841e-02 | counted_stable |
| wall_midband_eta_mean_real | taylor_c0 | 1 | false | 2.959540e-06 | 1.948878e-03 | stress_diagnostic_within_relative_tol |
| wall_midband_eta_mean_real | taylor_c1 | 1 | false | 2.690491e-05 | 1.295376e+00 | stress_diagnostic_exceeds_relative_tol |

Current-consistency gate decision: STOP-FLAGGED: no independent FV-current-vs-static-projection gate is counted. The live comparison was not a clean closing diagnostic on this truncated Step-8c CPU ladder, so it is not manufactured into pass_checks.

## WP3 Physical-Target Reachability

Physical d ln R_tr / R_target / epsilon_eta is blocked here by the open matter/gauge-to-wall source S_eta^(psi,A), cited above at compact lines 1383-1455 and 1377-1381. Per the parent-status decision and the 2026-06-14 consult record, that source is Path-A material and is deferred to the Step-9 verdict as blocked/deferred, not passed or failed.
Reachability split: R_norm, chi_Q, and N_Q=chi_Q^-1 are S_eta-independent source-map quantities and remain intact. R_pole, P_2, and P_4 are low-frequency-response-gated through D2/D4/N2/N4; Step 8c does not extract them and does not over-credit them as safe WP1 readouts.

## Counted Checks

- minimum_three_refinement_levels: PASS
- all_driven_solves_converged: PASS
- driven_phasor_current_non_null: PASS
- static_branch_current_null_contrast: PASS
- independent_gauss_residual_decreases_under_refinement: PASS
- low_frequency_samples_overdetermine_fit: PASS
- surrogate_coefficients_order_converge_under_refinement: PASS
- omega_sampling_stability: PASS

## Identity Checks

These are excluded from `passed` and carry `_not_a_physics_gate` suffixes.
- fv_divergence_theorem_closure_not_a_physics_gate: PASS - Finite-volume bookkeeping identity: net outward reconstructed current flux equals the volume integral of the same FV divergence.
- source_balance_closes_on_finest_grid_not_a_physics_gate: PASS - The source-balance residual is the continuity projection of the converged driven-solve residual. Since operator_exchange = -(number_divergence + static_projection), number_divergence cancels algebraically in balance_number = number_divergence - source_number; closing to the solve floor follows from driven LU convergence and is not an independent conservation measurement.
- source_balance_residual_decreases_not_a_physics_gate: PASS - The source-balance residual is the continuity projection of the converged driven-solve residual. Since operator_exchange = -(number_divergence + static_projection), number_divergence cancels algebraically in balance_number = number_divergence - source_number; closing to the solve floor follows from driven LU convergence and is not an independent conservation measurement.

## Asserted Checks

These are excluded from `passed` and carry `_not_a_physics_gate` suffixes.
- target_blind_surrogate_forcing_only_not_a_physics_gate: PASS - The RHS is p2_driven_surrogate_forcing; no reference packet or target value enters.
- physical_export_permitted_is_false_not_a_physics_gate: PASS - Step 8c writes only run manifests, a JSON table, and a report; it imports no firewalled model.
- matter_gauge_to_wall_source_deferred_not_a_physics_gate: PASS - The return source S_eta^(psi,A) remains open and is not supplied.
- maxwell_temporal_truncation_labelled_not_a_physics_gate: PASS - Retains only the diagonal -Z(w)*omega^2*A_N self-term.  It drops the omega-linear A0<->Ai temporal cross-couplings and the gauge-term temporal pieces; the A0 lane's -Z*omega^2 coefficient is an engineering-smoke placeholder, not the correct N=0 gauge-term reduction.  This temporal truncation is separate from the Step-8a spatial Ar scalarization and separate from the open S_eta^(psi,A) source.
- surrogate_response_coefficients_only_not_a_physics_gate: PASS - Fitted coefficients belong to predeclared CAP/operator-methodology functionals.

## Provenance

Machine-readable table: `software/stage1_solver/runs/step8c_conservation_response/step8c_conservation_response/step8c_table.json`.
Run with: `PYTHONPATH=software/stage1_solver/src timeout 600 python -m stage1_solver.step8c_harness`.

