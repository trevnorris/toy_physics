# Step 8b Driven P2 Tangent And Absorber

Overall engineering gate: PASS_WITH_RESOLUTION_LIMIT
Config hash: `865f59abd2372007`
Diagnostics digest: `0e63e7a1a935231f`

**Scope:** engineering-smoke, target-blind driven tangent. Raw magnitudes are diagnostics only; no physical packet, no extraction map, and no fitted response coefficients.

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

## Driven Operator

- representation: complex128 packed P2 fields with the same six lanes as Step 8a; the real Step-8a Jacobian is complexified without changing its entries.
- matter_time_term: The +/-hbar*omega BdG doublet is represented in real/imag lanes as deltaR_re=+i*hbar*omega*psi_imag and deltaR_im=-i*hbar*omega*psi_real.
- maxwell_time_term: Retains only the diagonal -Z(w)*omega^2*A_N self-term.  It drops the omega-linear A0<->Ai temporal cross-couplings and the gauge-term temporal pieces; the A0 lane's -Z*omega^2 coefficient is an engineering-smoke placeholder, not the correct N=0 gauge-term reduction.  This temporal truncation is separate from the Step-8a spatial Ar scalarization and separate from the open S_eta^(psi,A) source.
- wall_time_term: -mu_eta*omega^2*eta.
- absorber: Generic smooth additive -i*sigma(w) upper-exit CAP using a normalized polynomial/rational ramp inspired by the Manolopoulos rational form, but not the wavelength-sized pole profile, not the Manolopoulos transmission-free profile, and not assigned an analytic reflection bound; reflection is only the measured domain-extension contamination metric.

The omega->0, CAP-off sanity check gives relative=0.000000e+00. This is a weak static-path consistency check that is zero by construction when omega=0 and CAP is off; it is not an omega/CAP physics certification.
Matter sign dispersion sanity check: hbar*omega=1.233701e+00, E(k)=1.233701e+00, relative=0.000000e+00; the old negated sign would give relative=2.000000e+00.

## Frequencies

Drive frequencies are target-blind numerical probes selected from near-static, intermediate, and wall-propagating regimes.
| omega | wall_effective_k | propagating_wall_lane | total_response_l2 | wall_eta_l2 |
| --- | --- | --- | --- | --- |
| 5.000000e-02 | 0.000000e+00 | false | 1.350151e-02 | 2.385953e-03 |
| 1.500000e+00 | 0.000000e+00 | false | 1.926428e-02 | 3.065973e-03 |
| 6.000000e+00 | 5.126032e+00 | true | 1.898246e-02 | 5.073129e-04 |

## MMS

### mms_p2_driven_frequency_and_cap_matter

Matter BdG time term in real/imag lanes plus a generic additive CAP; the differential carrier is the P2 scalar operator with a nonzero anomalous h'(rho_bg)*psi_bg^2*delta_psi* block.
Source: compact lines 1406-1418 and boundary methods sections 3, 5.
Forcing: SymPy applied i*hbar*d_t to the manufactured doublet (delta_psi,delta_psi*)=(x+i*y, x-i*y) with exp(-i*omega*t), then transformed back to real/imag lanes.  The exact carrier also includes an independently evaluated anomalous g*delta_psi* block and the additive -i*sigma CAP coefficient.
| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 2.053793e-02 | - | 2.449827e+00 |
| nr_48_nw_40 | 4.166667e-02 | 5.291820e-03 | 1.956455e+00 | 2.449501e+00 |
| nr_96_nw_80 | 2.083333e-02 | 1.333202e-03 | 1.988869e+00 | 2.449419e+00 |

### mms_p2_driven_frequency_and_cap_maxwell

Temporal Maxwell engineering-smoke diagonal -Z(w)*omega^2 retained-lane self-term plus additive CAP.  This temporal truncation is additional to the Step-8a spatial Ar scalarization caveat.
Source: compact lines 1423-1434 and boundary methods sections 3, 5.
Forcing: SymPy evaluated d_t^2 exp(-i*omega*t)=-omega^2 and applied only the retained diagonal Z(w)*d_t^2 A_N self-term plus -i*sigma*A_N.  The omega-linear A0<->Ai terms and gauge temporal pieces are explicitly not included.
| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 5.239697e-02 | - | 1.012412e+00 |
| nr_48_nw_40 | 4.166667e-02 | 1.345334e-02 | 1.961519e+00 | 1.011684e+00 |
| nr_96_nw_80 | 2.083333e-02 | 3.135771e-03 | 2.101072e+00 | 1.011501e+00 |

### mms_p2_driven_frequency_and_cap_wall

Wall mu_eta*d_t^2 term and additive CAP on the modal wall operator.
Source: compact lines 1441-1451 and boundary methods sections 3, 5.
Forcing: SymPy evaluated d_t^2 exp(-i*omega*t)=-omega^2, giving the mu_eta*d_t^2 eta term as -mu_eta*omega^2*eta, plus -i*sigma*eta.
| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nw_32 | 5.000000e-02 | 1.829135e-02 | - | 5.624587e+00 |
| nw_64 | 2.500000e-02 | 4.627190e-03 | 1.982953e+00 | 5.625632e+00 |
| nw_128 | 1.250000e-02 | 1.160276e-03 | 1.995668e+00 | 5.625890e+00 |

## Gate A And Well-Posedness

Static background JVP gate-A re-check: relative=1.757570e-12, absolute=7.657467e-11. This is internal consistency only; omega/CAP physics is certified by MMS and fidelity review.
Driven assembled operator FD check: relative=1.723284e-16, absolute=1.427194e-14.
Complex small-grid singular values: sigma_min=1.106258e+00, condition=1.426338e+02. The configured wellposedness floor 1.0e-08 is a loose numerical floor, not a sharp spectral margin.

## Driven Solve Convergence

Resolution-limited fallback is active: the clean all-non-floor-observable gate at bar 1.450 FAILS.  The scoped convergence claim uses the predeclared bar 1.350; this is a reduced engineering claim, not a clean all-observable convergence pass.
Clean-gate failing observables: total_response_l2, matter_interior_response_l2, spatial_gauge_response_l2, wall_eta_l2.
Scoped gated observables: total_response_l2, scalar_gauge_response_l2, wall_lower_endpoint_eta_abs.
Scoped exclusions reported as under-resolved/drifting: matter_interior_response_l2, spatial_gauge_response_l2, wall_eta_l2.
Scoped rationale: Fallback for the omega=6 resolution-limited ladder: total response plus the scalar-gauge and lower-wall diagnostics are counted; lane diagnostics that still show under-resolved drift are reported but not claimed converged.

Wall-lane wavelength resolution for omega=6 shows why the individual lane diagnostics are resolution-limited on this CPU-bounded ladder.
| grid | nw | dw | wall_effective_k | wavelength | points_per_wavelength |
| --- | --- | --- | --- | --- | --- |
| p2_driven_l0_nr_4_nw_4 | 4 | 4.000000e-01 | 5.126032e+00 | 1.225741e+00 | 3.064351e+00 |
| p2_driven_l1_nr_8_nw_8 | 8 | 2.000000e-01 | 4.986188e+00 | 1.260118e+00 | 6.300590e+00 |
| p2_driven_l2_nr_16_nw_16 | 16 | 1.000000e-01 | 4.980272e+00 | 1.261615e+00 | 1.261615e+01 |

| level | grid | omega | dof | background_residual | driven_residual_linf | converged |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | p2_driven_l0_nr_4_nw_4 | 6.000000e+00 | 84 | 3.193170e-10 | 1.952362e-17 | true |
| 1 | p2_driven_l1_nr_8_nw_8 | 6.000000e+00 | 328 | 8.157156e-13 | 1.212598e-16 | true |
| 2 | p2_driven_l2_nr_16_nw_16 | 6.000000e+00 | 1296 | 1.691574e-12 | 5.103744e-16 | true |

| label | finest_grid | finest_value | last_observed_order | richardson_estimate | finest_error_estimate | floor_or_null_label |
| --- | --- | --- | --- | --- | --- | --- |
| total raw response L2 | p2_driven_l2_nr_16_nw_16 | 1.979944e-02 | 1.439501e+00 | 1.952078e-02 | 2.786528e-04 | reduced-order convergence |
| matter perturbation interior L2 | p2_driven_l2_nr_16_nw_16 | 1.895065e-02 | 3.629876e-01 | 1.691394e-02 | 2.036703e-03 | reduced-order convergence |
| A0 perturbation L2 | p2_driven_l2_nr_16_nw_16 | 7.850697e-04 | 3.637905e+00 | 7.832874e-04 | 1.782249e-06 | expected-order convergence |
| spatial gauge perturbation L2 | p2_driven_l2_nr_16_nw_16 | 6.349520e-04 | 5.454744e-01 | 6.654439e-04 | 3.049187e-05 | reduced-order convergence |
| wall eta L2 | p2_driven_l2_nr_16_nw_16 | 4.880059e-04 | -7.122167e-01 | - | - | drifts |
| lower-wall eta absolute value | p2_driven_l2_nr_16_nw_16 | 3.147658e-04 | 2.432552e+00 | 3.338700e-04 | 1.910428e-05 | expected-order convergence |
| driven_residual_linf | p2_driven_l2_nr_16_nw_16 | 5.103744e-16 | -1.935362e+00 | 5.103744e-16 | 0.000000e+00 | null diagnostic |

## Reflection

Domain-extension contamination metric on a standalone reduced 1D wall_s_eta_operator proxy for the retained wall lane: solve the proxy wall frequency operator on a base domain and on an aligned extended domain, then measure the fixed-window weighted L2 difference divided by the extended-domain signal.  This does not characterize matter, Maxwell, or full six-lane coupled-operator reflection.
The coordinate window, source, frequency, grid spacing, wall coefficients, and norm are fixed before absorber on/off is compared; the reflecting control uses the identical metric with sigma=0.
Generic smooth additive -i*sigma(w) upper-exit CAP using a normalized polynomial/rational ramp inspired by the Manolopoulos rational form, but not the wavelength-sized pole profile, not the Manolopoulos transmission-free profile, and not assigned an analytic reflection bound; reflection is only the measured domain-extension contamination metric.
The reflection floor multiplier is 8.0 times the Step-4 raw-field floor; it is a conditioning-motivated free choice, not a derived bound.
| case | reflection_coefficient | interior_l2_change | interior_signal_l2_reference | effective_k | fit_cells |
| --- | --- | --- | --- | --- | --- |
| absorbed | 1.155166e-03 | 1.080063e-05 | 9.349853e-03 | 5.128776e+00 | 14 |
| reflecting_control | 8.131772e-01 | 2.923444e-02 | 3.595089e-02 | 5.128776e+00 | 14 |
Target-blind floor: 2.782001e-03; control contrast: 7.039485e+02.

## Counted Checks

- static_background_jacobian_gate_a_reverified: PASS
- driven_operator_internal_fd_check: PASS
- omega_zero_recovers_static_operator: PASS
- matter_time_sign_dispersion_sanity: PASS
- new_frequency_and_cap_terms_mms_order: PASS
- driven_tangent_solves_converged: PASS
- complex_driven_operator_wellposed: PASS
- scoped_surrogate_observable_convergence: PASS
- reflection_metric_below_target_blind_floor: PASS
- reflecting_control_has_teeth: PASS
- propagating_drive_present: PASS

## Asserted Checks

These are excluded from `passed` and carry `_not_a_physics_gate` suffixes.
- target_blind_surrogate_forcing_only_not_a_physics_gate: PASS - The driven RHS extends p2_surrogate_forcing with smooth numerical source data.
- physical_export_permitted_is_false_not_a_physics_gate: PASS - Step 8b writes only run manifests/tables and a report; it does not import the firewalled model.
- matter_gauge_to_wall_source_deferred_not_a_physics_gate: PASS - The compact S_eta^(psi,A) source remains open and is not invented here.
- maxwell_radial_lane_scalarization_labelled_not_a_physics_gate: PASS - The Step-8a spatial Maxwell Ar lane remains componentwise and is not the full vector-harmonic l=2 reduction.
- maxwell_temporal_truncation_labelled_not_a_physics_gate: PASS - Retains only the diagonal -Z(w)*omega^2*A_N self-term.  It drops the omega-linear A0<->Ai temporal cross-couplings and the gauge-term temporal pieces; the A0 lane's -Z*omega^2 coefficient is an engineering-smoke placeholder, not the correct N=0 gauge-term reduction.  This temporal truncation is separate from the Step-8a spatial Ar scalarization and separate from the open S_eta^(psi,A) source.
- raw_response_vs_omega_not_fitted_not_a_physics_gate: PASS - The omega table reports raw magnitudes only and performs no coefficient extraction.

## Provenance

Machine-readable table: `software/stage1_solver/runs/step8b_driven_absorber/step8b_driven_p2_absorber/p2_driven_table.json`.
Run with: `PYTHONPATH=software/stage1_solver/src timeout 600 python -m stage1_solver.step8b_harness`.
Forward note for Step 8c: the OPEN matter/gauge-to-wall source remains deferred, and physical extraction/export gates remain untouched.

