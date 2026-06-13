# Step 6 Conservation Diagnostics

Overall engineering gate: PASS
Config hash: `684c53f5562074e4`
Diagnostics digest: `9752e5eb227f21f1`

**Stationary framing:** no time-marching drift loop is run. Conservation is reported as conservative finite-volume bulk divergence plus stationary budget closure, with explicit sponge accounting.

**Discipline:** engineering smoke, placeholder parameters, target-blind, no field-to-coefficient export, and no physical packet.

## Method

Levels: `[(6, 4), (12, 8), (24, 16)]` with refinement ratio `2`.
Interior window: `{'r_min': 0.0, 'r_max': 1.2, 'w_min': 0.2, 'w_max': 1.0}`.
Number current: `stage1_solver.coupled_branch._matter_number_current`.
Divergence: area-weighted FV face fluxes from operators.tensor_vector_face_fluxes, operators.tensor_weighted_gradient_fluxes, and operators.tensor_flux_divergence.
Gauss reconstruction: independent center-gradient reconstruction of E=-grad(A0), averaged to closed-surface faces and multiplied by Z; differs from the solver's tensor_weighted_gradient_fluxes operator.
Measures: grid.cell_volumes, grid.radial_face_areas, grid.radial_shell_volumes.
Sponge accounting: boundary_sponge_profile_torch and sponge_*_strength config fields.
EOS: U=(K/4)*rho^5; h=(5K/4)*rho^4 remains in physics.quintic_enthalpy.
Virial/Pohozaev identity: omitted; no independent term-by-term parent-action derivation included.

## Solve Rows

| mode | level | grid | dof | final_residual_linf | converged | message |
| --- | --- | --- | --- | --- | --- | --- |
| sponge_off | 0 | nr_6_nw_4 | 121 | 1.068892e-08 | true | continuation completed |
| sponge_off | 1 | nr_12_nw_8 | 481 | 9.691578e-09 | true | continuation completed |
| sponge_off | 2 | nr_24_nw_16 | 1921 | 9.375697e-09 | true | continuation completed |
| sponge_on | 0 | nr_6_nw_4 | 121 | 1.557096e-11 | true | continuation completed |
| sponge_on | 1 | nr_12_nw_8 | 481 | 3.353096e-12 | true | continuation completed |
| sponge_on | 2 | nr_24_nw_16 | 1921 | 3.836931e-13 | true | continuation completed |

## Local Conservation Residuals

The mass, charge-current, and energy-flux sectors are null diagnostics on this isotropic stationary branch because the spatial current and spatial gauge field vanish by symmetry. The current-carrying conservation test is deferred to the WP3 tangent in step 8.

| mode | grid | sector | interior_l2_relative | interior_linf_relative | transport_signal_l2 | label | observed_order |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_on | nr_24_nw_16 | number | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | null diagnostic | - |
| sponge_on | nr_24_nw_16 | charge | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | null diagnostic | - |
| sponge_on | nr_24_nw_16 | energy | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | null diagnostic | - |

## Gauss-Law Closure

This is an independent reconstruction check of the localized Gauss law: the closed-surface flux uses center-gradient `E=-grad(A0)` values averaged to faces, not the solver's weighted-gradient face flux. The comparison is `surface flux of Z F^{i0}` versus `mu0*q*int rho dV` on nested interior control volumes.

| mode | level | grid | surface | surface_flux | enclosed_mu0_charge | absolute_residual | relative_residual | observed_order | observed_order_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_off | 0 | nr_6_nw_4 | nested_surface_0 | 2.773652e-02 | 3.198151e-02 | -4.244991e-03 | 1.327326e-01 | - | - |
| sponge_off | 0 | nr_6_nw_4 | nested_surface_1 | 6.214916e-02 | 6.895434e-02 | -6.805181e-03 | 9.869112e-02 | - | - |
| sponge_off | 1 | nr_12_nw_8 | nested_surface_0 | 3.092049e-02 | 3.067569e-02 | 2.447935e-04 | 7.980047e-03 | - | - |
| sponge_off | 1 | nr_12_nw_8 | nested_surface_1 | 6.818302e-02 | 6.688398e-02 | 1.299039e-03 | 1.942227e-02 | - | - |
| sponge_off | 2 | nr_24_nw_16 | nested_surface_0 | 3.041634e-02 | 3.035031e-02 | 6.602921e-05 | 2.175569e-03 | 4.425760e+00 | - |
| sponge_off | 2 | nr_24_nw_16 | nested_surface_1 | 6.671528e-02 | 6.637924e-02 | 3.360444e-04 | 5.062493e-03 | 2.464720e+00 | - |
| sponge_on | 0 | nr_6_nw_4 | nested_surface_0 | 3.331778e-02 | 3.968641e-02 | -6.368628e-03 | 1.604738e-01 | - | - |
| sponge_on | 0 | nr_6_nw_4 | nested_surface_1 | 7.416819e-02 | 8.546095e-02 | -1.129277e-02 | 1.321395e-01 | - | - |
| sponge_on | 1 | nr_12_nw_8 | nested_surface_0 | 3.920973e-02 | 3.841801e-02 | 7.917231e-04 | 2.060812e-02 | - | - |
| sponge_on | 1 | nr_12_nw_8 | nested_surface_1 | 8.598800e-02 | 8.335046e-02 | 2.637539e-03 | 3.164396e-02 | - | - |
| sponge_on | 2 | nr_24_nw_16 | nested_surface_0 | 3.829671e-02 | 3.808895e-02 | 2.077542e-04 | 5.454446e-03 | 3.206302e+00 | - |
| sponge_on | 2 | nr_24_nw_16 | nested_surface_1 | 8.350721e-02 | 8.282492e-02 | 6.822897e-04 | 8.237735e-03 | 2.102167e+00 | - |

Operator-flux Maxwell residual closure:

This second table is not independent: it restates the discrete Gauss law as the integrated Maxwell-a0 Newton residual over the same control volumes. It is retained only as a solver-floor check.

| mode | level | grid | surface | absolute_residual | relative_residual | observed_order | observed_order_note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_on | 0 | nr_6_nw_4 | nested_surface_0 | -7.881196e-14 | 1.985868e-12 | - | solver-floor diagnostic |
| sponge_on | 0 | nr_6_nw_4 | nested_surface_1 | -4.801715e-15 | 5.618606e-14 | - | solver-floor diagnostic |
| sponge_on | 1 | nr_12_nw_8 | nested_surface_0 | -7.459311e-15 | 1.941618e-13 | - | solver-floor diagnostic |
| sponge_on | 1 | nr_12_nw_8 | nested_surface_1 | -9.617307e-15 | 1.153840e-13 | - | solver-floor diagnostic |
| sponge_on | 2 | nr_24_nw_16 | nested_surface_0 | -1.318390e-16 | 3.461344e-15 | - | solver-floor diagnostic |
| sponge_on | 2 | nr_24_nw_16 | nested_surface_1 | -8.326673e-17 | 1.005334e-15 | - | solver-floor diagnostic |

## Budgets

Budget closure here is a finite-volume consistency identity, not an independent physical conservation balance: `net outward flux = interior sponge absorbed + local balance residual` is exact to roundoff because the same conservative FV face fluxes generate both the divergence and the surface term. The interior sponge term is structurally zero because the configured sponge support lies outside the interior window; the real sponge measurement is the outer-layer absorbed amount in the Sponge Accounting table.

| mode | grid | sector | net_outward_flux | interior_local_balance_residual | closure_absolute | closure_relative | observed_order |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_off | nr_24_nw_16 | number | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |
| sponge_off | nr_24_nw_16 | charge | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |
| sponge_off | nr_24_nw_16 | energy | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |
| sponge_on | nr_24_nw_16 | number | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |
| sponge_on | nr_24_nw_16 | charge | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |
| sponge_on | nr_24_nw_16 | energy | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | - |

## Sponge Accounting

| mode | grid | sector | interior_absorbed | outer_layer_absorbed | total_absorbed | interior_zero |
| --- | --- | --- | --- | --- | --- | --- |
| sponge_off | nr_24_nw_16 | number | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |
| sponge_off | nr_24_nw_16 | charge | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |
| sponge_off | nr_24_nw_16 | energy | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |
| sponge_on | nr_24_nw_16 | number | 0.000000e+00 | 8.464724e-02 | 8.464724e-02 | true |
| sponge_on | nr_24_nw_16 | charge | 0.000000e+00 | 2.962654e-02 | 2.962654e-02 | true |
| sponge_on | nr_24_nw_16 | energy | 0.000000e+00 | 8.468377e-02 | 8.468377e-02 | true |

The sponge-on/off split isolates the configured sponge terms. The interior absorbed amount is zero in every sector, which is the conservation-side restatement of Step 5's interior-zero sponge property.

## Pass Checks

- minimum_three_levels: PASS
- all_solves_converged: PASS
- two_nested_gauss_surfaces: PASS
- independent_gauss_decreases_under_refinement: PASS
- independent_gauss_observed_order_at_least_one: PASS
- operator_maxwell_residual_positive_below_solver_floor: PASS
- sponge_on_absorbs_outer_layer: PASS
- null_current_sectors_labelled: PASS

Identity checks, reported but not counted as physics gates:
- fv_budget_identity_roundoff_not_a_physics_gate: PASS
- sponge_support_excludes_interior_window_not_a_physics_gate: PASS

## Limitations

Mass/charge spatial-current residuals and the energy-flux divergence are null/floor diagnostics here, not current-carrying conservation tests. The isotropic real-amplitude branch has no phase winding, no transverse vector potential, and therefore no spatial transport to conserve. The current-carrying test is intentionally deferred to step 8.

## Reproduction

```bash
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.conservation_harness
timeout 600 env PYTHONPATH=software/stage1_solver/src pytest software/stage1_solver/tests/test_stage1_solver.py
```

Machine-readable table: `software/stage1_solver/runs/step6_conservation_diagnostics/step6_conservation_diagnostics/conservation_diagnostics_table.json`.
Target-blindness: no benchmark targets, extraction coefficients, or export packet paths are imported by the Step 6 diagnostics.
Export guard: the external physical export guard remains untouched; no physical packet is emitted.

