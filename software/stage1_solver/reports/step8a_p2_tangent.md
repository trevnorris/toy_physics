# Step 8a P2 Tangent Operator

Overall engineering gate: PASS
Config hash: `a907c392e617de95`
Diagnostics digest: `2d11fb6cf245372c`

**Scope:** engineering-smoke, target-blind, static omega=0 grouped-real P2 tangent only. The observables below are raw-field numerical surrogates, not physical WP3 observables; there is no response-packet export and no extraction map.

## Transliteration Sources

- linearized_skeleton: notes/moving_throat_pde_program_compact.md lines 1383-1455
- wall_modal_split: notes/moving_throat_pde_program_compact.md lines 1298-1308
- p2_degeneracy: notes/moving_throat_pde_program_compact.md lines 1340-1364
- open_backreaction: notes/moving_throat_pde_program_compact.md lines 1377-1381
- wall_to_matter: notes/moving_throat_pde_program_compact.md lines 1075-1089
- bulk_angular_reduction: Step-8a directive bulk angular reduction paragraph
- maxwell_scalarization: Step-8a engineering-smoke retained-lane scalarization; the full vector-harmonic Maxwell l=2 reduction is deferred.

## Pinned Vs Open

Pinned in this build: centrifugal matter/scalar angular reduction, live `[K_eta+6T_Omega]` wall term, and the Step-3 ratio-confinement wall-to-matter `delta V_conf` linearized with respect to the displaced radius `R=R0+eta`.

Engineering-smoke retained Maxwell angular block: exact for scalar lanes `A0` and `Aw`; componentwise scalarization for the radial vector lane `Ar`, omitting the `-2 A_r/r^2` vector-Laplacian self-term and the truncated transverse vector-harmonic lanes. The full vector-harmonic `l=2` Maxwell reduction is a standard derivable downstream task and is deferred.

Genuinely open and deferred: the physical matter/gauge-to-wall `S_eta^(psi,A)` source in compact lines 1377-1381; Step 8a uses a target-blind surrogate `f_ext`.

## MMS For New Operator Pieces

### mms_p2_centrifugal_tensor_laplacian

P2 scalar/matter bulk angular reduction: tensor Laplacian plus -l(l+1)u/r^2 with l=2.
Source: moving_throat_pde_program_compact.md lines 1298-1308 and Step-8a bulk angular reduction paragraph.
Forcing: SymPy evaluated d_r^2 u + 2 r^-1 d_r u + d_w^2 u - l(l+1)u/r^2 for a regular r^2 manufactured P2 coefficient.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 3.852349e-02 | - | 3.831033e+00 |
| nr_48_nw_40 | 4.166667e-02 | 9.925992e-03 | 1.956455e+00 | 3.830051e+00 |
| nr_96_nw_80 | 2.083333e-02 | 2.500718e-03 | 1.988869e+00 | 3.829797e+00 |
| nr_192_nw_160 | 1.041667e-02 | 6.264030e-04 | 1.997179e+00 | 3.829733e+00 |

Checks:
- observed_order: PASS
- final_error: PASS

### mms_p2_localized_maxwell_angular

Localized Maxwell retained-lane l=2 engineering-smoke scalarization: existing H=Z radial/w operator plus componentwise l(l+1)Z(w)A_N/r^2.
Source: Step-8a retained-lane scalarization. Exact for scalar lanes A0 and Aw; Ar is not the full vector-harmonic l=2 reduction.
Forcing: SymPy evaluated the continuum H=Z Maxwell operator and then added the scalarized l(l+1)Z(w)A_N/r^2 angular term for each retained radial coefficient.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 5.239697e-02 | - | 1.105269e+00 |
| nr_48_nw_40 | 4.166667e-02 | 1.345334e-02 | 1.961519e+00 | 1.104459e+00 |
| nr_96_nw_80 | 2.083333e-02 | 3.135771e-03 | 2.101072e+00 | 1.104256e+00 |
| nr_192_nw_160 | 1.041667e-02 | 7.339654e-04 | 2.095036e+00 | 1.104205e+00 |

Checks:
- observed_order: PASS
- final_error: PASS

## Tangent Jacobian Gate

Gate A central-FD check: relative=1.760876e-12, absolute=7.672023e-11, epsilon=1.000000e-05.

Static l=2 well-posedness on the small grid: sigma_min=5.040231e-01, condition=3.113853e+02.

## Static Tangent Solves

| level | grid | dof | background_residual | iterations | final_residual_norm | static_residual_linf | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | p2_static_l0_nr_4_nw_4 | 84 | 1.280842e-08 | 1 | 1.214306e-17 | 1.214306e-17 | true |
| 1 | p2_static_l1_nr_8_nw_8 | 328 | 1.014022e-08 | 1 | 4.510281e-17 | 4.510281e-17 | true |
| 2 | p2_static_l2_nr_16_nw_16 | 1296 | 9.513920e-09 | 1 | 2.428613e-16 | 2.428613e-16 | true |
| 3 | p2_static_l3_nr_32_nw_32 | 5152 | 9.338841e-09 | 1 | 2.121567e-15 | 2.121567e-15 | true |

## Surrogate Observable Convergence

These are target-blind raw-field functionals of `delta_u`; they are not physical P2 response observables.

| label | finest_grid | finest_value | last_observed_order | richardson_estimate | finest_error_estimate | floor_or_null_label |
| --- | --- | --- | --- | --- | --- | --- |
| total raw response L2 | p2_static_l3_nr_32_nw_32 | 1.583991e-02 | 2.491622e+00 | 1.584067e-02 | 7.576916e-07 | expected-order convergence |
| matter perturbation interior L2 | p2_static_l3_nr_32_nw_32 | 1.065270e-02 | 2.146991e+00 | 1.065483e-02 | 2.135729e-06 | expected-order convergence |
| A0 perturbation L2 | p2_static_l3_nr_32_nw_32 | 4.499283e-03 | 1.983393e+00 | 4.498506e-03 | 7.775965e-07 | expected-order convergence |
| spatial gauge perturbation L2 | p2_static_l3_nr_32_nw_32 | 5.078055e-03 | 2.235645e+00 | 5.082442e-03 | 4.386445e-06 | expected-order convergence |
| wall eta L2 | p2_static_l3_nr_32_nw_32 | 3.046083e-03 | 2.005527e+00 | 3.045910e-03 | 1.736092e-07 | expected-order convergence |
| lower-wall eta absolute value | p2_static_l3_nr_32_nw_32 | 2.327929e-03 | 1.986551e+00 | 2.327414e-03 | 5.146302e-07 | expected-order convergence |
| static tangent residual Linf | p2_static_l3_nr_32_nw_32 | 2.121567e-15 | -3.247928e+00 | 2.121567e-15 | 0.000000e+00 | null diagnostic |

## Counted Checks

- tangent_matches_fd_jacobian: PASS
- p2_centrifugal_mms_order: PASS
- p2_maxwell_angular_mms_order: PASS
- static_tangent_solves_converged: PASS
- p2_static_operator_wellposed: PASS
- surrogate_observable_orders: PASS

## Asserted Checks

These are excluded from `passed` and carry `_not_a_physics_gate` suffixes.
- target_blind_surrogate_forcing_only_not_a_physics_gate: PASS - The Step-8a RHS is generated by p2_surrogate_forcing and reads no target, reference packet, or extraction map.
- physical_export_permitted_is_false_not_a_physics_gate: PASS - Step 8a emits no physical packet and does not import the firewalled physical model.
- matter_gauge_to_wall_source_deferred_not_a_physics_gate: PASS - Compact lines 1377-1381 leave S_eta^(psi,A) open; Step 8a does not invent it.
- wall_constitutive_placeholders_labelled_not_a_physics_gate: PASS - step-8a engineering-smoke placeholders only; target-blind; not a physical response packet

## Provenance And Reproduction

Machine-readable table: `software/stage1_solver/runs/step8a_p2_tangent/step8a_p2_tangent_operator/p2_tangent_table.json`.
Run with: `PYTHONPATH=software/stage1_solver/src timeout 600 python -m stage1_solver.p2_tangent_harness`.
Forward note for Step 8b: if the absorber becomes a differential stencil, it must be included in the MMS continuum forcing before certification.

