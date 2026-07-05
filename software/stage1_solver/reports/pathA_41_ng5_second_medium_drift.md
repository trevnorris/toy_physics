SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu})

# pathA_41 NG5 SECOND_MEDIUM_DRIFT Gate

Primary verdict: `SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu})`.
Active irreducible set (computed): `['rho_B0', 'chi_c', 'C_hu']`.
Lineage finding: `NO_OVERCOUNT_ROUTE_A_PENDING`.
Engine agreement: `ENGINE_AGREE` over `201` compared scalar/audit quantities.

## Interpretation

- `ONE_CANDIDATE_MEDIUM_4D_TO_3D_REDUCTION_INCOMPLETE`.
- Physical meaning: The drift is not a separate substance; it is three unreduced 3D-brane-surface parameters: compression rho_B0, chi_c and embedding-mixing C_hu.
- Location closure: every row lives in `['4D bulk', '3D brane surface', 'throat/embedding seam']`; no fourth arena = `True`.
- Honest caveat: The brane is a postulated shear-supporting ordered phase; whether the one medium yields it and whether the three reductions close rather than no-go is genuinely open.

Named future reduction routes:
- `compression-sector 4D->3D reduction` targets `['rho_B0', 'chi_c']`: `DEFERRED_NOT_REGISTERED`.
- `embedding-overlap reduction` targets `['C_hu']`: `DEFERRED_NOT_REGISTERED`.

Reduction status:
- `rho_br`: `REGISTERED_PENDING(Route-A)`
- `mu_R`: `REGISTERED_PENDING(Route-A)`
- `rho_B0`: `NOT_YET_REGISTERED`
- `chi_c`: `NOT_YET_REGISTERED`
- `C_hu`: `NOT_YET_REGISTERED`

## Sim-Deferred And Calibrated Rows

- sim-deferred: `rho_br` -> `Route-A` -> re-adjudicate on solve
- sim-deferred: `mu_R` -> `Route-A` -> re-adjudicate on solve
- sim-deferred: `c_E` -> `Route-A` -> re-adjudicate on solve
- calibrated: `Q_E` -> `CALIBRATED_ANCHOR`
- calibrated: `ell` -> `CALIBRATED_GEOMETRY_INPUT`
- calibrated: `b` -> `CALIBRATED_GEOMETRY_INPUT`
- calibrated: `M_h` -> `CALIBRATED_GEOMETRY_INPUT`

## Origin Ledger

| p | incidence | route eval | origin | location |
|---|---|---|---|---|
| `rho` | base GNLS substrate | `` | `BASE_SUBSTRATE` | `4D bulk` |
| `K` | base GNLS EOS substrate | `` | `BASE_SUBSTRATE` | `4D bulk` |
| `m` | base GNLS constituent mass | `` | `BASE_SUBSTRATE` | `4D bulk` |
| `a` | T0 polar substrate length | `` | `BASE_SUBSTRATE` | `4D bulk` |
| `ell_g` | pathA_35 confinement width | `` | `ACCEPTED_GEOMETRY_SUBSTRATE` | `throat/embedding seam` |
| `g_ell(w)` | codim-1 confinement profile | `` | `ACCEPTED_PROFILE_GIVEN_ell_g` | `throat/embedding seam` |
| `varrho_br[rho]` | closed pathA_25 density-smectic layer inertia | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `Sigma_n[rho]` | closed pathA_25 density-smectic layer support | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `delta_Sigma[rho]` | closed pathA_25 density-smectic layer measure | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `rho_br` | active pathA_35 shear-surface brane inertia | `Route-A:True:VALID_REGISTERED_ROUTE` | `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED` | `3D brane surface` |
| `mu_R` | active pathA_35 shear-surface modulus | `Route-A:True:VALID_REGISTERED_ROUTE` | `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED` | `3D brane surface` |
| `c_E` | electric throat dynamic Green speed | `Route-A:True:VALID_REGISTERED_ROUTE` | `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED` | `throat/embedding seam` |
| `c_gamma` | light/shear speed | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `DEPENDENT` | `3D brane surface` |
| `rho_B0` | C5 compression density amplitude | `Future-Compression-4D-to-3D:False:REJECTED_RESULT_STATUS_PROMISSORY_ONLY` | `IRREDUCIBLY_INDEPENDENT` | `3D brane surface` |
| `chi_c` | C5 compression susceptibility | `Future-Compression-4D-to-3D:False:REJECTED_RESULT_STATUS_PROMISSORY_ONLY` | `IRREDUCIBLY_INDEPENDENT` | `3D brane surface` |
| `B_eff` | derived C5 density modulus | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `DEPENDENT_ON_IRREDUCIBLE` | `3D brane surface` |
| `C_hu` | embedding h/u_L mixing overlap | `Future-Embedding-Overlap:False:REJECTED_RESULT_STATUS_PROMISSORY_ONLY` | `IRREDUCIBLY_INDEPENDENT` | `throat/embedding seam` |
| `Q_E` | electric throat source magnitude | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `CALIBRATED_ANCHOR` | `throat/embedding seam` |
| `ell` | throat/wall profile scale | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `CALIBRATED_GEOMETRY_INPUT` | `throat/embedding seam` |
| `b` | compact throat source half-separation/form factor scale | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `CALIBRATED_GEOMETRY_INPUT` | `throat/embedding seam` |
| `M_h` | h-sector zero-mode normalization/mass coefficient | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `CALIBRATED_GEOMETRY_INPUT` | `throat/embedding seam` |
| `K_h` | h-sector stiffness | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `DEPENDENT` | `throat/embedding seam` |
| `q_h` | electric throat source projection | `None:False:NO_REGISTERED_ROUTE_FOR_PARAM` | `DEPENDENT` | `throat/embedding seam` |
| `c_L1` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `c_L2` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `A_R` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `k_R` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `lambda_Cdiv` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `chi_Cpin` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `J_Pu` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `kappa_Pu` | pathA_25 density-smectic driver | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `lambda_Pu` | pathA_35 parity-repaired P-u coupling | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `Omega_w` | pathA_35 bare u_w gap scale | `` | `OUT_OF_ACTIVE_NG5` | `3D brane surface` |
| `lambda_N` | pathA_38 wall-internal potential coefficient | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |
| `lambda_tau` | pathA_38 wall-internal tau mass coefficient | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |
| `Nu` | moving-source normalization | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |
| `a_T` | transverse moving-source amplitude | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |
| `a_Tp` | second transverse moving-source amplitude | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |
| `a_L` | longitudinal moving-source amplitude | `` | `OUT_OF_ACTIVE_NG5` | `throat/embedding seam` |

## Controls

| control | fired | transition |
|---|---:|---|
| `AB_delete_registry` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu}) |
| `route_blank_Route_A` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu}) |
| `route_field_blank_Route_A_target_blind` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu}) |
| `route_field_blank_Route_A_missing_objects` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu}) |
| `calibration_ablation_Q_E` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,Q_E}) |
| `irreducible_synthetic` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,xi_active}) |
| `reducible_derived_synthetic` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) |
| `contradiction` | `True` | SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu}) -> NO_GO(cone-lock-feedback) |
| `residual_multiplier_ablation` | `True` | lineage SAME when residual=1; lineage DIFFERENT when residual=Xi_residual |
| `route_eval_recorded_for_all_active_rows` | `True` | all active production rows carry RouteEvaluation records |

## Dual-Engine Split

- SymPy and Mathematica independently compute MLT dimension closure, residual lineage states, RouteEvaluation validity, origin classification, active irreducibles, and control transitions.
- The contradiction control is gated on a recomputed pathA_40 `freedom_tie` UNSAT result, not on a typed no-go flag.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_41_ng5_second_medium_drift_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_41_ng5_second_medium_drift.wl
timeout 600 python3 software/stage1_solver/tools/pathA_41_ng5_second_medium_drift_sympy.py --compare
```

Controls fired: `AB_delete_registry, route_blank_Route_A, route_field_blank_Route_A_target_blind, route_field_blank_Route_A_missing_objects, calibration_ablation_Q_E, irreducible_synthetic, reducible_derived_synthetic, contradiction, residual_multiplier_ablation, route_eval_recorded_for_all_active_rows`.