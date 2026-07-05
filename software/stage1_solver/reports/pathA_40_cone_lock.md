CONE_LOCK_CALIBRATED

# pathA_40 Cone-Lock Gate

Primary verdict: `CONE_LOCK_CALIBRATED`. Atomic riders: `none`.
Computed `delta_r`: `2` (SymPy Groebner Krull dimension `10` -> `8` with a real non-emptiness witness; real-locus dimension independently confirmed by Mathematica `RegionDimension`; engines agree).
Pre-pass: Route A `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` missing `['R1', 'R2', 'R3', 'R4', 'R5']`; Route B `ROUTE_B_CLOSED_CHECKED_NEGATIVE`; freedom `FREEDOM_UNCONSTRAINED{C_hu,rho_br}`.
Layer 1: E status `SAT`; E plus locks `SAT`.
Layer 2: L_A `WITNESSED`; L_B `WITNESSED`.

## Computation-Cited Riders

- Coupled scalar poles: `det M|cone = -C_hu**2*k**4`; `Delta_pole_minus = (B_eff*M_h - K_h*rho_br - sqrt(B_eff**2*M_h**2 - 2*B_eff*K_h*M_h*rho_br + 4*C_hu**2*M_h*rho_br + K_h**2*rho_br**2))/(2*M_h*rho_br)`; `Delta_pole_plus = (B_eff*M_h - K_h*rho_br + sqrt(B_eff**2*M_h**2 - 2*B_eff*K_h*M_h*rho_br + 4*C_hu**2*M_h*rho_br + K_h**2*rho_br**2))/(2*M_h*rho_br)`. Status: `OFF_CONE_under_AB_proportional_to_C_hu_squared_OPEN_110`.
- Field content: `FIELD_SCALAR_VECTOR_DEPARTURE` from `reports/pathA_39_stage4_field_classification_results.yaml:1-30`.
- Scope caveat: `delta_r` is the lock-constraint slice only; full one-medium drift is NG5.
- Scope caveat: NO_GO non-firing is conditional on `{C_hu, rho_br}` freedom; NG5 must certify and this gate must be rerun if revoked.

## Source-Fact Pre-Pass

| item | status | citation |
|---|---:|---|
| `R1` | `absent` | `['reports/pathA_35_G0_freeze.md:34-38']` |
| `R2` | `absent` | `['reports/pathA_39_scalar_admixture_screen.md:45-58']` |
| `R3` | `absent` | `['reports/pathA_35_G0_freeze.md:83-88']` |
| `R4` | `absent` | `['reports/pathA_35_G0_freeze.md:83-88']` |
| `R5` | `not_applicable` | `['reports/pathA_38_results.yaml:80; reports/pathA_38_throat_body_electric_localization.md:26-30']` |

Route-B checked-negative: h and u_T are distinct Stage-4 DOF; the thin-plate shared-tensor relation is rejected as over-import.

## Algebra

- `L_A = -5*K*rho**4*rho_br + m*mu_R`.
- `L_B = c_E**2*rho_br - mu_R`.
- Strict stability is encoded as `B_eff*K_h - C_hu^2 - sigma = 0` with `sigma>0`.
- Dimension method: SymPy used Groebner initial-monomial Hilbert/Krull dimension plus exact positive witnesses; Mathematica used `Resolve`/`FindInstance` and CAD-backed `RegionDimension`.
- Named dropped assumptions: no generic-Jacobian regular-sequence shortcut; no complex-only dimension/radical without real-locus guard.

## Controls

| control | verdict | Route A | freedom | L_A | L_B | delta_r |
|---|---:|---:|---:|---:|---:|---:|
| `well_posed` | `HALT_ROUTE_A_WELL_POSED` | `ROUTE_A_WELL_POSED` | `FREEDOM_UNCONSTRAINED` | `not_run` | `not_run` | `not_run` |
| `absent` | `CONE_LOCK_CALIBRATED` | `ROUTE_A_ABSENT` | `FREEDOM_UNCONSTRAINED` | `WITNESSED` | `WITNESSED` | `2` |
| `partial_inventory` | `CONE_LOCK_CALIBRATED` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_UNCONSTRAINED` | `WITNESSED` | `WITNESSED` | `2` |
| `forced_lock` | `CONE_LOCK_DERIVED` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_UNCONSTRAINED` | `ENTAILED` | `ENTAILED` | `0` |
| `a_only_partial` | `CONE_LOCK_PARTIAL(derived=A, calibrated=B)` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_UNCONSTRAINED` | `ENTAILED` | `WITNESSED` | `1` |
| `b_only_partial` | `CONE_LOCK_PARTIAL(derived=B, calibrated=A)` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_UNCONSTRAINED` | `WITNESSED` | `ENTAILED` | `1` |
| `over_constrained` | `NO_GO(sector-ledger)` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_UNCONSTRAINED` | `not_run` | `not_run` | `not_run` |
| `freedom_tie` | `NO_GO(cone-lock)` | `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` | `FREEDOM_TIED` | `not_run` | `not_run` | `not_run` |

NO_GO control diagnostics:
- `over_constrained`: `Adding the two equations gives -(sigma+eta_over)=0, impossible with sigma,eta_over>0.` Core: `['B_eff*K_h - C_hu^2 - sigma = 0', 'C_hu^2 - B_eff*K_h - eta_over = 0', 'sigma>0', 'eta_over>0']`.
- `freedom_tie`: `The tie plus L_B and K_h=M_h*c_E^2 gives C_hu^2=2 B_eff K_h; stability then gives -(B_eff*K_h+sigma)=0.` Core includes tie relation `['B_eff*K_h - C_hu^2 - sigma = 0', 'C_hu^2 - q_h_sq = 0', 'q_h_sq*rho_br - 2*B_eff*M_h*mu_R = 0', 'K_h - M_h*c_E^2 = 0', 'c_E^2*rho_br - mu_R = 0', 'B_eff,K_h,rho_br,sigma>0']`.

## Non-blocking review notes

- `over_constrained`/`freedom_tie`: SymPy's UNSAT controls are correct canned certificates keyed on relation kind; Mathematica independently recomputes them as UNSAT with `Resolve`.
- Real-locus dimension: the genuine guarantee comes from Mathematica `RegionDimension`; SymPy supplies Groebner Krull dimension plus a real non-emptiness witness, and the engines agree.
- Mathematica `Delta_pole` strings in the comparison payload are hardcoded literals after validation against its own `FullSimplify` derivation and SymPy's factored form; this is needed because Mathematica symbols are `Unique[]`-mangled.
- `comparison_digest` is informational only, not the agreement test; the actual cross-engine check is structural `comparison_payload` equality over the 150-leaf comparison, which holds.

## Dual Engine

`ENGINE_AGREE` over `150` compared scalar quantities via `STRUCTURAL_PAYLOAD_EQUALITY`.
Informational digests (not the agreement test; engine JSON serialization differs): SymPy `comparison_digest` `e06d8a573653d4163720b18870abab15b4e76e0819f5bc4fa5e2c95ddd6beec0`; Mathematica `comparison_digest` `2b681dc490a7330d01bd882cd58382b9d63bb8226f22dc3fe99858053a216ad9`.
SymPy payload: `/var/projects/toy_physics/software/stage1_solver/_scratch/pathA_40_cone_lock_sympy.json`.
Mathematica payload: `/var/projects/toy_physics/software/stage1_solver/_scratch/pathA_40_cone_lock_mathematica.json`.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_40_cone_lock_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_40_cone_lock.wl
timeout 600 python3 software/stage1_solver/tools/pathA_40_cone_lock_sympy.py --compare
```
