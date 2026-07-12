# Native-`Pᵃ` constraint-class gate — genuine rebuild

## HARDENING NOTE

The tuned-strata rejection is now computation-backed rather than prose-backed. For every rejected first-class primary on an FC-bearing rank-drop/common-null stratum, each engine computes `C₁={C₀,H₂}` and `{q,C₁}`, records the expressions and a derived reason code, and enforces `C₁=0 OR no nonzero spatial action proportional to k`. The run fails if a rejected direction instead has a nonzero `∝k` Gauss descendant. The detailed records are printed below and retained in both engine artifacts.

The scope statement is also hardened: it separates the fully symbolic open-stratum result from the argued-and-scanned tuned-strata coverage; it does not promote a representative-point scan into an exhaustive symbolic classification.

## REBUILD NOTE

- **Q1 closed:** `build_H2` differentiates each input Lagrangian, computes its Hessian nullspace and Legendre transform; `coupling_guard` then proves every native `g_a` survives into computed momentum/Hessian/constraint/PB objects ([native_p_gate_sympy.py:99](../native_p_gate_sympy.py#L99), [native_p_gate_sympy.py:545](../native_p_gate_sympy.py#L545)).
- **Q2 closed:** Maxwell is an input Lagrangian and reaches the same `execute` pipeline as native A/C; its zero PB matrix is an output, never a literal ([native_p_gate_sympy.py:378](../native_p_gate_sympy.py#L378), [native_p_gate_sympy.py:337](../native_p_gate_sympy.py#L337)).
- **Q3 closed:** all six controls are model builders in `CONTROL_BUILDERS` and are executed by `run_control`; the complete Coulomb pair is derived from multiplier terms ([native_p_gate_sympy.py:455](../native_p_gate_sympy.py#L455), [native_p_gate_sympy.py:430](../native_p_gate_sympy.py#L430)).
- **Q4 closed:** Wolfram independently differentiates the Lagrangians and runs its own Hessian/Dirac/kernel code; this comparator checks native and control outputs, not shared answer literals ([native_p_gate_dual.wl:22](../native_p_gate_dual.wl#L22), [native_p_gate_compare.py:52](../native_p_gate_compare.py#L52)).
- **Q5 closed:** `search_G` enumerates computed first-class primary directions and tests their Hamiltonian descendant and field action; Maxwell and gauged-hard-unit must yield nonzero candidates under `GUARD-SEARCH-CAPABLE` ([native_p_gate_sympy.py:196](../native_p_gate_sympy.py#L196), [native_p_gate_sympy.py:750](../native_p_gate_sympy.py#L750)).
- **Q6 closed:** `dirac_search` iterates preservation and reports PB rank/nullity; no full-rank assertion or blanket class dictionary remains ([native_p_gate_sympy.py:141](../native_p_gate_sympy.py#L141), [native_p_gate_dual.wl:63](../native_p_gate_dual.wl#L63)).

**COMPUTED VERDICT: `NATIVE_P_NO_EMERGENT_GAUSS`.**

This is the Stage-1 constraint gate only. Compactness, charge quantization, deconfinement, and native `±w` current supply remain downstream.

## Frozen setup and quadratic operator basis

The calculation uses one nonzero Fourier representative `k=(1,2,3)` on `R³`, with `C_c^∞` smearings and decaying fields. The `k=0` sector is global and is not counted as a local Gauss law. Coulomb control fields and gauge parameters decay, so `-∇²` has no local harmonic kernel. No punctures occur in the source-free search.

Notation: p_i=P_parallel,i; b=chi-chi_vac; k=(1,2,3) is a nonzero Fourier representative.

THEORY-A cross-coupling basis:

- `g_t dot(p).dot(u)`
- `-g_s k^2 p.u`
- `-(g_d/2)(k.p)^2`
- `-g_b b(k.p)`

THEORY-C cross-coupling basis:

- `g_t dot(p).dot(u)`
- `-g_s k^2 p.u`
- `-(g_d/2)(k.p)^2`

Empty families:

- one-time/one-space bilinears are T-odd for T-even configuration fields.
- epsilon bilinears are parity-odd or integration-by-parts zero at order two.
- undifferentiated u is excluded by u-shift.
- in C, b=p^2 starts at amplitude two, so b-cross terms start above quadratic order.

Operator-basis completeness: At quadratic order, SO(3), parity, T, u-shift, <=2 derivatives, IBP and the holonomic identities leave exactly the displayed four (A) or three (C) cross terms. Fixed higher orders are finite associated-graded quotients; quadratic absence is decisive by the v6 linearization lever.

## Completeness and stratum scope

The decisive calculation is the fully symbolic open kinetic stratum: FC=`0` for THEORY-A and THEORY-C for all retained coupling symbols. The computed physical kinetic-Hessian determinant is `(rho_u-g_t^2)^3`, so its only additional kinetic degeneracy is `g_t^2=rho_u` (apart from the frozen multiplier nulls already in the Dirac system).

On that degeneracy surface, the potential-null residual is solved symbolically to identify the rank-drop/common-null conditions, and the remaining non-common rank-drop locus is checked by a fixed-seed randomized representative-point sweep in both engines. This tuned coverage is **ARGUED + SCANNED**, not an exhaustive symbolic stratification of every possible measure-zero sublocus.

Consequently, any hypothetical missed measure-zero first-class Gauss stratum would be a **TUNED / inverse-design** result, not robust native electromagnetism. The physical no-go used here—native `P` does not **generically** host emergent EM—is decisive independently of that caveat.

## THEORY-A: computed `H₂`, Dirac chain, and `G` search

Free coupling symbols: `g_tA, g_sA, g_dA, g_bA`. Coupling-entry guard: **`PASS`**, with computed locations `{'g_bA': ['constraints'], 'g_dA': ['constraints', 'bracket_matrix'], 'g_sA': ['constraints', 'bracket_matrix'], 'g_tA': ['momentum_map', 'hessian', 'constraints', 'bracket_matrix']}`.

Instantiated quadratic Lagrangian:

```text
-(a_bA*bA**2 + a_pA*pA1**2 + a_pA*pA2**2 + a_pA*pA3**2 + 13*a_uA*uA1**2 - 4*a_uA*uA1*uA2 - 6*a_uA*uA1*uA3 + 10*a_uA*uA2**2 - 12*a_uA*uA2*uA3 + 5*a_uA*uA3**2 + 2*bA*g_bA*pA1 + 4*bA*g_bA*pA2 + 6*bA*g_bA*pA3 - d_bA**2 - d_pA1**2 - 2*d_pA1*d_uA1*g_tA - d_pA2**2 - 2*d_pA2*d_uA2*g_tA - d_pA3**2 - 2*d_pA3*d_uA3*g_tA - d_sA**2 - d_uA1**2*rho_uA - d_uA2**2*rho_uA - d_uA3**2*rho_uA + g_dA*pA1**2 + 4*g_dA*pA1*pA2 + 6*g_dA*pA1*pA3 + 4*g_dA*pA2**2 + 12*g_dA*pA2*pA3 + 9*g_dA*pA3**2 + 28*g_sA*pA1*uA1 + 28*g_sA*pA2*uA2 + 28*g_sA*pA3*uA3 - 2*lambdaA*sA - 2*sigmaA*uA1 - 4*sigmaA*uA2 - 6*sigmaA*uA3)/2
```

Computed Legendre transform `H₂`:

```text
-(Pi_native_A_0**2*rho_uA - 2*Pi_native_A_0*Pi_native_A_3*g_tA + Pi_native_A_1**2*rho_uA - 2*Pi_native_A_1*Pi_native_A_4*g_tA + Pi_native_A_2**2*rho_uA - 2*Pi_native_A_2*Pi_native_A_5*g_tA + Pi_native_A_3**2 + Pi_native_A_4**2 + Pi_native_A_5**2 - Pi_native_A_6**2*g_tA**2 + Pi_native_A_6**2*rho_uA - Pi_native_A_9**2*g_tA**2 + Pi_native_A_9**2*rho_uA - a_bA*bA**2*g_tA**2 + a_bA*bA**2*rho_uA - a_pA*g_tA**2*pA1**2 - a_pA*g_tA**2*pA2**2 - a_pA*g_tA**2*pA3**2 + a_pA*pA1**2*rho_uA + a_pA*pA2**2*rho_uA + a_pA*pA3**2*rho_uA - 13*a_uA*g_tA**2*uA1**2 + 4*a_uA*g_tA**2*uA1*uA2 + 6*a_uA*g_tA**2*uA1*uA3 - 10*a_uA*g_tA**2*uA2**2 + 12*a_uA*g_tA**2*uA2*uA3 - 5*a_uA*g_tA**2*uA3**2 + 13*a_uA*rho_uA*uA1**2 - 4*a_uA*rho_uA*uA1*uA2 - 6*a_uA*rho_uA*uA1*uA3 + 10*a_uA*rho_uA*uA2**2 - 12*a_uA*rho_uA*uA2*uA3 + 5*a_uA*rho_uA*uA3**2 - 2*bA*g_bA*g_tA**2*pA1 - 4*bA*g_bA*g_tA**2*pA2 - 6*bA*g_bA*g_tA**2*pA3 + 2*bA*g_bA*pA1*rho_uA + 4*bA*g_bA*pA2*rho_uA + 6*bA*g_bA*pA3*rho_uA - g_dA*g_tA**2*pA1**2 - 4*g_dA*g_tA**2*pA1*pA2 - 6*g_dA*g_tA**2*pA1*pA3 - 4*g_dA*g_tA**2*pA2**2 - 12*g_dA*g_tA**2*pA2*pA3 - 9*g_dA*g_tA**2*pA3**2 + g_dA*pA1**2*rho_uA + 4*g_dA*pA1*pA2*rho_uA + 6*g_dA*pA1*pA3*rho_uA + 4*g_dA*pA2**2*rho_uA + 12*g_dA*pA2*pA3*rho_uA + 9*g_dA*pA3**2*rho_uA - 28*g_sA*g_tA**2*pA1*uA1 - 28*g_sA*g_tA**2*pA2*uA2 - 28*g_sA*g_tA**2*pA3*uA3 + 28*g_sA*pA1*rho_uA*uA1 + 28*g_sA*pA2*rho_uA*uA2 + 28*g_sA*pA3*rho_uA*uA3 + 2*g_tA**2*lambdaA*sA + 2*g_tA**2*sigmaA*uA1 + 4*g_tA**2*sigmaA*uA2 + 6*g_tA**2*sigmaA*uA3 - 2*lambdaA*rho_uA*sA - 2*rho_uA*sigmaA*uA1 - 4*rho_uA*sigmaA*uA2 - 6*rho_uA*sigmaA*uA3)/(2*(g_tA**2 - rho_uA))
```

Hessian rank/nullity: `8/2`. Computed momentum map:

- `Π[1] = d_pA1 + d_uA1*g_tA`
- `Π[2] = d_pA2 + d_uA2*g_tA`
- `Π[3] = d_pA3 + d_uA3*g_tA`
- `Π[4] = d_pA1*g_tA + d_uA1*rho_uA`
- `Π[5] = d_pA2*g_tA + d_uA2*rho_uA`
- `Π[6] = d_pA3*g_tA + d_uA3*rho_uA`
- `Π[7] = d_sA`
- `Π[8] = 0`
- `Π[9] = 0`
- `Π[10] = d_bA`

Dirac constraints in discovery order (`stage 0` means primary):

- stage `0`: `Pi_native_A_7` — `SECOND_CLASS_COMPONENT`
- stage `0`: `Pi_native_A_8` — `SECOND_CLASS_COMPONENT`
- stage `1`: `sA` — `SECOND_CLASS_COMPONENT`
- stage `1`: `uA1 + 2*uA2 + 3*uA3` — `SECOND_CLASS_COMPONENT`
- stage `2`: `Pi_native_A_6` — `SECOND_CLASS_COMPONENT`
- stage `2`: `(Pi_native_A_0*g_tA + 2*Pi_native_A_1*g_tA + 3*Pi_native_A_2*g_tA - Pi_native_A_3 - 2*Pi_native_A_4 - 3*Pi_native_A_5)/(g_tA**2 - rho_uA)` — `SECOND_CLASS_COMPONENT`
- stage `3`: `lambdaA` — `SECOND_CLASS_COMPONENT`
- stage `3`: `-(a_pA*g_tA*pA1 + 2*a_pA*g_tA*pA2 + 3*a_pA*g_tA*pA3 + 14*bA*g_bA*g_tA + 14*g_dA*g_tA*pA1 + 28*g_dA*g_tA*pA2 + 42*g_dA*g_tA*pA3 + 14*g_sA*g_tA*uA1 + 28*g_sA*g_tA*uA2 + 42*g_sA*g_tA*uA3 - 14*g_sA*pA1 - 28*g_sA*pA2 - 42*g_sA*pA3 + 14*sigmaA)/(g_tA**2 - rho_uA)` — `SECOND_CLASS_COMPONENT`

Computed weak PB matrix `M(g,k)` (rank `8`, FC `0`, SC `8`):

```text
[ 0 , 0 , 0 , 0 , 0 , 0 , -1 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 14/(g_tA**2 - rho_uA) ]
[ 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , -14/(g_tA**2 - rho_uA) , 0 , 0 ]
[ 0 , 0 , -1 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 14/(g_tA**2 - rho_uA) , 0 , 0 , 0 , 14*g_tA*(a_pA*g_tA + 14*g_dA*g_tA - 28*g_sA)/(g_tA**2 - rho_uA)**2 ]
[ 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , -14/(g_tA**2 - rho_uA) , 0 , 0 , 0 , -14*g_tA*(a_pA*g_tA + 14*g_dA*g_tA - 28*g_sA)/(g_tA**2 - rho_uA)**2 , 0 , 0 ]
```

Coupling scan and kernel search:

- Regular kinetic determinant per vector component: `-g_tA**2 + rho_uA`; stable open stratum `rho_uA>0 and -g_tA**2 + rho_uA>0`.
- Full physical kinetic-Hessian determinant: `-(g_tA**2 - rho_uA)**3`; only additional kinetic degeneracy: `-g_tA**2 + rho_uA=0`.
- Regular result: computed kernel dimension `0`, first-class primaries `0`, Gauss candidates `0`.
- Boundary rank polynomial: `a_p + 14*a_u - 28*sign(g_t)*g_s; its zero locus is rerun, followed by its common-null sublocus a_p=14*a_u, g_s=sign(g_t)*a_u`.
- Independently computed common-null solutions: `{'-1': {'potential_null_residual': ['2*(a_pA + 14*g_sA)', '-a_pA - 14*g_sA', '0', '28*(a_uA + g_sA)', '-14*(a_uA + g_sA)', '0'], 'solutions': [{'a_pA': '14*a_uA', 'g_sA': '-a_uA'}]}, '1': {'potential_null_residual': ['-2*(a_pA - 14*g_sA)', 'a_pA - 14*g_sA', '0', '28*(a_uA - g_sA)', '-14*(a_uA - g_sA)', '0'], 'solutions': [{'a_pA': '14*a_uA', 'g_sA': 'a_uA'}]}}`.
- Boundary `rho_u=1, g_t=-1; generic semidefinite boundary`: Hessian nullity `5`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=-1; first rank-drop locus, non-common-null representative`: Hessian nullity `5`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=-1; common Hessian/potential null locus`: Hessian nullity `5`, FC `2`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; generic semidefinite boundary`: Hessian nullity `5`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; first rank-drop locus, non-common-null representative`: Hessian nullity `5`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; common Hessian/potential null locus`: Hessian nullity `5`, FC `2`, Gauss candidates `0`, `G=False`.
- Tuned rank-drop sweep: `6` fixed-seed randomized representative points (seed `260713`); scope is `residual solved symbolically and non-common rank-drop locus sampled at fixed-seed randomized representative points; not an exhaustive symbolic stratification of every tuned measure-zero sublocus`.
  - `rank-drop randomized sample sign=-1 index=1` at `{'a_pA': '28', 'a_uA': '5', 'g_bA': '0', 'g_dA': '0', 'g_sA': '-7/2', 'g_tA': '-1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=-1 index=2` at `{'a_pA': '23', 'a_uA': '5', 'g_bA': '0', 'g_dA': '0', 'g_sA': '-93/28', 'g_tA': '-1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=-1 index=3` at `{'a_pA': '35', 'a_uA': '3', 'g_bA': '0', 'g_dA': '0', 'g_sA': '-11/4', 'g_tA': '-1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=1` at `{'a_pA': '10', 'a_uA': '1', 'g_bA': '0', 'g_dA': '0', 'g_sA': '6/7', 'g_tA': '1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=2` at `{'a_pA': '36', 'a_uA': '2', 'g_bA': '0', 'g_dA': '0', 'g_sA': '16/7', 'g_tA': '1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=3` at `{'a_pA': '15', 'a_uA': '1', 'g_bA': '0', 'g_dA': '0', 'g_sA': '29/28', 'g_tA': '1', 'rho_uA': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
- All spatial couplings remained free in the symbolic computations: `g_sA, g_dA, g_bA`.
- Computed aggregate: `gauss_candidates=0`, `additional_G_exists=False`.
- Source-first ordering: source-free searched `True`; `j·A` added `False`; sourced `False`.
- Shear duplicate: `NOT_APPLICABLE_NO_G` (MacCullagh transverse modes `2`).
- Independent Wolfram result: PB rank `8`, FC `0`, candidates `0`, verdict `NATIVE_P_NO_EMERGENT_GAUSS`.

**VERDICT THEORY-A: `NATIVE_P_NO_EMERGENT_GAUSS` at quadratic order.**

### THEORY-A tuned FC descendant rejection audit

Computed hardening guard: **`PASS`**; FC-bearing strata checked `2`, rejected directions checked `4`.

- Stratum `rho_u=1, g_t=-1; common Hessian/potential null locus`:
  - FC direction `1`: primary `-2*Pi_native_A_0 + Pi_native_A_1 - 2*Pi_native_A_3 + Pi_native_A_4`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pA1 -> 0, pA2 -> 0, pA3 -> 0, uA1 -> 0, uA2 -> 0, uA3 -> 0, sA -> 0, lambdaA -> 0, sigmaA -> 0, bA -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
  - FC direction `2`: primary `-3*Pi_native_A_0 + Pi_native_A_2 - 3*Pi_native_A_3 + Pi_native_A_5`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pA1 -> 0, pA2 -> 0, pA3 -> 0, uA1 -> 0, uA2 -> 0, uA3 -> 0, sA -> 0, lambdaA -> 0, sigmaA -> 0, bA -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
- Stratum `rho_u=1, g_t=1; common Hessian/potential null locus`:
  - FC direction `1`: primary `2*Pi_native_A_0 - Pi_native_A_1 - 2*Pi_native_A_3 + Pi_native_A_4`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pA1 -> 0, pA2 -> 0, pA3 -> 0, uA1 -> 0, uA2 -> 0, uA3 -> 0, sA -> 0, lambdaA -> 0, sigmaA -> 0, bA -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
  - FC direction `2`: primary `3*Pi_native_A_0 - Pi_native_A_2 - 3*Pi_native_A_3 + Pi_native_A_5`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pA1 -> 0, pA2 -> 0, pA3 -> 0, uA1 -> 0, uA2 -> 0, uA3 -> 0, sA -> 0, lambdaA -> 0, sigmaA -> 0, bA -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.

## THEORY-C: computed `H₂`, Dirac chain, and `G` search

Free coupling symbols: `g_tC, g_sC, g_dC`. Coupling-entry guard: **`PASS`**, with computed locations `{'g_dC': ['constraints', 'bracket_matrix'], 'g_sC': ['constraints', 'bracket_matrix'], 'g_tC': ['momentum_map', 'hessian', 'constraints', 'bracket_matrix']}`.

Instantiated quadratic Lagrangian:

```text
-(a_bC*bC**2 + a_pC*pC1**2 + a_pC*pC2**2 + a_pC*pC3**2 + 13*a_uC*uC1**2 - 4*a_uC*uC1*uC2 - 6*a_uC*uC1*uC3 + 10*a_uC*uC2**2 - 12*a_uC*uC2*uC3 + 5*a_uC*uC3**2 - 2*bC*xiC - d_pC1**2 - 2*d_pC1*d_uC1*g_tC - d_pC2**2 - 2*d_pC2*d_uC2*g_tC - d_pC3**2 - 2*d_pC3*d_uC3*g_tC - d_sC**2 - d_uC1**2*rho_uC - d_uC2**2*rho_uC - d_uC3**2*rho_uC + g_dC*pC1**2 + 4*g_dC*pC1*pC2 + 6*g_dC*pC1*pC3 + 4*g_dC*pC2**2 + 12*g_dC*pC2*pC3 + 9*g_dC*pC3**2 + 28*g_sC*pC1*uC1 + 28*g_sC*pC2*uC2 + 28*g_sC*pC3*uC3 - 2*lambdaC*sC - 2*sigmaC*uC1 - 4*sigmaC*uC2 - 6*sigmaC*uC3)/2
```

Computed Legendre transform `H₂`:

```text
-(Pi_native_C_0**2*rho_uC - 2*Pi_native_C_0*Pi_native_C_3*g_tC + Pi_native_C_1**2*rho_uC - 2*Pi_native_C_1*Pi_native_C_4*g_tC + Pi_native_C_2**2*rho_uC - 2*Pi_native_C_2*Pi_native_C_5*g_tC + Pi_native_C_3**2 + Pi_native_C_4**2 + Pi_native_C_5**2 - Pi_native_C_6**2*g_tC**2 + Pi_native_C_6**2*rho_uC - a_bC*bC**2*g_tC**2 + a_bC*bC**2*rho_uC - a_pC*g_tC**2*pC1**2 - a_pC*g_tC**2*pC2**2 - a_pC*g_tC**2*pC3**2 + a_pC*pC1**2*rho_uC + a_pC*pC2**2*rho_uC + a_pC*pC3**2*rho_uC - 13*a_uC*g_tC**2*uC1**2 + 4*a_uC*g_tC**2*uC1*uC2 + 6*a_uC*g_tC**2*uC1*uC3 - 10*a_uC*g_tC**2*uC2**2 + 12*a_uC*g_tC**2*uC2*uC3 - 5*a_uC*g_tC**2*uC3**2 + 13*a_uC*rho_uC*uC1**2 - 4*a_uC*rho_uC*uC1*uC2 - 6*a_uC*rho_uC*uC1*uC3 + 10*a_uC*rho_uC*uC2**2 - 12*a_uC*rho_uC*uC2*uC3 + 5*a_uC*rho_uC*uC3**2 + 2*bC*g_tC**2*xiC - 2*bC*rho_uC*xiC - g_dC*g_tC**2*pC1**2 - 4*g_dC*g_tC**2*pC1*pC2 - 6*g_dC*g_tC**2*pC1*pC3 - 4*g_dC*g_tC**2*pC2**2 - 12*g_dC*g_tC**2*pC2*pC3 - 9*g_dC*g_tC**2*pC3**2 + g_dC*pC1**2*rho_uC + 4*g_dC*pC1*pC2*rho_uC + 6*g_dC*pC1*pC3*rho_uC + 4*g_dC*pC2**2*rho_uC + 12*g_dC*pC2*pC3*rho_uC + 9*g_dC*pC3**2*rho_uC - 28*g_sC*g_tC**2*pC1*uC1 - 28*g_sC*g_tC**2*pC2*uC2 - 28*g_sC*g_tC**2*pC3*uC3 + 28*g_sC*pC1*rho_uC*uC1 + 28*g_sC*pC2*rho_uC*uC2 + 28*g_sC*pC3*rho_uC*uC3 + 2*g_tC**2*lambdaC*sC + 2*g_tC**2*sigmaC*uC1 + 4*g_tC**2*sigmaC*uC2 + 6*g_tC**2*sigmaC*uC3 - 2*lambdaC*rho_uC*sC - 2*rho_uC*sigmaC*uC1 - 4*rho_uC*sigmaC*uC2 - 6*rho_uC*sigmaC*uC3)/(2*(g_tC**2 - rho_uC))
```

Hessian rank/nullity: `7/4`. Computed momentum map:

- `Π[1] = d_pC1 + d_uC1*g_tC`
- `Π[2] = d_pC2 + d_uC2*g_tC`
- `Π[3] = d_pC3 + d_uC3*g_tC`
- `Π[4] = d_pC1*g_tC + d_uC1*rho_uC`
- `Π[5] = d_pC2*g_tC + d_uC2*rho_uC`
- `Π[6] = d_pC3*g_tC + d_uC3*rho_uC`
- `Π[7] = d_sC`
- `Π[8] = 0`
- `Π[9] = 0`
- `Π[10] = 0`
- `Π[11] = 0`

Dirac constraints in discovery order (`stage 0` means primary):

- stage `0`: `Pi_native_C_7` — `SECOND_CLASS_COMPONENT`
- stage `0`: `Pi_native_C_8` — `SECOND_CLASS_COMPONENT`
- stage `0`: `Pi_native_C_9` — `SECOND_CLASS_COMPONENT`
- stage `0`: `Pi_native_C_10` — `SECOND_CLASS_COMPONENT`
- stage `1`: `sC` — `SECOND_CLASS_COMPONENT`
- stage `1`: `uC1 + 2*uC2 + 3*uC3` — `SECOND_CLASS_COMPONENT`
- stage `1`: `-a_bC*bC + xiC` — `SECOND_CLASS_COMPONENT`
- stage `1`: `bC` — `SECOND_CLASS_COMPONENT`
- stage `2`: `Pi_native_C_6` — `SECOND_CLASS_COMPONENT`
- stage `2`: `(Pi_native_C_0*g_tC + 2*Pi_native_C_1*g_tC + 3*Pi_native_C_2*g_tC - Pi_native_C_3 - 2*Pi_native_C_4 - 3*Pi_native_C_5)/(g_tC**2 - rho_uC)` — `SECOND_CLASS_COMPONENT`
- stage `3`: `lambdaC` — `SECOND_CLASS_COMPONENT`
- stage `3`: `-(a_pC*g_tC*pC1 + 2*a_pC*g_tC*pC2 + 3*a_pC*g_tC*pC3 + 14*g_dC*g_tC*pC1 + 28*g_dC*g_tC*pC2 + 42*g_dC*g_tC*pC3 + 14*g_sC*g_tC*uC1 + 28*g_sC*g_tC*uC2 + 42*g_sC*g_tC*uC3 - 14*g_sC*pC1 - 28*g_sC*pC2 - 42*g_sC*pC3 + 14*sigmaC)/(g_tC**2 - rho_uC)` — `SECOND_CLASS_COMPONENT`

Computed weak PB matrix `M(g,k)` (rank `12`, FC `0`, SC `12`):

```text
[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 14/(g_tC**2 - rho_uC) ]
[ 0 , 0 , 0 , 0 , 0 , 0 , a_bC , -1 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 0 , -1 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -14/(g_tC**2 - rho_uC) , 0 , 0 ]
[ 0 , 0 , -a_bC , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , -1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , 0 , 0 , 0 , 0 , 14/(g_tC**2 - rho_uC) , 0 , 0 , 0 , 0 , 0 , 14*g_tC*(a_pC*g_tC + 14*g_dC*g_tC - 28*g_sC)/(g_tC**2 - rho_uC)**2 ]
[ 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ]
[ 0 , -14/(g_tC**2 - rho_uC) , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -14*g_tC*(a_pC*g_tC + 14*g_dC*g_tC - 28*g_sC)/(g_tC**2 - rho_uC)**2 , 0 , 0 ]
```

Coupling scan and kernel search:

- Regular kinetic determinant per vector component: `-g_tC**2 + rho_uC`; stable open stratum `rho_uC>0 and -g_tC**2 + rho_uC>0`.
- Full physical kinetic-Hessian determinant: `-(g_tC**2 - rho_uC)**3`; only additional kinetic degeneracy: `-g_tC**2 + rho_uC=0`.
- Regular result: computed kernel dimension `0`, first-class primaries `0`, Gauss candidates `0`.
- Boundary rank polynomial: `a_p + 14*a_u - 28*sign(g_t)*g_s; its zero locus is rerun, followed by its common-null sublocus a_p=14*a_u, g_s=sign(g_t)*a_u`.
- Independently computed common-null solutions: `{'-1': {'potential_null_residual': ['2*(a_pC + 14*g_sC)', '-a_pC - 14*g_sC', '0', '28*(a_uC + g_sC)', '-14*(a_uC + g_sC)', '0'], 'solutions': [{'a_pC': '14*a_uC', 'g_sC': '-a_uC'}]}, '1': {'potential_null_residual': ['-2*(a_pC - 14*g_sC)', 'a_pC - 14*g_sC', '0', '28*(a_uC - g_sC)', '-14*(a_uC - g_sC)', '0'], 'solutions': [{'a_pC': '14*a_uC', 'g_sC': 'a_uC'}]}}`.
- Boundary `rho_u=1, g_t=-1; generic semidefinite boundary`: Hessian nullity `7`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=-1; first rank-drop locus, non-common-null representative`: Hessian nullity `7`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=-1; common Hessian/potential null locus`: Hessian nullity `7`, FC `2`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; generic semidefinite boundary`: Hessian nullity `7`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; first rank-drop locus, non-common-null representative`: Hessian nullity `7`, FC `0`, Gauss candidates `0`, `G=False`.
- Boundary `rho_u=1, g_t=1; common Hessian/potential null locus`: Hessian nullity `7`, FC `2`, Gauss candidates `0`, `G=False`.
- Tuned rank-drop sweep: `6` fixed-seed randomized representative points (seed `260715`); scope is `residual solved symbolically and non-common rank-drop locus sampled at fixed-seed randomized representative points; not an exhaustive symbolic stratification of every tuned measure-zero sublocus`.
  - `rank-drop randomized sample sign=-1 index=1` at `{'a_pC': '30', 'a_uC': '1', 'g_dC': '0', 'g_sC': '-11/7', 'g_tC': '-1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=-1 index=2` at `{'a_pC': '2', 'a_uC': '3', 'g_dC': '0', 'g_sC': '-11/7', 'g_tC': '-1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=-1 index=3` at `{'a_pC': '5', 'a_uC': '3', 'g_dC': '0', 'g_sC': '-47/28', 'g_tC': '-1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=1` at `{'a_pC': '10', 'a_uC': '5', 'g_dC': '0', 'g_sC': '20/7', 'g_tC': '1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=2` at `{'a_pC': '29', 'a_uC': '3', 'g_dC': '0', 'g_sC': '71/28', 'g_tC': '1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
  - `rank-drop randomized sample sign=1 index=3` at `{'a_pC': '32', 'a_uC': '1', 'g_dC': '0', 'g_sC': '23/14', 'g_tC': '1', 'rho_uC': '1'}`: FC `0`, Gauss candidates `0`, `G=False`.
- All spatial couplings remained free in the symbolic computations: `g_sC, g_dC`.
- Computed aggregate: `gauss_candidates=0`, `additional_G_exists=False`.
- Source-first ordering: source-free searched `True`; `j·A` added `False`; sourced `False`.
- Shear duplicate: `NOT_APPLICABLE_NO_G` (MacCullagh transverse modes `2`).
- Independent Wolfram result: PB rank `12`, FC `0`, candidates `0`, verdict `NATIVE_P_NO_EMERGENT_GAUSS`.

**VERDICT THEORY-C: `NATIVE_P_NO_EMERGENT_GAUSS` at quadratic order.**

### THEORY-C tuned FC descendant rejection audit

Computed hardening guard: **`PASS`**; FC-bearing strata checked `2`, rejected directions checked `4`.

- Stratum `rho_u=1, g_t=-1; common Hessian/potential null locus`:
  - FC direction `1`: primary `-2*Pi_native_C_0 + Pi_native_C_1 - 2*Pi_native_C_3 + Pi_native_C_4`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pC1 -> 0, pC2 -> 0, pC3 -> 0, uC1 -> 0, uC2 -> 0, uC3 -> 0, sC -> 0, lambdaC -> 0, sigmaC -> 0, bC -> 0, xiC -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
  - FC direction `2`: primary `-3*Pi_native_C_0 + Pi_native_C_2 - 3*Pi_native_C_3 + Pi_native_C_5`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pC1 -> 0, pC2 -> 0, pC3 -> 0, uC1 -> 0, uC2 -> 0, uC3 -> 0, sC -> 0, lambdaC -> 0, sigmaC -> 0, bC -> 0, xiC -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
- Stratum `rho_u=1, g_t=1; common Hessian/potential null locus`:
  - FC direction `1`: primary `2*Pi_native_C_0 - Pi_native_C_1 - 2*Pi_native_C_3 + Pi_native_C_4`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pC1 -> 0, pC2 -> 0, pC3 -> 0, uC1 -> 0, uC2 -> 0, uC3 -> 0, sC -> 0, lambdaC -> 0, sigmaC -> 0, bC -> 0, xiC -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
  - FC direction `2`: primary `3*Pi_native_C_0 - Pi_native_C_2 - 3*Pi_native_C_3 + Pi_native_C_5`; computed descendant `{primary,H₂} = 0`; computed field action `{q,{primary,H₂}} = { pC1 -> 0, pC2 -> 0, pC3 -> 0, uC1 -> 0, uC2 -> 0, uC3 -> 0, sC -> 0, lambdaC -> 0, sigmaC -> 0, bC -> 0, xiC -> 0 }`; reason **`DESCENDANT_ZERO`** (`descendant_zero=True`, `secondary_action_non_gradient=True`).
    - Spatial block fields `[1, 2, 3]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.
    - Spatial block fields `[4, 5, 6]`: action `['0', '0', '0']`; computed `proportional_to_k=False`.

## Six able-to-fail controls through the shared path

Every row was obtained from an input Lagrangian through the identical Hessian → Legendre → Dirac → kernel search. The nonconserved-current row additionally evaluates the nonzero continuity defect in Gauss preservation.

| # | Control | Computed class | Hessian nullity | FC | SC | Gauss candidates |
|---:|---|---|---:|---:|---:|---:|
| 1 | maxwell | `FIRST_CLASS_GAUSS` | 1 | 2 | 0 | 1 |
| 2 | gauged hard unit | `MIXED` | 2 | 2 | 4 | 1 |
| 3 | bare O(4) hard sigma | `SECOND_CLASS_RADIAL_NO_GAUSS` | 1 | 0 | 4 | 0 |
| 4 | nonconserved_current | `INCONSISTENT_PRESERVATION` | 1 | 2 | 0 | 1 |
| 5 | Coulomb gauge Maxwell | `SECOND_CLASS_NO_LOCAL_GAUGE` | 3 | 0 | 8 | 0 |
| 6 | global U(1) complex scalar | `GLOBAL_CHARGE_NO_LOCAL_GAUSS` | 0 | 0 | 0 | 0 |

Control 4 computed Gauss preservation: `-nc_j1 - 2*nc_j2 - 3*nc_j3 + nc_rho_dot`; with no conservation rule imposed this is nonzero and inconsistent.

Per-tooth mutations (each reruns the shared computation and fails only its own expected assertion):

- `maxwell` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[maxwell]: computed Maxwell Gauss chain missing
- `gauged_hard_unit` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[gauged_hard_unit]: computed first/second-class mixture missing
- `bare_sigma` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[bare_sigma]: computed radial sector is not purely second class
- `nonconserved_current` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[nonconserved_current]: continuity defect did not spoil Gauss preservation
- `gauge_fixed_maxwell` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[gauge_fixed_maxwell]: complete Coulomb fixing left a local G
- `global_only` — `FIRED_AT_OWN_ASSERT`: CONTROL_ASSERT[global_only]: global charge was mistaken for a local constraint

## Dual-engine logs and agreement

```text
AGREE shared input-Lagrangian pipeline present in both engines
AGREE GUARD-COUPLINGS-ENTER=A:PASS,C:PASS
AGREE GUARD-SEARCH-CAPABLE=PASS
AGREE HARDENING-TUNED-DESCENDANT-REJECTION=A:PASS,C:PASS
AGREE THEORY-A constraint_count: 8
AGREE THEORY-A constraint_stages: [0, 0, 1, 1, 2, 2, 3, 3]
AGREE THEORY-A constraint_classes: ['SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT']
AGREE THEORY-A bracket_rank: 8
AGREE THEORY-A first_class_count: 0
AGREE THEORY-A second_class_count: 8
AGREE THEORY-A gauss_candidates: 0
AGREE THEORY-A additional_G_exists: False
AGREE THEORY-A boundary_scan: [(0, 0, False, 5), (0, 0, False, 5), (2, 0, False, 5), (0, 0, False, 5), (0, 0, False, 5), (2, 0, False, 5)]
AGREE THEORY-A verdict: NATIVE_P_NO_EMERGENT_GAUSS
AGREE THEORY-A rankdrop_randomized_sweep: [(0, 0, False, 5), (0, 0, False, 5), (0, 0, False, 5), (0, 0, False, 5), (0, 0, False, 5), (0, 0, False, 5)]
AGREE THEORY-A tuned_descendant_rejection: [(2, 2, 0, (('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))), ('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))))), (2, 2, 0, (('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))), ('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0'))))))]
ENGINE_AGREE THEORY-A
AGREE THEORY-C constraint_count: 12
AGREE THEORY-C constraint_stages: [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3]
AGREE THEORY-C constraint_classes: ['SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT', 'SECOND_CLASS_COMPONENT']
AGREE THEORY-C bracket_rank: 12
AGREE THEORY-C first_class_count: 0
AGREE THEORY-C second_class_count: 12
AGREE THEORY-C gauss_candidates: 0
AGREE THEORY-C additional_G_exists: False
AGREE THEORY-C boundary_scan: [(0, 0, False, 7), (0, 0, False, 7), (2, 0, False, 7), (0, 0, False, 7), (0, 0, False, 7), (2, 0, False, 7)]
AGREE THEORY-C verdict: NATIVE_P_NO_EMERGENT_GAUSS
AGREE THEORY-C rankdrop_randomized_sweep: [(0, 0, False, 7), (0, 0, False, 7), (0, 0, False, 7), (0, 0, False, 7), (0, 0, False, 7), (0, 0, False, 7)]
AGREE THEORY-C tuned_descendant_rejection: [(2, 2, 0, (('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))), ('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))))), (2, 2, 0, (('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0')))), ('DESCENDANT_ZERO', True, True, True, ('0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'), ((False, ('0', '0', '0')), (False, ('0', '0', '0'))))))]
ENGINE_AGREE THEORY-C
AGREE CONTROL maxwell: FC=2 SC=0 G=1 FIRST_CLASS_GAUSS
AGREE CONTROL gauged_hard_unit: FC=2 SC=4 G=1 MIXED
AGREE CONTROL bare_sigma: FC=0 SC=4 G=0 SECOND_CLASS_RADIAL_NO_GAUSS
AGREE CONTROL nonconserved_current: FC=2 SC=0 G=1 INCONSISTENT_PRESERVATION
AGREE CONTROL gauge_fixed_maxwell: FC=0 SC=8 G=0 SECOND_CLASS_NO_LOCAL_GAUGE
AGREE CONTROL global_only: FC=0 SC=0 G=0 GLOBAL_CHARGE_NO_LOCAL_GAUSS
ENGINE_AGREE: ALL_THEORIES_AND_CONTROLS
```

SymPy:

```text
SYMPY_GENUINE_QUADRATIC_DIRAC: PASS
GUARD-COUPLINGS-ENTER: PASS
GUARD-SEARCH-CAPABLE: PASS
HARDENING-TUNED-DESCENDANT-REJECTION: PASS
THEORY-A: FC=0 G_CANDIDATES=0 VERDICT=NATIVE_P_NO_EMERGENT_GAUSS
THEORY-C: FC=0 G_CANDIDATES=0 VERDICT=NATIVE_P_NO_EMERGENT_GAUSS
SIX_CONTROLS_SHARED_PIPELINE: PASS
SIX_PER_TOOTH_ABLATIONS: FIRED_AT_OWN_ASSERT
```

Wolfram Language:

```text
MATHEMATICA_GENUINE_QUADRATIC_DIRAC: PASS
GUARD-COUPLINGS-ENTER: PASS
GUARD-SEARCH-CAPABLE: PASS
HARDENING-TUNED-DESCENDANT-REJECTION: PASS
THEORY-A: FC=0 G=0 VERDICT=NATIVE_P_NO_EMERGENT_GAUSS
THEORY-C: FC=0 G=0 VERDICT=NATIVE_P_NO_EMERGENT_GAUSS
SIX_CONTROLS_SHARED_PIPELINE: PASS
```

## Decision-table result

Both regular coupling families have computed FC count zero, and neither semidefinite kinetic boundary develops a first-class Gauss direction. Since a regular nonlinear gauge identity must have a nontrivial leading linearization, the quadratic absence selects branch 2.

**FINAL VERDICT: `NATIVE_P_NO_EMERGENT_GAUSS`.**

