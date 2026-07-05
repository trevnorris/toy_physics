SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED

# pathA_42 Charge-Coupled Scalar Gate

Computed overall verdict: `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`.
The headline is computed from the five channel states and guard predicates by the directive first-match rules.

## Five-Channel Table

| channel | state | subtags | adjudication |
|---|---:|---|---|
| `h_EP` | `EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR` | `none` | Mass-channel safe on the decoupled floor only; full decoupled-floor EP safety also needs universality=EARNED_SAFE. |
| `radiation` | `SIM_GATED` | `CHERENKOV_DEFERRED` | Current ledger branch is exponent-3; magnitude remains SIM_GATED on the free h-sector kinetic normalization. |
| `universality` | `SIM_GATED_REQUIRED_UNIVERSALITY` | `none` | q_h/Q_E universality is required but not earned from the current b/ell ledger. |
| `u_L_EP` | `SIM_GATED` | `MIXED_SCALAR_EP_RISK` | u_L couples to charge and mass; C_hu mixing transfers mass coupling into charge-sourced eigenmodes. |
| `preferred_frame` | `SIM_GATED` | `none` | c_E=c_gamma is calibrated, not earned; large-c_E hiding would carry preferred-frame cost. |

## EP And Mixing

`h_EP` is `EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR` only for the Gate-L mass-coupled channel on the decoupled floor.
Full decoupled-floor EP safety is not reported because it additionally requires `EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR AND universality=EARNED_SAFE`.
C_hu*q_h*q_M term reintroduces a mass coupling when mixing is nonzero.

## Radiation

Setup: nonrelativistic point-charge dipole; scalar and EM far-zone fields at c_E/c_gamma; scalar flux vs Poynting flux.
Derived flux powers: scalar `d**2*omega**4*q_h**2/(M_h*c_E**3)`, EM `Q_E**2*d**2*omega**4/c_gamma**3`.
The current ledger selects the `exponent-3 bare-fixed/current branch`: exponent `3`.
The exponent-`1` branch is conditional on a future pinned-`K_h` fact; the default pinned-`K_h` test is `False`.
At `c_E=c_gamma` without pinned `K_h`, radiation remains `SIM_GATED`; with pinned `K_h`, control K reaches `FALSIFIABLE_TENSION`.
Magnitude status: `SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED`. Subtag: `CHERENKOV_DEFERRED`.

## Static Exchange

`S = [[B_eff, C_hu], [C_hu, K_h]]`, determinant `B_eff*K_h - C_hu**2`.
`A_qq = (B_eff*q_h**2 - 2*C_hu*q_L*q_h + K_h*q_L**2)/(B_eff*K_h - C_hu**2)`
`A_qm = (B_eff*m_h*q_h - C_hu*m_h*q_L - C_hu*q_M*q_h + K_h*q_L*q_M)/(B_eff*K_h - C_hu**2)`
`A_mm = (B_eff*m_h**2 - 2*C_hu*m_h*q_M + K_h*q_M**2)/(B_eff*K_h - C_hu**2)`
With `q_L=0`, nonzero mixing gives `-C_hu*q_M*q_h/(B_eff*K_h - C_hu**2)`, so `MIXED_SCALAR_EP_RISK` is computed.

## HARD_WALL And Guard A

`HARD_WALL`: `h_time_kinetic_parent_action` is absent; `M_h`, `c_E`, `K_h`, radiation magnitude, and EP magnitude are not pinned.
Guard A refuses all laundering fixtures: `M_h_from_N0`=LAUNDERING_REFUSED, `M_h_from_K_parallel`=LAUNDERING_REFUSED, `c_E_from_c_gamma_Green`=LAUNDERING_REFUSED, `K_h_from_N0_cgamma2`=LAUNDERING_REFUSED.
Direct assembled-result injection: `LAUNDERING_REFUSED guarded_numeric_paths=P_h/P_EM,M_h`.
Forged local flag: `LAUNDERING_REFUSED`; wired check `LAUNDERING_REFUSED guarded_numeric_paths=P_h_over_P_EM`.
Guard A bypass regressions: `tuple_wrapped_guarded_number` -> `LAUNDERING_REFUSED guarded_numeric_paths=assembled_results_fixture.power_ratio.0`, `string_embedded_guarded_number` -> `LAUNDERING_REFUSED guarded_numeric_paths=assembled_results_fixture.power_ratio`.
Guard A is a denylist scoped to M_h, c_E, K_h, P_h/P_EM, EP magnitude, and residue floor; unrelated new numeric fields are out of scope.
The guarded residue-floor attempt is also refused by the serialization guard; no numeric scalar magnitude is emitted.

## Preferred Frame, u_L, Gate-L

Preferred-frame state: `SIM_GATED`. Large c_E/c_gamma radiation hiding would trigger PREFERRED_FRAME_TENSION.
`u_L`/mixing state: `SIM_GATED`, magnitude `SIM_GATED via a_L AND C_hu`.
Gate-L connection: same embedding-direction family, distinct ledger object; ungapped bending family triggers Gate-L fifth-force failure.

## Controls

| control | status | transition | key result |
|---|---:|---|---|
| `A` | `FIRED` | `serialization guard wired; production emits no guarded number` | fixtures, bypass regressions, and direct output injection -> LAUNDERING_REFUSED |
| `B` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> FALSIFIABLE_TENSION(channel=h_EP)` | EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR -> FIFTH_FORCE_TRIGGERED |
| `C` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | q_L=0 with C_hu*q_h*q_M nonzero -> MIXED_SCALAR_EP_RISK |
| `D` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | radiation exponent changes under corrupted Green/flux speed |
| `E` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | SIM_GATED -> PREFERRED_FRAME_TENSION; NATURALLY_HIDDEN unreachable |
| `F` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> FALSIFIABLE_TENSION(channel=universality)` | SIM_GATED_REQUIRED_UNIVERSALITY -> FALSIFIABLE_TENSION |
| `G` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> NO_GO_CONSISTENCY` | strict stability C_hu**2 < B_eff*K_h -> STABILITY_VIOLATED |
| `H` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> NO_GO_CONSISTENCY` | static-Coulomb match perturbed -> h_EP=NO_GO |
| `I` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | q_M sign flip -> A_qm and signed A_mm/q_M projection flip |
| `J` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | bare-fixed exponent 3 -> pinned-K_h exponent 1 under wrong normalization fixture |
| `K_without_pinned_Kh` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | c_E=c_gamma without pinned K_h -> radiation=SIM_GATED |
| `K` | `FIRED` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED -> FALSIFIABLE_TENSION(channel=radiation)` | pinned K_h plus c_E=c_gamma -> radiation=FALSIFIABLE_TENSION |

## Deletion Sensitivity

| deleted/stubbed channel | recomputed verdict | changed |
|---|---|---:|
| `h_EP` | `NO_GO_CONSISTENCY` | `True` |
| `radiation` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | `False` |
| `universality` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | `False` |
| `u_L_EP` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | `False` |
| `preferred_frame` | `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED` | `False` |
| `radiation+universality+u_L_EP+preferred_frame` | `NATURALLY_HIDDEN` | `True` |

## Dual Engine

`ENGINE_AGREE` over `22` compared leaves.

Run commands:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_42_charge_coupled_scalar_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_42_charge_coupled_scalar.wl
timeout 600 python3 software/stage1_solver/tools/pathA_42_charge_coupled_scalar_sympy.py --compare
```
