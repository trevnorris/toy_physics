DENSITY_PORT_HOSTED

# PathA-43 Density Quadrupole Port

The remediated gate computes a hosted reduced density-port structure, not a
literal throat magnitude. Both engines produce the same two-port numerator
from `(q2, Phi2)`:

`N0_den = I25^2 Xi_Q^2 c_s^4 rho (eta_phi varpi_q2 + eta_q lambda_c)^2 / (a^7 (lambda_c^2 - varpi_Phi2 varpi_q2)^2)`.

SymPy uses the full inverse of the real `2x2` static operator. Mathematica uses
the independent Green-function/DtN route, eliminating `q2` before substituting
into the `Phi2` equation. The agreement is symbolic; neither engine consumes
the other's numerator or booleans.

Hardening folded here: `source_map_complete` now requires every symbol in
`N0_den.free_symbols` to be present in the source map with a non-empty
provenance tag set. A missing host symbol or a host symbol with an empty tag set
fails origin and verdicts `FAIL_NOT_DENSITY_DERIVED`.

## Computed Source Graph

The taint set is computed from `N0_den.free_symbols`, using the code's
symbol-to-source map:

| symbol | tags |
| --- | --- |
| `I25` | `continuity_interface`, `pathA_32_wall` |
| `Xi_Q` | `continuity_interface` |
| `a` | `pathA_29_bulk`, `pathA_32_wall` |
| `c_s` | `pathA_29_bulk` |
| `eta_phi` | `continuity_interface` |
| `eta_q` | `continuity_interface` |
| `lambda_c` | `continuity_interface`, `pathA_29_bulk`, `pathA_32_wall` |
| `rho` | `pathA_29_bulk` |
| `varpi_Phi2` | `pathA_29_bulk` |
| `varpi_q2` | `pathA_32_wall` |

Computed taint set:
`{continuity_interface, pathA_29_bulk, pathA_32_wall}`.

Every production host symbol listed above has a non-empty tag set, so the
baseline remains `DENSITY_PORT_HOSTED`.

`vector_host_symbols` is empty. `VECTOR_SYMBOLS` is referenced directly by the
host guard and expression-level ablation.

## Continuity Dependency

`I25` is no longer a trusted standalone placeholder. The engines build the l=2
moment symbol through a continuity-moment constructor whose lineage must contain
the PathA-29 operator, the `M0` moment, the `D1_i` moment, and the `Q2_m =
Integral(Y2_m_star*S_leak d3x)` construction. The validator checks structural
tokens, not a shared `valid` flag.

The literal value of the `Y2*` moment, `Xi_Q`, `eta_q`, `eta_phi`, and
`lambda_c` remains `SIM_DEFERRED`. What is hosted here is the reduced coupling
structure, vector-freedom, continuity ancestry, dimensions, a-scaling, DtN sign,
and closure overlay.

## Checks

| check | result |
| --- | --- |
| origin_derivation_ok | PASS |
| nonzero | PASS |
| dimension `[N0_den]=L^-1 M` | PASS |
| a-scaling provenance, `P0_phys` a-power `-5` | PASS |
| radiative sign, outgoing `+i z^5/27` | PASS |
| vector-independence expression ablation | PASS, delta `0` |

Port picture: `ii two-port(q2,Phi2)`.

## Closure

`P0_physical = (c_s/a)^2 N0_den/D0`.

`Kbar4 - 4 Kbar2^2/Kbar0 = 0`, and the calibrated gamma residual to
`2G/(5c^5)` is `0`. The constants `G`, `2/5`, and `54/5` remain
`CALIBRATED`.

## Controls

| control | recomputed verdict |
| --- | --- |
| vector_hosted | FAIL_NOT_DENSITY_DERIVED |
| relabel_rig | FAIL_NOT_DENSITY_DERIVED |
| fake_continuity | FAIL_NOT_DENSITY_DERIVED |
| ablation_isolation | FAIL_NOT_DENSITY_DERIVED |
| attack2_continuity_corruption | FAIL_NOT_DENSITY_DERIVED |
| attack5_vector_injection | FAIL_NOT_DENSITY_DERIVED |
| zero_coupling | FAIL_PORT_VANISHES |
| dimensional | FAIL_PORT_MALFORMED(dimensional) |
| sign | FAIL_PORT_MALFORMED(sign) |
| scaling | FAIL_PORT_MALFORMED(scaling) |
| deferred_scalar | PORT_INCONCLUSIVE_SIM_DEFERRED |
| deferred_scalar_proven_converse | DENSITY_PORT_HOSTED |
| free_carrier_dimension_corruption | DENSITY_PORT_HOSTED |
| provenance_less_rider | FAIL_NOT_DENSITY_DERIVED |

`free_carrier_dimension_corruption` now multiplies in a dimensionless,
a-neutral `free_carrier` with explicit
`pathA_34_dimensionless_free_carrier` provenance, preserving the PathA-34
sourced-vs-free dimensional lesson. `provenance_less_rider` multiplies in the
same carrier with an empty tag set; both engines set `source_map_complete` to
`FAIL` and recompute `FAIL_NOT_DENSITY_DERIVED`.

External source-mutated copies were also checked under `timeout 600`:
corrupting the `Q2_m` construction and injecting `Omega_U/Omega_W` into the
production density coupling both recomputed `FAIL_NOT_DENSITY_DERIVED`.

Engine stdout evidence:

- `software/stage1_solver/reports/pathA_43_density_quadrupole_port_sympy.txt`
- `software/stage1_solver/reports/pathA_43_density_quadrupole_port_mathematica.txt`
