# S9 formal coverage

Updated 2026-09-11. The original Lean pilot is complete within its stated smooth,
constant-coefficient, whole-spacetime plane-wave setting. This does not close
every claim in the [S9 ledger record](../../steps/S9_light_requires_shear.md).

| Claim or obligation | Status |
|---|---|
| Supplied D=3 action to integrated first variation and local PDE | Proved for smooth backgrounds and smooth compact test variations; finite total action is not required. |
| Plane-wave reduction and complete census within that ansatz | Proved: two transverse propagating directions and one longitudinal static direction under positive coefficients and nonzero wavevector. |
| Agreement with the later S10 baseline | Exact D=3 identities are proved for the action, operator, PDE, integrated stationarity and mode spaces. |
| Arbitrary-D extension, dimensions and control actions | Covered by the later S10 libraries; see their coverage map. These extensions do not certify the original S9 CAS emissions. |
| P2: scalar GNLS/Madelung has no transverse propagating mode in the intended linearized regime | Not formalized. This remains the principal separate mathematical gap in the S9 ledger argument. |
| Original S9 CAS scripts, parser/comparator, exports and prose | Not formally connected to the Lean objects. Historical computational evidence remains distinct. |
| General PDE solution completeness beyond the plane-wave ansatz | Not proved; requires an additional analytic setting and result. |
| Finite boundaries, interfaces, curved domains and weak solutions | Outside this pilot's setting. |
| Microscopic origin of the action, bulk shear-freeness and confinement | Supplied physical premises or separate work; not consequences of this action proof. |

The static longitudinal result does not remove that degree of freedom: its
restoring stiffness vanishes in the supplied action. A theorem about that action
does not independently justify choosing it.

The [proof report](RESULT.md) and [verification record](VERIFICATION.txt) give
the precise S9 statements and 21 selected audits. The later
[S10 coverage map](../s10/COVERAGE.md) records the additional formal coverage
and the remaining integration work. The [combined checkpoint](../CHECKPOINT.md)
records the state being committed.
