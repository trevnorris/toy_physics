# S9 formal coverage

Updated 2026-09-15 (closure recorded 2026-09-16 UTC). The original Lean pilot is complete within its stated smooth,
constant-coefficient, whole-spacetime plane-wave setting. This does not close
every claim in the [S9 ledger record](../../steps/S9_light_requires_shear.md).

## Bounded completion contract

This is the user-authorized S9 closure task under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Status: **complete within C1–C4**;
both independent fidelity reviews are CLEAR. The reviewed revision, evidence and
dispositions are recorded in [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md).
It reuses the original action proof and S10 specialization rather than rebuilding
their CAS bridge. The finite obligations are:

| Item | Claim, evidence and completion requirement |
|---|---|
| C1 — Existing action and coverage | Reuse the supplied D=3 curl action, its integrated variational theorem, and the exhaustive plane-wave classification for `rho > 0`, `mu > 0`, `k ≠ 0`: at `omega = 0` the kernel is `span{k}` of dimension 1; at nonzero real frequency on `omega² = (mu/rho)|k|²` it is `k`-orthogonal of dimension 2; elsewhere it is zero. Reuse the exact S10-to-S9 action, operator, PDE and stationarity identities. Record the original S9 SymPy/Wolfram constructor, parameter/coordinate map, sign and normalization; check that compact connection without importing every CAS emission into Lean. |
| C2 — Scalar phase velocity | From the supplied Madelung law `v = (hbar/m) grad(theta)`, differentiate a smooth, single-valued scalar phase `theta0 + epsilon (A cos(k·x − omega t) + B sin(k·x − omega t))` about a constant phase, in the uniform, nonzero-density, vortex-free regime where that law is applicable. Prove the full pair of linear velocity quadratures lies in `span{k}`; for `k ≠ 0`, its intersection with the transverse space is zero. Include nonzero longitudinal examples and the `k = 0` limit. This is a necessary kinematic constraint for such GNLS perturbations; no GNLS dispersion or universal statement about emergent excitations is asserted. |
| C3 — Mutation controls | Record passing counterparts and isolated mathematical failures for load-bearing action/variation signs, mode counts/domain assumptions, and the new phase-gradient identification and exclusion. Reuse relevant S10 controls only through the documented exact specialization. Keep source provenance and diagnostics; syntax/import/resource failures do not count. |
| C4 — Verification and fidelity | Run appropriate sequential Lean checks and axiom audits without admissions/custom physics axioms. Obtain two independent non-author reviews of the fixed statements, compact source connection, cases and controls; resolve findings and record the reviewed source hashes. |

The action coefficients are positive real constants; coordinates are ordered
`(t,x,y,z)` and the displacement is a real three-vector. Time and space have
their ledger units; `mu/rho` has squared-velocity units. The phase is
dimensionless, `hbar/m` has units length²/time, and its spatial gradient has
units inverse length. Unit assignments and the Madelung law are supplied
physical identifications, not additional kernel-derived physics.

For C1 the three frequency cases are exhaustive and disjoint on the stated
domain; the dimensions concern entire amplitude kernels, not isolated vectors
or determinant multiplicities. The two signs of the nonzero frequency share
the same squared-frequency cone. For C2, all real quadrature coefficients and
frequencies are allowed: imposing a scalar evolution equation can only restrict
that family. Positive uniform density supplies the interpretation of phase
velocity; the algebraic proof does not manufacture a density hypothesis it
does not use. Physical `hbar,m > 0` ensure the longitudinal family is nonempty.

Exclusions: deriving the supplied actions or Madelung law, proving the full
Bogoliubov/GNLS dispersion relation or number of dynamical branches, arbitrary
PDE/Fourier completeness, vortices, interfaces, confinement, microscopic photon
claims, S11 work, and certification of the original export/parser/comparator
pipeline. C1's source connection is a checked translation boundary, not a
kernel proof that Python or Wolfram executed correctly. Stop when C1–C4 close.

| Claim or obligation | Status |
|---|---|
| Supplied D=3 action to integrated first variation and local PDE | Proved for smooth backgrounds and smooth compact test variations; finite total action is not required. |
| Plane-wave reduction and complete census within that ansatz | Proved: two transverse propagating directions and one longitudinal static direction under positive coefficients and nonzero wavevector. |
| Agreement with the later S10 baseline | Exact D=3 identities are proved for the action, operator, PDE, integrated stationarity and mode spaces. |
| Arbitrary-D extension, dimensions and control actions | Covered by the later S10 libraries; see their coverage map. These extensions do not certify the original S9 CAS emissions. |
| Narrow scalar-phase part of P2 | `Madelung.lean` proves that both velocity quadratures obtained from the actual supplied phase-gradient law are longitudinal, with no nonzero transverse polarization. The S9 rebuild, 32 axiom audits, nine rejected mutations and five positive controls pass; both independent fidelity reviews are CLEAR. The broader spin-1/GNLS branch claim remains outside this theorem. |
| Original S9 CAS action/operator | Compact source connection recorded in [FIDELITY.md](FIDELITY.md): exact symbolic tests of the native SymPy MAIN action and both modal routes pass; Wolfram's matching native constructors are source-inspected and pinned. Both independent reviews clear this connection at that stated level. |
| Original parser/comparator, exports and broader ledger prose | Not kernel-certified by this closure task. Historical computational evidence remains distinct. |
| General PDE solution completeness beyond the plane-wave ansatz | Not proved; requires an additional analytic setting and result. |
| Finite boundaries, interfaces, curved domains and weak solutions | Outside this pilot's setting. |
| Microscopic origin of the action, bulk shear-freeness and confinement | Supplied physical premises or separate work; not consequences of this action proof. |

The static longitudinal result does not remove that degree of freedom: its
restoring stiffness vanishes in the supplied action. A theorem about that action
does not independently justify choosing it.

The [proof report](RESULT.md) and original [verification record](VERIFICATION.txt)
give the precise pilot statements and its 21 selected audits. The new closure
instrument is [_measurements/S9_lean_contract_check.py](../../_measurements/S9_lean_contract_check.py).
It adds the scalar-phase audits and reproducible mutation records. The later
[S10 coverage map](../s10/COVERAGE.md) records the additional formal coverage
and the remaining integration work. The [combined checkpoint](../CHECKPOINT.md)
is historical; [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md) records this bounded S9
closure. No further Lean coverage is required by C1–C4.
