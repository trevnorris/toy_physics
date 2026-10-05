# Homogeneous S11: statement fidelity and compact source connection

Author: Codex. This record identifies the finite H1–H4 contract in
[COVERAGE.md](COVERAGE.md). Local verification passed and both independent
fidelity reviews returned CLEAR. [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md)
records closure and dispositions. This completes H1–H4, not the whole S11 step.

## Supplied action, coordinates and retained setting

The actual selected action is

`L = rho/2 |u_t|² - mu/2 S_curl(du) - B/2 (div u)²`,

where `S_curl(J) = (1/2) sum_ij (J i.succ j - J j.succ i)²` and
`div(J) = sum_i J i.succ i`. In three dimensions `S_curl = |curl u|²` by the
existing S10/S9 specialization. This is the homogeneous MAIN action in
[the shared physics](../../directives/S11_SHARED_PHYSICS.md), not a derivation
of the material constants or of the physical choice of fields.

The Lean jet uses time row 0 and spatial rows `i.succ`; its displacement
component is the second index. The native SymPy constructor's `G[i,j]` is
`partial_(x_i) u_j`, and its velocity placeholders are `partial_t u_j`.
Its position-space functions list spatial arguments before time, whereas the
Lean spacetime vector lists time first. The identification uses named derivatives,
not an assumed common argument order. Both real plane-wave phases are
`k·x - omega t`.

Map `rho_br`/`rhoBr` to `rho`, `mu_R`/`muR` to `mu`, and
`B_comp`/`bComp` to `B`. The physical domain is D=3, positive constant
`rho,mu,B`, and nonzero real `k`; frequencies and amplitudes are real. The
integrated theorem allows smooth nondecaying backgrounds and uses smooth
compact variations on all spacetime. Finite background total action is not
required. Boundary conditions and variable coefficients are not supplied here.
The `B=0` and `k=0` results are separately stated boundary cases.

The displacement has length units; `rho` is inertia per D-volume and `mu,B`
share energy-per-D-volume units. Thus both squared-frequency branches scale as
inverse time squared, and `B/rho` has squared-speed units. These unit and
material identifications are inherited physical premises, not fresh dimensional
theorems or a certification of all dimensional output fields.

## Action proof and native matrix normalization

[Action.lean](S11Homogeneous/Action.lean) defines the unreduced combined density.
`lagrangian_split` identifies it exactly with S10's curl density plus its
divergence-only density with **zero inertia**. The existing variation theorems
then differentiate the actual sum. Kinetic energy is counted once.

`relativeAction_split` proves the same identity for the finite integrated
density change, using the existing integrability results before splitting the
integral. `relativeAction_deriv_eq_eulerLagrange` differentiates that integral;
`actionStationary_iff_eulerLagrange` uses compact tests and smoothness to obtain
the local PDE. The local expression is the sum of the two already derived
Euler–Lagrange expressions, not a separately supplied field equation.
`actionStationary_planeWave_iff` connects arbitrary compact-test stationarity
to the full modal kernel.

Write `K=|k|²` and `z=omega²`. Lean's action-derived operator is

`E = (rho z - mu K) I + (mu - B) k kᵀ`.

The [native SymPy builder](../../scripts/S11_stray_longitudinal_sympy_audit.py)
constructs MAIN in `package_build`, using `stiffness_densities`. Its
`route_a_matrix` follows the opposite EL sign convention from Lean and yields
`M_A = -E`. Its `route_b_matrix` differentiates the phase-averaged density
twice and yields `M_B = E/2`. The selected downstream matrix is `M_B`.
The [Wolfram builder](../../mathematica/S11_stray_longitudinal_mathematica_audit.wl)
uses the same MAIN terms in `stiffnessBlueprint` and `kineticRecords`, the same
opposite EL convention in `eulerExpressions`, and the same normalized cosine
phase integral in `averagedLagrangian` before constructing `matrixB`.

These signs and factors are retained, not discarded by comparing eigenvalues.
Inside Lean, `phaseDensity_eq` and `phaseAverage_eq` reuse the two existing
phase-density proofs and the actual normalized integral over `0..2 pi`.
The resulting average is half the modal action. At the handwritten-reference
boundary, `lagrangian_split` and the existing `S10Pilot.modalAction_eq` and
`S10Controls.modalAction_eq` expand that action to

`[(rho z - mu K) |a|² + (mu - B) (k·a)²] / 2 = dot a (E a) / 2`.

This is the algebraic identification used by the checker and inspected by both
reviewers; it is not a separate new S11 theorem. SymPy's `period_average` uses
the rewrite `sin(phase)² -> 1/2` on this quadratic density. Lean and Wolfram
use the normalized integral. For MAIN the density is proportional to sin²,
so the rewrite gives that integral's value, checked by an exact symbolic
comparison. Multiplying E by the nonzero
scalars -1 or 1/2 preserves its kernel; it would not preserve resolvent or
residue normalization. No such spectral-response equivalence is claimed.

The compact checker (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_source_check.py`) executes
only selected original SymPy definitions. It reads the needed inherited scalar
values from their actual records in `S10_exports.py`; it does not execute that
whole export file or a CAS audit. An explicit unused sentinel fills the invariant
census input expected by the generic constructor, and the checker verifies that
MAIN does not depend on it. No invariant census is performed or certified.

Exact symbolic zero residuals in the source-check record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_source_checks.json`)
identify the action, both matrices, phase-average normalization and the native
D3 determinant of the selected matrix `M_B`,
`(rho z-mu K)² (rho z-B K)/8`. Reversing the compression sign
changes both the density and operator; omitting the phase-average half is also
detected. The reference formulas are handwritten translations of the Lean
definitions. The checker does not parse Lean or prove its own translation.

The Wolfram leg is source inspection, not a fresh run. Seven literal anchors
check selected constructor and normalization lines; the entire file is hashed
for the packet revision. The anchors do not guard every line. All source hashes
are provenance, not automatic clearance of future edits. Both reviewers
inspected the current curl-normalization body as well as the anchored lines;
the anchors do not separately guard that body. The completed reviews and
comparison with the fixed packet establish fidelity at this revision.
Existing full-audit, parser, registry and component-census debts are not waived.

The recorded `route_A_stripped_phase` is the original SymPy routine's convention
label, not an independently computed phase factor. The actual phase identity
is proved by `eulerLagrange_planeWave`; the native matrix is checked separately.
Likewise, the source record's PASS and sentinel-absence flags summarize the
preceding successful assertions; those flags alone supply no extra evidence.

## Exhaustive kernel contract

For `rho,mu,B > 0`, `k ≠ 0`, set `T=mu K/rho` and `L=B K/rho`.
The general classification in [Spectrum.lean](S11Homogeneous/Spectrum.lean)
has no genericity or coordinate-chart assumption.

| Case | Entire kernel | D3 dimension |
|---|---|---|
| `B ≠ mu`, `omega²=T` | `transverseSpace k` | 2 |
| `B ≠ mu`, `omega²=L` | `longitudinalSpace k` | 1 |
| `B=mu`, `omega²=T=L` | the entire amplitude space | 3 |
| frequency squared unequal to both T and L | `{0}` | 0 |

`frequency_coincidence_iff` proves exactly `B=mu` on the nonzero-wavevector,
nonzero-inertia domain. The conditional expression in `modeSpace_classification`
is exhaustive and includes the merged case; `kernel_census_three` counts those
full spaces. These are geometric dimensions, not numbers of null vectors sampled
or numbers of distinct roots. The native determinant factorization separately
gives algebraic multiplicities two and one, or three at coincidence. This
contract does not add a Lean polynomial-root-multiplicity theorem.

`positive_frequencies` supplies positive square roots of both branches. Both
frequency signs are included by the classification in `omega²`. The zero vector
belongs to every kernel but is not counted as a propagating mode. `transverse_operator`
holds for every B: compression vanishes on the entire transverse space.
`longitudinal_operator` is independent of mu. No inference about these properties
is drawn from a single vector or a potentially mixed basis at coincidence.

`zero_compression_action`, `zero_compression_operator` and
`zero_compression_stationarity` identify B=0 with the existing S10 model exactly.
Its static longitudinal line is therefore recovered. `zero_wavevector_iff`
separately says that at k=0 stationarity holds iff `omega=0` or the amplitude
is zero; the nonzero-k dimension formula is not used there.

## Kinematic threshold

[Threshold.lean](S11Homogeneous/Threshold.lean) uses the supplied bulk sound
relation. On the longitudinal branch, `phase_matching` identifies

`k_w² = K (B/(rho c_s²)-1) = omega²/c_s²-K`.

This premise is Q11 of `directives/S11_SHARED_PHYSICS.md` (line 880 onward),
with `c_s0`/`cs0` mapped to Lean's `cs`. The original SymPy `q11_objects`
(line 1672 onward) forms `root = c_s0² (K + kwSquared)` and solves for the
normal square. The Wolfram bulk record (lines 2064–2091) supplies the same
dispersion and substitutes each root before solving. Substituting the proven
longitudinal root `B K/rho` gives the displayed formula. This is source
inspection of Q11, not execution by the compact H1 checker or a fresh native
Q11 run. The supplied section explicitly lacks an interface condition or
interaction operator; none is added here.

`threshold_classification` proves all three sign equivalences for positive
rho and c_s and nonzero k: below `B=rho c_s²` the value is negative; at the
locus it is zero; above it it is positive. These are evanescent, grazing and
propagating normal-wavevector possibilities. Neither existence of a propagating
channel nor evanescence establishes actual coupling, leakage or a bound state.
Those require the separate interface and spectral-boundary problem.

## Mutation map and stopping boundary

The verification instrument (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_contract_check.py`)
runs one Lean process at a time with `-j1 -M4096`, warnings as errors, and low
scheduling priority. It records commands, source and output hashes, generated
control statements and diagnostics. Source-mutation rejections must expose a
mathematical identity in the named declaration; concrete false statements must
reduce to `False`. Syntax/import/time/resource failures do not count.

| Load-bearing claim | Control |
|---|---|
| Correct compression density and inherited variational combination | Reverse its sign in the actual density, retaining `lagrangian_split` and the derivative proofs. The unmodified module is the positive control; the native checker independently detects both density and matrix changes. |
| Correct phase-average normalization | Change the proved half-factor to one, retaining the actual normalized integral. The unmodified theorem and native averaged-action residual provide passing evidence. |
| Complete separate and coincident kernels | Three concrete count pairs assert 2/1/3 versus false 3/2/2, using the general kernel-census theorem at an admissible nonzero wavevector. |
| Compression leaves the transverse branch unchanged | A transverse mode at B=9 remains stationary on the shear cone; its negation must fail. The universal result is `transverse_operator`. |
| Longitudinal frequency depends on B, not mu | At rho=mu=1, B=4, a longitudinal amplitude is not stationary at omega=1. Its false stationary counterpart must fail. The on-branch count control uses omega=2. |
| No off-root mode | At omega=3 with roots squared 1 and 4, a nonzero transverse amplitude is not stationary; the contrary statement must fail. |
| B=0 recovery | A nonzero static longitudinal mode passes; its negation must fail. The action and integrated-stationarity identities establish the general limit. |
| k=0 boundary | A nonzero static amplitude is stationary and a nonzero-frequency amplitude is not. Both have contrary-statement mutants. |
| Three threshold signs and k-domain | Explicit below/on/above-locus examples and opposite sign/equality statements; at k=0 the normal square is zero even away from the coefficient locus, so the nonzero-k guard is essential. |

The validated [verification record](VERIFICATION.txt) contains three module
builds, 29 standard-axiom audits, 15 mathematical rejections and 13 explicit
positive controls. The two source mutants fail on false identities: the sign
mutant would equate `B (div J)²/2` with its negative (take B=1, div J=1);
the phase mutant would equate the modal action with twice itself (take
rho=mu=1, B=4, omega=2, k=(0,0,1), a=(1,0,0), giving modal action 3/2).
Every concrete false-statement mutant leaves the goal `False` in its named
control declaration. These are mathematical failures, not instrument errors.
The verified proof/dependency hashes and source-check hashes still match the
live files. Both independent fidelity reviews are complete, and H4 closes with
the dispositions in [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md).
The historical “pure trace” sentence is corrected in the S11 step record:
longitudinal rank-one gradients are symmetric and curl-free but generally have
a nonzero symmetric-traceless part. This clarification neither changes the
selected action nor establishes uniqueness of an invariant basis.

SO(D)/O(D) invariant counting, total-divergence quotients, other action packages,
S11b interfaces and S11c scattering are excluded. No files or jobs in the active
S11c calculation chain are modified by this contract. Stop after H1–H4; do not
grow a transcript-wide CAS bridge or infer a confinement theorem.
