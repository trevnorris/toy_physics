# D2 odd-invariant dynamics: fidelity boundary

E1–E4 in [DYNAMICS_COVERAGE.md](DYNAMICS_COVERAGE.md), governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Local verification
passed; see [DYNAMICS_VERIFICATION.txt](DYNAMICS_VERIFICATION.txt). Both independent
fidelity reviews returned CLEAR; see
[DYNAMICS_FIDELITY_REVIEW.md](DYNAMICS_FIDELITY_REVIEW.md).

## Object and conventions

The native D2 Q9 V6 output is the actual row-major gradient polynomial
`P=(G11+G22)(G12-G21)`. `Gij=partial_i u_j`; the first index is the derivative
index. The already reviewed `S11Invariants.oddPairing` is precisely this form.
The new `density_identity` identifies the supplied density with `-beta*P/2`.
That coefficient is XFORM_EXTRA minus MAIN in `package_build`, not an inferred
normalization of a one-dimensional invariant space.

Beta is an arbitrary constant real coefficient with the supplied stiffness
units. It may be positive, negative or zero. This proves a conditional result
for the supplied displacement field and density, without identifying a new
physical field or fixing beta's value. Lean coordinates are `(t,x1,x2)`;
native functions have arguments `(x1,x2,t)`. The spatial index mapping is
`G j i = J j.succ i`, and the phase is `k·x-omega*t`.

Put `d=partial_1 u_1+partial_2 u_2` and
`c=partial_1 u_2-partial_2 u_1`. The actual derivative-defined momenta are

| Derivative coordinate | Momentum for u1 | Momentum for u2 |
|---|---|---|
| t | 0 | 0 |
| x1 | -beta c/2 | -beta d/2 |
| x2 | beta d/2 | -beta c/2 |

Lean's variational derivative is
`E=(beta/2)*(partial_1 c-partial_2 d, partial_1 d+partial_2 c)`.
`eulerLagrange_eq` proves this expression from those momenta. Its form does not
require commuting mixed derivatives. The native positive-divergence convention
is `E_native=-E`. The compact SymPy comparison simplifies mixed derivatives of
smooth functions; that simplification is CAS evidence, not a new Lean theorem.

## First variation and the precise non-null conclusion

`relativeAction` integrates the pointwise density change on R³ for a smooth
background and a smooth compactly supported variation. It does not subtract
infinite background actions. `relative_density_integrable` and the proved
quadratic expansion justify differentiation of the actual integral.
`relativeAction_deriv_eq_eulerLagrange` uses the S10 compact-support integration
by parts theorem. The fundamental lemma, continuity and full-support Lebesgue
measure give `actionStationary_iff_eulerLagrange`.

`witnessField=(cos x1,0)` is smooth and its E at the origin is `(0,-beta/2)`.
For every nonzero beta, `exists_nonzero_firstVariation` supplies an admissible
compact test field with nonzero first variation. `variationallyNull_iff` proves
that vanishing for **all** smooth backgrounds and compact tests occurs exactly
at beta=0. Some particular backgrounds may remain stationary for nonzero beta.
This is the bounded meaning of “not boundary-only”; no general classification
of divergence densities or finite-domain boundary conditions is supplied.

## Plane-wave mixing and exhaustive cases

For `r=(-k2,k1)`, the unaveraged modal action is
`-beta*(k·a)*(r·a)/2` and its amplitude gradient is
`M a = -beta*((r·a) k + (k·a) r)/2`.
Both the actual scalar derivative `modal_variation` and the independently
differentiated position-space identity `eulerLagrange_planeWave` identify M.
The latter holds for every real frequency, wavevector, amplitude and point.

`polarization_decomposition` covers every D2 amplitude when k is nonzero.
`modalOperator_longitudinal` and `modalOperator_transverse` prove
`M k = -beta*|k|²*r/2` and `M r = -beta*|k|²*k/2`.
These vectors need not be normalized: the cross dot products are both
`-beta*|k|⁴/2`. `mixing_iff` proves they are nonzero exactly when beta and k
are both nonzero. `mixing_cases` covers the following mutually incompatible
conditions explicitly:

1. beta=0, with arbitrary k: M is zero.
2. beta≠0 and k=0: M is zero.
3. beta≠0 and k≠0: both cross-sector matrix elements are nonzero.

The zero intersection belongs to case 1. At k=0 there is no nonzero directional
polarization frame. These statements concern the odd increment only; no claim
about the roots, positivity or kernels of the complete XFORM_EXTRA package is
made here.

## Compact native check

`_measurements/S11_lean_dynamics_source_check.py` AST-selects the original native
helpers for D2 Q9, package construction, V5 and the two modal routes. It does
not import the production driver, emit outputs or invoke any S11c job. It reads
P_D from `compute_q9(2)`, uses the actual XFORM_EXTRA-minus-MAIN action, and
selects the odd V5 combination using its actual coordinates in the native V1
RREF basis. It checks exact symbolic identities, not sampled equality.

| Native object | Correspondence with Lean |
|---|---|
| Q9 P_D | `oddPairing(spatialGradient J)` |
| XFORM_EXTRA minus MAIN | `lagrangian beta J = -beta P/2` |
| Native local equation | `-eulerLagrange` |
| Native route A (cosine stripped) | `-modalOperator` |
| Native period-averaged action | `modalAction/2` |
| Native route B Hessian | `modalOperator/2` |

These normalizations are recorded separately; equality of kernels would not
identify their signs or factors. The source check intentionally changes the
native action's sign and half-factor in memory and requires failed symbolic
identities. Its five explicit locus checks evaluate the Python translation of the Lean-side
target matrix, including positive/negative beta and the zero intersections.
They do not separately reevaluate the native routes at those points: native
correspondence comes from the exact symbolic identities. Exhaustiveness is
proved by Lean, not inferred from those finite examples.

The native SymPy `period_average` helper uses the algebraic substitution
`sin²→1/2`; the comparison checks that stated normalization and does not
formalize an integration routine inside Lean.

The Wolfram link is limited to source inspection of its actual Q9 V5,
`stiffnessBlueprint` XFORM_EXTRA coefficient, action subtraction, equation convention
and period average. Its source hash and exact anchors are recorded. No Wolfram
execution, production-export agreement, or comparator clearance is claimed.

## Reuse, review and stopping point

The implementation directly reuses the S10 coordinate-derivative, smoothness,
integration-by-parts and test-field lemmas. `S11OddDynamics.Variation` repeats
the action-specific relative-action argument for the odd density; it does not
apply a density-generic variation theorem. The completed S9/S10,
H1–H4 and I1–I4 records retain their historical hashes; the appended Lake library
declaration is part of the new snapshot rather than a rewrite of those records.

The local verification instrument records source/dependency hashes, imported
local olean hashes, sequential module builds, standard-axiom audits, isolated
mathematical mutations and their passing counterparts. Instrument failures and
timeouts are not accepted rejections. The user authorized the fixed packet's
transfer, both independent non-author reviews returned CLEAR, and their findings
are resolved in the closure record. E1–E4 is complete; higher-dimensional,
full-spectrum, interface and S11c work remains excluded.

Four S10 dependency olean hashes changed during the fresh sequential rebuild,
while their proof-source hashes remained unchanged. The E1–E4 report binds the
current objects to that build. Historical records retain their original object
hashes and are not represented as freshly rerun homogeneous/invariant reviews.
