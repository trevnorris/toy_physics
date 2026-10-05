# D2 invariant-space fidelity record

Status: I1–I4 complete. Local verification PASS; both independent reviews CLEAR.
See [INVARIANT_FIDELITY_REVIEW.md](INVARIANT_FIDELITY_REVIEW.md) for the reviewed
revision and finding dispositions. This record belongs to the separately authorized I1–I4
contract in [INVARIANT_NEXT.md](INVARIANT_NEXT.md), under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

## Intended object

The specification is [S11_SHARED_PHYSICS.md, Q9](../../directives/S11_SHARED_PHYSICS.md).
The ambient vector space is all homogeneous real quadratic forms in the four
entries of an arbitrary real 2×2 matrix `G`. The source convention is row-major
`g_1,g_2,g_3,g_4` in SymPy and `g1,g2,g3,g4` in Wolfram; the physical substitution
is `G_ij = partial_i u_j`. Entries are independent at the census stage. There
are no positivity, symmetry, tracelessness, nonzero-entry, or rank assumptions.
The domain has no dimensional units or material-parameter normalization here.
Equality means equality of quadratic polynomials on all matrices. It is not
an Euler–Lagrange or total-divergence quotient.

The action is `G ↦ R G Rᵀ`. SO invariance quantifies over all real `R` with
`RᵀR=I` and `det R=1`; O invariance drops the determinant restriction. The
reflection is `diag(-1,1)`, identical to the native Q9 choice. Reflection-odd
means eigenvalue -1 for this reflection **within the SO-invariant space**.

## Formal object and completeness route

[Quadratic.lean](S11Invariants/Quadratic.lean) uses Mathlib's actual
`QuadraticForm ℝ Mat`, rather than assuming a list of candidate monomials is
complete. The coordinates and inverse are

- `t=G11+G22`, `s=G12-G21`, `x=G11-G22`, `y=G12+G21`;
- `G11=(t+x)/2`, `G12=(s+y)/2`, `G21=(y-s)/2`, `G22=(t-x)/2`.

The associated bilinear form proves the exhaustive ten-coefficient
representation. Polynomial coefficients are unique. The general
`coefficient_action` identity proves that images of monomials stored in rows
act on polynomial coefficient columns through the transposed matrix.

[Rotation.lean](S11Invariants/Rotation.lean) proves every proper orthogonal
matrix has the form `[[a,-b],[b,a]]` with `a²+b²=1`. Under this action, `t,s`
are fixed and `(x,y)` transform by the double-angle matrix. Orthogonal matrices
have determinant ±1; the negative-determinant part is handled by composing
with the specified reflection. No sampled group replaces either quantifier. The native infinitesimal generator
`[[0,1],[-1,0]]` has the opposite sign to the derivative of this rotation
parameterization at the identity; multiplying a generator by -1 leaves its
zero constraint unchanged.

[Classification.lean](S11Invariants/Classification.lean) derives necessary
coefficient constraints using the proper rotations `(a,b)=(0,1)` and
`(3/5,4/5)`. These finite evaluations are used only for the upper bound.
The converse proves each remaining form invariant under **every** proper
rotation. The result is precisely

`A t² + B s² + C (x²+y²) + D t s`.

O invariance is equivalent to `D=0`. Reflection oddness within this space is
equivalent to `A=B=C=0`. [Census.lean](S11Invariants/Census.lean) gives linear
equivalences with `ℝ⁴`, `ℝ³`, and `ℝ`, respectively, and proves the even/odd
spaces span the SO space and have zero intersection. Thus the dimensions are
4, 3, and 1 and the decomposition is unique, including the zero form.
The recorded build and axiom audit pass; see
[INVARIANT_VERIFICATION.txt](INVARIANT_VERIFICATION.txt). Both independent source fidelity reviews returned CLEAR.

## Native discrepancy and bounded correction

At checkpoint `d80d02d2`, SymPy `compute_q9` stored each generator's monomial
images as rows and computed their right nullspace. This has the wrong
coefficient-action orientation. The frozen
probe report (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_q9_orientation_probe.json`)
records the original source SHA-256 and the concrete failure. Its companion
probe script is historical reconnaissance against that source revision; its
old expected-failure assertions should not be rerun against the repaired
source as a success test.

The minimal native change transposes **each generator block before stacking**.
No expected dimension or candidate basis is inserted into the engine. The
existing reflection difference block is diagonal in the monomial basis for
`diag(-1,1)` and is unaffected. Wolfram already uses `Transpose[actionRows]`;
its source is unchanged. This generic coefficient-action correction affects
all dimensions using this routine, but this increment measures and proves the
D2 census only. It does not certify D3–D5 production results.

The compact native check (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_invariant_source_check.py`)
executes five original helpers with their real coordinate convention, omitting
registry initialization and production execution. It compares the entire D2
V1, V2 and V6 row spaces to the four/three/one forms above using exact rational
RREF and rank identities. An in-memory wrong-orientation control retains
counts 4/3 but must fail the span comparison. The
record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_invariant_source_checks.json`)
contains the instrument and source hashes and exact polynomials. Its numeric
residual is an independent check of the explicitly supplied polynomial; it is
not extracted from the in-memory mutant basis. The frozen original probe
separately records that polynomial as an actual pre-repair native basis element.
The full-span comparison remains a separate test of the actual native output.

This is a tested SymPy-to-formal-definition correspondence, not a Lean proof
of the Python interpreter or SymPy's implementation. The Wolfram connection
is source inspection of the same action and block orientation; no fresh
Wolfram execution or cross-engine production clearance is claimed. The
automated Wolfram guard checks only for the transpose substring; both
independent reviewers inspected its actual per-generator `Table` placement.
That source inspection supplies the contextual fidelity link.
[Controls.lean](S11Invariants/Controls.lean) also expresses the actual wrong
polynomial `G11² + G12 G21 + G22²`, with the proper rotation and matrix giving
values 1 and 481/625, and the non-invariance conclusion.

## Closure and scope boundary

Local builds, 46 standard-axiom audits, nine mathematical mutations and eight
positive controls passed. Live sources, compiled outputs and the 15 dependency
checkout revisions match the recorded evidence. See
the validation record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_invariant_validation.json`).
Two independent non-author fidelity reviews inspected the approved fixed
packet and native correction and returned CLEAR with no blockers. All live
shared-source and proof hashes matched the reviewed snapshot before these
documentation-only closure updates. The bounded I1–I4 contract is complete.

The closed H1–H4 evidence remains historical at its reviewed snapshot; its
source hashes are not relabeled as current after this Q9 repair and the Lake
target addition. H1–H4 excluded Q9 and recorded MAIN's independence from it.
No homogeneous theorem source has changed. Broader production reruns, Q9 V5
Euler–Lagrange classification, PD-package consequences, D3–D5 completeness,
comparator/export reconciliation, and S11b/S11c calculations remain outside
I1–I4 and are not reported as completed here.
