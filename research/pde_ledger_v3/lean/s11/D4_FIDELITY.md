# D4 invariant classification: object and fidelity boundary

Fidelity record for [D4_COVERAGE.md](D4_COVERAGE.md). Canonical proof builds,
the axiom audit, all controls and compact native span identification passed.
Claude and Grok independently returned CLEAR with no blockers. The bounded
contract is complete; [D4_FIDELITY_REVIEW.md](D4_FIDELITY_REVIEW.md) records
review provenance, limits and optional finding dispositions.

## Supplied object and conventions

`Mat = Matrix (Fin 4) (Fin 4) ℝ`, with all 16 entries independent. The gradient
convention is `G_ij = partial_i u_j`; coordinates are row-major. No symmetric
gradient restriction or equation of motion is used. A quadratic density is a
`QuadraticForm ℝ Mat`, so constants and linear terms are outside this claim.
`SOInvariant` and `OInvariant` quantify over every real matrix R with
`RᵀR = I`, with the additional `det R = 1` requirement for SO. The action is
exactly `G ↦ R G Rᵀ`. Reflection is `diag(-1,1,1,1)`.

The coefficient order is `(a,b,c,d)` for `(tr G)²`, `tr(G²)`, `tr(G Gᵀ)`, P.
Coefficients range over all reals, including zero and negative values. The
normalization is not an energy prefactor: this increment classifies Q itself
and does not introduce `L = -Q/2` or any variation of Q. There are no field
regularity, boundary or wave-vector assumptions in this algebraic statement.

The orientation `(0,1,2,3)` fixes

```
P = (G01-G10)(G23-G32) - (G02-G20)(G13-G31)
    + (G03-G30)(G12-G21).
```

Lean proves the determinant identity
`P(R G Rᵀ) = det(R) P(G)` for every real R, without assuming invertibility.
Together with trace identities, this supplies full-group sufficiency. The
reflection-odd submodule means the minus eigenspace of this particular
reflection **inside** the SO-invariant quadratic forms; it does not include
arbitrary non-invariant odd forms.

## Completeness and the native connection

The associated bilinear form gives all 136 coefficients of an arbitrary
quadratic form on 16 entries. The generator selects 132 independent necessary
constraints from quarter turns in the XY, YZ and ZW planes and an XY rotation
with cosine/sine `(3/5,4/5)`. Lean proves each rotation is proper and each
selected equation follows from the full invariance hypothesis. Four blocks
of 33 equations feed the rational reconstruction in `invariant_polynomial`,
which reduces every invariant to the four displayed forms. Uniqueness and
explicit linear equivalences establish dimensions of the actual submodules,
not counts of a supplied candidate list. The generator is independent of the
native Q9 routine; its generated proofs have compiled in Lean. The proof uses
the checked equations and reconstruction, not an assumed generator rank.

`_measurements/S11_lean_d4_source_check.py` executes only the five selected
native helper definitions needed for D4 Q9, without running the production
driver. `_measurements/S11_lean_d4_source_checks.json` records source hashes
and exact complete-span comparisons:

| Native object | Contract identification |
|---|---|
| Q9 V1 | Span of the four displayed forms, dimension 4. |
| Q9 V2 | Span of the first three forms, dimension 3. |
| Q9 V6 | Span of P, dimension 1; its actual reflection operator is checked in the native V1 basis. |
| Native `PD_POLY` | Exactly P, with coefficient +1 on `G01 G23`. |
| Fully summed `epsilon_ijkl Gij Gkl` | Exactly `2 P`, with `epsilon_0123 = +1`. |

All span comparisons check exact row spaces, not just dimensions. The
same-count/wrong-span control replaces P with a single-entry square; the
historical wrong generator orientation is also tested only in memory.
Neither control edits the native source. Wolfram evidence here is source
inspection of the generator-block transpose, not a new engine execution.
The proof-to-native translation check is ordinary Python/SymPy evidence and
is not described as kernel-certified.

The governing `directives/S11_SHARED_PHYSICS.md` §7 defines native `P_D` as the
sum of the computed V6 basis in its emitted normalization. It does not define
`P_D` by a fixed epsilon formula. The check above follows that rule and records
the factor of two explicitly; no native output or coefficient is rescaled.

## Verification, review and limits

Verification run5 passed all 44 checks: thirteen module builds and the audit
root, fourteen mathematical rejections and sixteen positives. The 49 selected
axiom declarations use only `propext`, `Classical.choice` and `Quot.sound`.
Run5 reused run4's successful canonical builds with full local source,
pin/generator and input/output object guards; every control ran afresh.

The controls cover all four omission claims, the 4/3/1 census, reflection
restrictions, orientation normalization and sign, trace-form normalization,
and a non-invariant single-entry square. Twelve false statements reduce to
`False`. The two source mutations leave explicitly false polynomial identities:
one changes a trace-square cross coefficient from 2 to 3; the other changes a
trace-of-square cross coefficient from 2 to -2. Their canonical identities and
all sixteen positive controls pass. No timeout or instrument failure is
accepted as mathematical evidence.

`_measurements/S11_lean_d4_validation.json` records correspondence of all 44
logs, fourteen compiled D4 objects, proof sources, generator, instruments,
native sources and pins, including clean tracked source at the fifteen pinned
dependency-package revisions. [D4_VERIFICATION.txt](D4_VERIFICATION.txt)
records the bounded tactic repairs, certificate splitting and tested timeout
cleanup. They changed proof execution, not the theorem, necessary equations
or reconstruction. The final determinant-sign control repair changed only
rewrite order; all canonical sources and objects remained unchanged.

The user authorized the fixed packet with “Proceed”. Its repository manifest,
`_measurements/S11_lean_d4_review_packet.json`, is byte-identical to the
`MANIFEST.json` supplied to both reviewers. Both completed independent source
fidelity reviews; neither independently compiled Lean, reran instruments or
recomputed hashes. The author validated the archive/snapshot/transport and
live proof/evidence correspondence separately. Only closure documentation
changed afterward; no proof or instrument change required re-review. Completed
H/I/E/J/K proof files and historical evidence remain untouched.

This increment does not prove P is a divergence or has zero bulk variation.
That is the next separate contract. D5, general null Lagrangians, spectra,
interfaces, S11c calculations/exports and systematic CAS bridging are excluded.
