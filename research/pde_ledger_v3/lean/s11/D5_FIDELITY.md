# D5.1–D5.4 fidelity record

Local verification PASS. Both independent fidelity reviews are CLEAR; the
bounded D5.1–D5.4 contract is complete. See
[D5_FIDELITY_REVIEW.md](D5_FIDELITY_REVIEW.md) for dispositions and review limits.
Governed by [D5_COVERAGE.md](D5_COVERAGE.md) and
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

## Object and exhaustive proof

The object is a real homogeneous quadratic form on all real 5×5 matrices,
with the full conjugation action `G -> R G Rᵀ`. The gradient convention is
`G_ij=partial_i u_j`; the scalar coordinates are the 25 row-major entries. No
symmetric-gradient or transverse restriction is supplied. This increment
classifies densities, with no variation, integration or divergence quotient.
There is no additional dimensional rescaling or kinetic normalization.

The proof gives a unique coefficient representation in the basis `(tr G)^2`,
`tr(G^2)`, `tr(G Gᵀ)`, allowing arbitrary constant real coefficients. These are three distinct quadratic forms, not
three bulk responses. The proved census is SO/O/reflection-odd = 3/3/0.

`Quadratic` represents every actual `QuadraticForm` by its associated bilinear
form in the 25 coordinate matrices, giving all 325 unordered degree-two
monomials. This representation is a theorem, not an assumption about a
particular ansatz. Its `CoordinateAlgebra` and `BilinearExpansion` imports
separate finite definitions and checked row/tail identities into bounded
compilation units. `Rotation` proves the full-group predicates and the identity
that conjugation by `R` equals conjugation by `-R`. Orthogonal determinants are
±1; negation reverses their sign in dimension five. Thus SO and O invariance
coincide. This parity argument is separate from the completeness proof.

The independent finite generator selects 322 linearly independent necessary
equations from four adjacent quarter-turns and one rational plane rotation.
They are split into seventeen blocks of at most twenty equations. Each
equation is derived in Lean from quantified O invariance at an explicit
matrix. Seventeen reconstruction blocks retain all 322 coefficient targets
and their exact rational linear combinations; every necessary equation remains used.
The three remaining coordinates are coefficients 6, 29 and 25, giving
`a=c6/2`, `b=c29/2`, `c=c25`. `Forms` separately proves invariance under every
orthogonal matrix using trace identities; the selected rotations are never
used as a definition of full-group invariance.

`Classification` proves existence and uniqueness of coefficients and the zero
odd subspace. `Census` constructs an actual linear equivalence with `R^3` and
computes dimensions of the actual invariant submodules. `Controls` retains
zero/nonzero examples, three non-omission results, a non-invariant single-entry
quadratic, and nonzero form-normalization examples.

## Compact native link

`_measurements/S11_lean_d5_source_check.py` extracts only the original Q9,
coordinate and polynomial helpers. It calls the actual native `compute_q9(5)`
without importing the production module or running its driver. Its checks
compare the complete V1, V2 and V6 spans in the native monomial ordering,
the reflection operator, and the empty odd basis/zero `P_D`.

The same-count/wrong-span control replaces an invariant basis direction by
`G00^2`. A second control restores the historical coefficient-action orientation
in memory, leaving the engine file unchanged. The instrument must record the
observed counts and reject the wrong span; count agreement alone is insufficient.
Wolfram evidence is source inspection of `Transpose[actionRows]`, not execution.
The compact translation is outside Lean's kernel. Run2's native report confirms
the complete SO/O/odd spans 3/3/0, the reflection action and zero `P_D`.
Both wrong-span controls are detected; the old orientation retains the same
3/3/0 counts but fails span agreement. The native result is a tested artifact
connection, separate from the Lean proof and the independent fidelity reviews.

## Verification and controls

Recorded run9 built the first 43 canonical modules, including the complete
representation, all 322 necessary equations and their reconstruction,
full-group classification and census. Run10 reused those objects under the
full guards, then built `Controls` and the audit root. All 45 canonical objects
and 51 selected standard-axiom audits passed. Run11 reused those same objects
and passed all thirteen mathematical rejections and seventeen fresh positive
executions. Every local source, generator, instrument, dependency, input/output
object and recorded log matches the author validation.

The complete suite requires 75 check records: 44 module builds, one root build,
thirteen paired mathematical rejections and seventeen positive executions.
The 51 audits allow only `propext`, `Classical.choice` and `Quot.sound`, including
empty axiom sets. No admissions or custom physics axioms are accepted.

| Contract claim | Deliberately false paired statement |
|---|---|
| Dimensions and parity | SO dimension 4, O dimension 2, odd dimension 1; existence of a nonzero odd invariant; SO and O spaces unequal |
| Complete three-form span | Each of the three forms expressible without its own coefficient; `G00^2` invariant |
| Normalization | Trace-square witness 2 instead of 4; trace-of-square witness 1 instead of 2; Frobenius witness 2 instead of 1 |
| Reflection | Determinant +1 instead of −1 |

Each false statement must produce exactly one diagnostic in `contract_control`
with exactly one unsolved `False`. Its true counterpart must compile cleanly.
Four additional positives cover the zero form, a nonzero invariant, arbitrary
real coefficients and unique coefficients. These are paired statement controls;
no canonical source-replacement mutation is claimed. Compiler, syntax, import,
warning, timeout and memory failures are never mathematical rejections.

## Bounded proof repairs and recorded build lineage

All repairs preserved the supplied definitions, public classification claims,
325-entry coefficient witness and 322 necessary equations/reconstruction
certificate. They changed proof construction or compilation boundaries to stay
within one worker, `-j1 -M4096`, strict warnings and 600 seconds per process.
No resource or tactic limit was increased.

| Runs | Issue and resolution |
|---|---|
| 1–2 | Concrete reflection entries needed diagonal transpose/product reduction. The focused repair passed before the first formal suite. |
| 2–7 | A monolithic 625-term bilinear expansion and whole-sum normalization exceeded memory or simplifier limits. The final proof uses 24 bounded mixed rows, 25 tail identities and Mathlib `List.foldl1_eq_foldr1` for assembly. |
| 8 | Separate `CoordinateAlgebra`, `BilinearExpansion` and `Quadratic` builds passed, followed by all necessary-equation blocks. Combined reconstruction exhausted memory. |
| 9 | Seventeen reconstruction modules preserved the exact rational combinations and passed. Concrete vector/coercion reduction stopped `Controls`. |
| 10 | Explicit coordinate reduction fixed the canonical witnesses; all canonical builds/audits passed. The last paired mutant reduced to `oSpace ≠ oSpace`, failing the stricter requirement of literal `False`. The guard correctly refused it. |
| 11 | Only the paired SO/O control proof adds `ne_self_iff_false`. Both focused cases passed; canonical sources, generator and all 45 objects stayed unchanged. All thirty fresh control executions passed their intended outcomes. |

Detailed diagnostics, preflight comparisons and original sources remain in
`_scratch/S11_lean_d5/verification_run*/`. Focused checks never supply recorded
build reuse. Generator changes invalidated earlier reuse even when individual
Lean files were unchanged. Current reuse requires transitive local sources,
pins, generator, exact command, input objects and output object hashes. The
whole-process-group timeout cleanup is unchanged from its tested version.
An author validation attempt during the run10 repair overlapped an edit and
failed; it was discarded. Exact original bytes were restored and validated
before reapplying the repair. No failed validation was accepted.

## Preservation, trust boundary and completed reviews

The preservation manifest pins 383 historical proof/evidence files, two native
sources and 51 D4C/NP/T1/VC objects. D5 imports and rebuilds no earlier local
Lean module. Fifteen clean package pins and five direct Mathlib source/object
pairs are checked; the transitive Mathlib cache is the existing pinned baseline,
not a fresh rebuild in this increment. Shared lakefile/README additions belong
to D5. The separately validated installation/replay files remain unchanged.

The native span comparison and its translation are outside Lean's kernel.
No Wolfram execution, production driver, export regeneration or S11c change
is claimed. The proof concerns quadratic densities in dimension five before
any EL/divergence quotient, not D5 bulk responses or boundary effects.

Formal evidence is `_measurements/S11_lean_d5_contract_checks.json`; compact
native evidence is `_measurements/S11_lean_d5_source_checks.json`. The live-source
and build-lineage validation is `_measurements/S11_lean_d5_validation.json`;
[D5_VERIFICATION.txt](D5_VERIFICATION.txt) summarizes acceptance.

Claude and Grok independently cleared the user-approved fixed 64-file packet.
They reviewed statement fidelity and source conventions; neither independently
executed Lean/native checks or recomputed hashes. Author closure validation
rechecked all evidence and recorded only documentation changes. The frozen
snapshot, archive and transport remain unchanged.
No general-dimensional classification, general null-Lagrangian theorem,
interface/spectral/scattering/pole work or systematic CAS bridge is included.
This closes D5.1–D5.4 only. No new commit was requested or performed.
