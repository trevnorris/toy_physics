# D3 invariant classification — fidelity map

Scope: [J1–J4](D3_COVERAGE.md), authorized after checkpoint `9e66a534`.
Local verification and live source/object validation passed. Both independent
reviews returned CLEAR; J1–J4 is complete. See
[D3_VERIFICATION.txt](D3_VERIFICATION.txt) and
[D3_FIDELITY_REVIEW.md](D3_FIDELITY_REVIEW.md).

## Object and quantifiers

`S11D3Invariants.Quad` is mathlib's actual `QuadraticForm ℝ Mat`, where
`Mat = Matrix (Fin 3) (Fin 3) ℝ`. It is not a selected list of candidate
polynomials. `quadratic_representation` expands an arbitrary form through its
associated bilinear form into all 45 degree-two monomials in the nine entries.
The matrix-unit expansion is proved, so this representation does not assume
invariant-sector completeness.

Native Q9 and Lean both use `G_ij = ∂_i u_j`, row-major entries, and
`G ↦ R G Rᵀ`. Lean uses zero-based indices; native symbols are `g_1,…,g_9`.
`Orthogonal R` means `RᵀR=I`; `Proper R` adds `det R=1`. Both invariant
predicates quantify over every such R and every real matrix G.

This is a homogeneous degree-two density classification. There are no constant
or linear terms, field-equation restrictions, coefficient sign assumptions or
quotients by total divergences. The three coefficients have whatever units the
supplied density requires; no kinetic or elastic normalization is introduced.

## Basis and exact count

`invariantForm ![a,b,c]` is precisely

`a (tr G)² + b tr(G²) + c tr(G Gᵀ)`.

Here `tr(G²)` is the trace of a matrix product, **not** the square of the trace.
`tr(G Gᵀ)` is the sum of all nine squared entries. The source definitions and
application theorems identify these polynomials exactly, including every factor
of two on paired off-diagonal monomials. We do not claim that the native RREF
basis has these same rows; equality of their full spans is the required link.

The proved classification gives unique real coefficients for every SO(3)
invariant. Sufficiency under the full O(3) group follows from trace identities.
The actual invariant submodules therefore coincide, each with dimension three.
For `R₀=diag(-1,1,1)`, an SO-invariant form satisfying
`Q(R₀ G R₀ᵀ)=-Q(G)` must be zero. This identifies the odd submodule itself,
rather than inferring its size by subtracting two CAS counts.

## Finite certificate boundary

The independent proof generator uses two quarter-turn rotations and the
rational `(cos,sin)=(3/5,4/5)` rotation in a coordinate plane. Evaluating the
full invariance hypothesis on matrix units and pair sums yields 42 necessary
linear constraints on 45 coefficients. The generator checks their independence;
Lean does not need or separately prove that rank assertion. The emitted proof
derives each constraint from that hypothesis and checks explicit rational
linear combinations. Full-group sufficiency is proved separately; these are
not sampled witnesses offered as evidence of sufficiency.

Concrete test matrices use `Matrix.of` so matrix rewrites retain their intended
types. `conjugate_apply` gives the actual matrix-product entry formula. The
proof checks each image, then `polynomial_vec` expands the full quadratic
coordinate expression by definitional equality. These helpers address Lean
simplification behavior; they add no hypotheses or numerical assumptions.

The generator reads no native audit output. SymPy selects a finite algebraic
certificate, while Lean verifies its conclusions. The generator's own rank
assertion is not the completeness theorem: Lean checks the reconstruction of
every coefficient needed for classification. Its `--check` mode detects drift
between generated text and canonical sources.

## Actual native connection

`_measurements/S11_lean_d3_source_check.py` extracts five unchanged functions
from the native SymPy source and runs only `compute_q9(3)`. It checks the
row-major monomial order and equality of the entire V1/V2 RREF spans with the
three forms above, the 0×45 V6 basis, its reflection operator and zero `P_D`.
The explicit wrong-span control retains dimension three; an additional
in-memory control restores the old generator orientation. Neither may pass the
span comparison. Canonical native sources are not mutated.

This is a tested CAS identification, not a kernel-certified Python translator.
The correspondence between the three Lean trace formulas and the Python target
forms was independently inspected by both reviewers, including monomial order
and factors of two; there is no automatic extraction of Lean coefficients into
the Python check. Native V1 uses Lie-algebra generators, while the Lean theorem
quantifies over the full group. The compact check identifies the resulting
spaces by exact span equality, without formalizing the native algorithm.
The Wolfram connection is source inspection of its generator transpose only;
no Wolfram execution or dual-engine result is claimed for this increment.
Native source and instrument hashes are recorded. No production driver or
pinned export is read as an oracle or regenerated.

## Review and limits

The approved fixed 27-file packet includes source, contract, checks, controls
and provenance. Claude and Grok independently returned CLEAR with no blockers.
Both reviewed source and mathematics; neither rebuilt Lean, reran the checks or
recomputed hashes. The author separately revalidated the archive, snapshot,
transport, live sources, objects and dependency pins. The
[closure record](D3_FIDELITY_REVIEW.md) records identities, limits and optional
finding dispositions. No substantive changes or re-review were required.

No conclusion here concerns the complete dynamical spectrum, null-Lagrangian
equivalence classes, nonlinear mixing, interface response, or dimensions other
than three. The physical choice D=3 remains supplied by the ledger.
