# The cross-engine comparator — fix round 1. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: `scripts/S10_cross_engine_comparator.py` · its regenerated `scripts/out/*.out` ·
`directives/S10_cross_engine_comparator_directive.md` (⭐ **the exclusion decision lives there — fix it at
that level too**) · a new directive under `directives/`.
⛔ Out of scope: both engines, every `.wl`, every step record, every `.tex`.

⭐ **Two legs reviewed your build.** ⭐⭐ **The central property HOLDS**: an independent reimplementation found
**no row where you report agreement and it does not** — all 362 confirmed, every divergence conservative,
and every one of the 562 shared names carries exactly one verdict. ⛔ Both items below are about
**falsification power**, ⛔ not about false agreement.

---

## F1 · ⛔⛔ A non-canonical object still has to be compared. Excluding it absorbs wrong physics.

Any name ending `N6_NULLSPACE_BASIS` is dropped from value comparison as "defined only up to basis change."
⚠ **Measured — three mutations, and NOTHING moves**: substituting root 2's basis for root 1's; scaling one
component `k1/k2 → 2*k1/k2`; replacing `{k1/k3,k2/k3,1}` with `{1,0,0}`. ⛔ Not residual, not guard, not a
counter, not the exit code — on the third, the **only** changed line in the whole stdout is the echoed
operand.

⛔ **None of those is a basis change.** Measured against the emitted `N1_MATRIX`: `M·v` is `[0,0]` for the
emitted basis and for a legitimate `7×` rescale, and **nonzero for all three mutations** — the substituted
vectors are **not in the nullspace at all**.

⛔ **And the retained rows do not entail the excluded one.** Rank and nullity fix the **dimension** of the
nullspace, ⛔ never **which** subspace; under the third mutation `N2_NULLITY`, `N7_BASIS_COUNT` and every
other retained row still agree. ⇒ this is **exclude-for-reading-zero, ⛔ not exclude-for-entailment.**
⚠ It also discards real agreement: **15 of the 26** excluded rows agree exactly by the leg's own route.

**Decision: a non-canonical object is compared by a residual that is invariant under exactly the freedom it
has, and sensitive to everything else.** ⭐ Both engines already emit what that needs. ⛔ If you conclude no
such residual exists for some object, **name it and stop** — ⛔ do not restore a blanket exclusion.

⭐ **Fix the directive too.** The exclusion is stated there as a decision; a repaired script under an
unrepaired directive leaves the next builder the same instruction.

## F2 · The disagreement count is not a physics count, and nothing says so.

**164 disagreements decompose, measured:** ⭐ **23 naming-only** (`D↔braneDimension` on 13 rows,
`s↔coefficientScale` on 10 — the physics agrees), ⭐ **13 representational** (`Q1_EULER_LAGRANGE_SYSTEM`:
py a Matrix of expressions, wl a list of `expr == 0`), ⭐ **128 genuine content divergence**.
⛔ No printed counter separates them, so the step record cannot cite the number that matters.

**Decision: emit the decomposition.** ⛔ Do not state what any of the counts comes out as; emit them and let
them be read.

⚠ ⭐ **One entry in your D12 worklist is mis-typed and a naming pass would act on it wrongly:**
`PY=M_B / WL=quadraticFormRoute` is a **route-tag divergence, ⛔ not a symbol spelling** — consuming it as a
rename would bind two different descriptions. ⭐ Also: the `s` pair inference missed
`S10_XCOEF_SCALE_D3_Q3_DETERMINANT` (conservative, no verdict affected).

---

## Record, ⛔ do not fix

⭐ **The baseline `FINAL_GUARD` is already `FAIL` and the exit code is 1**, so ⛔ neither carries any ablation
signal — only per-row residuals do. ⭐ **Say so in your directive**, because a future regression bar keyed on
the exit code would be blind to every mutation above.
⭐ Of the 362 agreements, **215 are bare integers and 14 are empty containers**; ⭐ **133 are symbolic or
structured.** ⭐ State that — it is what the step record may claim.

## Constraints

- ⛔ Do not modify either engine or re-run them; the two transcripts are committed inputs.
- ⛔ Do not start `wolframscript`.
- ⛔ No hardcoded name→name pair table, no config file, no per-step file, no runner.
- ⭐ Report the complete `.out` diff.

## Acceptance

⛔⛔ **Each of the three mutations above must move a residual.** Show it, by running.
⛔⛔ **And a legitimate basis change or rescale must NOT move it.** ⚠ A check that fires on both is no better
than one that fires on neither.
⛔ **A test that passes on a weaker fix is not a test.**

## Deliverables

The changed files, the literal stdout, the complete diff, and every ablation script with its stdout at
named absolute paths **outside the repository**.
