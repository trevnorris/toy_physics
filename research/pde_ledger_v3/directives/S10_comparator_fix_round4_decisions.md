# The cross-engine comparator — fix round 4. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: `scripts/S10_cross_engine_comparator.py` · its regenerated `scripts/out/*.out` · a new directive
under `directives/`. ⛔ Out of scope: both engines, every `.wl`, every step record, every `.tex`.

⭐ Two legs reviewed round 3. ⭐⭐ **The constant/shadow repair holds and its mechanism is object-based, not a
list** — sized against the whole library it classifies hundreds of names correctly where the previous list
covered thirteen. ⭐ **The nullspace property still passes both halves on real payloads**, including a
symbolic invertible change of basis, and **round 2's repairs all still hold.**

⛔⛔ **This round fixes ONE defect and DELETES two guards. It adds nothing.**

⛔ Do not state, anywhere in the script or your directive, what any count, tally or population size comes
out as.

---

## F1 · ⛔⛔ A refuted pair must not classify a row as a naming difference. **This is blocking.**

⚠ **Measured on the REAL transcripts:** repoint shear stiffness to brane density in one emitted matrix, and
the run records it as a **naming-only** disagreement while the content population is unmoved. ⛔ The same
row prints that the pair is **transcript-refuted** and that the row is **naming-only**, together.

⚠ **Cause, and it is ordering, not a typo:** the disagreement kind is computed inside the per-row
comparison, and the refutation is not computed until after every comparison already exists. ⇒ ⛔ **the
classifier cannot see the refutation by construction.**

⭐ Round 3 got the **worklist** half right — a refuted pair is correctly refused from the rename worklist.
⛔ The **classification** half is what feeds the tallies a step record will cite, and a repointed physical
symbol is exactly what must not land in the naming population.

**Decision: a pair the transcripts refute supplies no naming evidence to any row's disagreement kind.**
⚠ Your own round-3 directive already says row-local substitution "does not alter a residual, guard, or
disagreement classification." ⭐ Make that true.

## F2 · Delete the two guard terms that cannot be nonzero.

⭐ **The content-decomposition residual**: the sub-populations are a complete, pairwise-disjoint case split
of the parent, so their sum equals it for **every** input.
⭐ **The shared-accounting residual**: the compared set is built as a disjoint union of exactly the pieces
it is then differenced against.

⚠ A leg attacked both with duplicated names, a one-sided duplicate, and a nullspace row whose matrix is
absent. ⛔ Neither can be made nonzero. ⇒ they read as measurements and carry no information.

**Decision: both go, and they come out of `FINAL_GUARD` with them.** ⛔ Do not replace them.
⭐ **Keep emitting the populations and the totals** — a reader can do the arithmetic, and that is where the
property was always established. ⚠ **This is my defect, not yours: I asked for both.**

---

## Record, ⛔ do not fix

1. ⛔ **The constant/shadow rule is applied to one engine only.** The other payload parser has no shadow
   protection, so a name that is a library callable there resolves to the library object while the first
   engine yields a free symbol. ⚠ Measured: the failure is **loud** — the row reports the residual as not
   computed and fails — ⛔ never a false agreement. ⭐ **But the standard spellings for two Greek letters are
   in that population**, so a future transcript using either gets **no cross-engine residual at all.**
   ⭐ Say which engine is protected and which is not.
2. ⛔ **Two sub-population selectors test the emitted NAME's suffix, not the residual.** ⭐ Measured clean on
   this input — every row so selected genuinely belongs there — ⛔ but a future transcript could route a row
   carrying real algebra away from the population that reports algebra.

## Constraints

- ⛔ Do not modify either engine or re-run them; the transcripts are committed inputs.
- ⛔ Do not start `wolframscript`. ⛔ No hardcoded name→name pair table, no config file, no runner.
- ⭐ Report the complete `.out` diff and account for every changed line.

## Acceptance

⛔⛔ **For F1, show it on the REAL transcripts**: repoint one physical symbol to another that both engines
distinguish elsewhere, and show the row leaves the naming population and enters the content population.
⭐ **And show a genuine spelling difference still classifies as naming-only** — ⚠ a fix that pushes
everything into content is no better than the defect.
⛔⛔ **Show the earlier rounds still hold**: the nullspace residual in both halves, the constant/shadow
split, derivative dependence sets, the computed disagreement kind. ⚠ **A round that quietly breaks an
earlier one has happened once already on this artifact.**
⭐ For F2, show by running that each deleted residual could not have been nonzero.

## Deliverables

The changed files, the literal stdout, the complete diff, and every ablation script with its stdout at
named absolute paths **outside the repository**.
