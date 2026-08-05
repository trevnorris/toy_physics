# S10 Mathematica engine — repair round 1

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
**Its specification, still binding:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`

⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**
⛔ **Do not read** `research/pde_ledger_v3/steps/` or `research/pde_ledger_v3/paper/`, and ⛔ do not reach
any earlier version of this engine via `git log` / `git show`.

⭐ Two independent review legs ablated the committed engine. **The core method survived both**: `N3`, `N7`,
the separately-constructed matrix routes, and per-package re-entry at the action are all correct — ⛔ **do
not "improve" any of them.** What follows are the defects.

⚠ After every fix, **re-run the engine** and confirm it still completes under **10 minutes** and exits `0`.

---

## ⛔⛔ R1 — BLOCKING: the matrix used downstream is built under an assumption about the very sign a later test determines

`periodAverage` (≈246–257) injects **`omegaSquared > 0`**, plus `GenerateConditions -> False`, into the
construction of `M_B` — and `M_B` is the matrix **everything downstream uses**, including the Q3 test whose
job is to determine the sign of each root.

⛔⛔ **A test cannot determine a sign that was assumed while building the object under test.** ⚠ Any root
this assumption excludes is excluded **silently and by construction**, and a closed-form result looks
equally healthy either way.

⭐ **Fix:** the period average must be taken with **only the §3 joint assumption set**, which ⛔ says nothing
about the sign of `omegaSquared`. If the averaging integral genuinely needs a nonzero-frequency condition
to evaluate, ⭐ **carry it as an emitted condition tag** rather than as a silent assumption, and ⛔ do not
let it reach the sign test.

⚠ **Related, and it is why this matters rather than being cosmetic:** because the spectrum is taken from
`M_B` alone and the `M_A − M_B` residual is coding-consistency-only, an error introduced in the
period-average path **flows into the printed spectrum with nothing cross-checking it.**

## ⛔⛔ R2 — BLOCKING: a nonzero exit on a search outcome

`Quit[2]` (≈1073) is reached whenever `Reduce` returns anything but literal `False` and `FindInstance`
then finds no point. ⚠ An ablation drove this and the run **died with 75% of the sweep unemitted.**

⛔ **A physics finding must EXIT 0.** Only an operational failure — an exception, an unevaluable
expression — may exit nonzero. ⭐ **Fix:** emit the failure to locate a point **as its own tag**, with the
`Reduce` output and the `FindInstance` call that failed, and **continue the sweep.**

⚠ Note also that `intersects = Not[SameQ[allowed, False]]` (≈993) is a **syntactic** test on `Reduce`'s
output, so an *undecided* `Reduce` is treated as an intersection. ⭐ Distinguish **decided-nonempty**,
**decided-empty**, and **undecided**, and emit which one occurred.

## ⛔⛔ R3 — BLOCKING: the dimension walker returns a clean, wrong answer for anything it does not recognise

The walk's fallback (≈639–644) assigns `{0,0,0}` to any unhandled head. ⇒ an unrecognised structure is
**silently dimensionless** and its tag is indistinguishable from a real dimension.

⭐ **Fix:** an unhandled head must produce an **explicit indeterminate marker** that ⛔ cannot be mistaken
for a computed dimension, and the walk must **emit the head it failed on**. ⚠ A wrong dimension that looked
like a meaningful signal is the exact defect that survived two review legs and a full ablation suite at the
previous step.

## ⛔⛔ R4 — BLOCKING: the homogeneity check is evaluated under the solution it was derived from

The dimension equations are formed by setting each action term equal to the energy-density target, solved
for the coefficient dimensions, and then the homogeneity constraints are evaluated **under that same
solution** (≈652–668, 760–766). ⇒ it is true by construction.

⚠ **Measured:** an ablation corrupting the dimension walk moved **59 dimension tags** and moved
**no `_HOMOGENEITY` tag at all.**

⭐ **Fix:** homogeneity must compare **independently obtained** per-term dimension vectors — the ones the
tree walk produces — against each other, ⛔ not against the constraint that defined them. ⭐ Emit the
per-term vectors **and** their pairwise differences, so a failure is readable.

## ⛔ R5 — a check that cannot fail

`Q6_*_COEFFICIENT_HOMOGENEITY`: a bare symbol has no addends, so the test returns `True` for **every
package and every ablation**. ⭐ Either make it compare something that can differ, or ⛔ **delete it** — an
unfailable check is worse than none, because it reads as coverage.

## ⛔⛔ R6 — the record of what ran is a restatement of the declaration

`WL_S10_RUN_PAIRS` and `WL_S10_SKIPPED_PAIRS` (≈1209–1213) are emitted **before** the sweep, and
`skippedRunPairs = Complement[declaredRunPairs, declaredRunPairs] ≡ {}`. ⚠ In an ablation that died inside
the **first** package, both tags reported a **complete sweep**.

⭐ **Fix:** accumulate the `(package, D)` pairs **as each one actually completes**, and emit both tags
**after** the sweep. ⭐ `SKIPPED` is then `declared ∖ completed`, computed from two genuinely different
lists.

## ⛔ R7 — emission conditional on a payload's value

`_Q5_SCALING_EXPONENT` (≈529–531) is emitted only when `purePowerTest` holds. ⛔ Corollary 4: a **missing**
tag is indistinguishable from *never computed*. ⭐ **Fix:** always emit it; where the ratio is not a pure
power, emit an explicit not-a-pure-power marker as the payload.

## ⛔⛔ R8 — 314 of 510 tags can vanish at exit 0

When the root list is empty, all of Q4 `N1`–`N7`, Q5, and the per-root Q6/Q8 tags are **absent**, with no
tag recording that the list was empty and a **clean exit 0**. ⭐ §7's strata path already does this
correctly (≈1048) — ⭐ **apply the same treatment to roots.**

## ⛔ R9 — the root count can silently under-report

`Select[..., FreeQ[#, omegaSquared] &]` (≈337–339) discards solutions `Solve` returned implicitly, with no
tag. ⇒ `_Q3_ROOT_COUNT` can mean less than its name. ⭐ **Fix:** emit the discarded set as its own tag, and
emit the count **before and after** the filter.

## ⛔ R10 — solver conditions are emitted and then discarded

`removeConditionalWrappers` (≈314–315) strips conditions before Q4/Q5/Q6/Q8 consume the roots, and they are
**never re-imposed**. ⭐ **Fix:** re-impose each stripped condition where its root is used, or ⭐ emit a tag
per root recording the condition that was dropped and that it was not re-imposed.

## ⛔ R11 — Q8b is incomplete and collides with the generic tag shape

Strata reuse `_Q3_ROOT_COINCIDENCE_LOCI` for a `{pair, boolean}` list where the generic path carries a
real-domain locus, and the four `_ROOTp_ROOTq_Q3_COINCIDENCE_*` tags are omitted on strata.
⭐ **Fix:** give stratum tags their own distinct names, and emit the **full** Q3 tag set on each stratum —
Q8b exists so that the stratum is analysed exactly as the generic case is.

## ⛔ R12 — the ansatz reintroduces the branch §3 removed

`Sqrt[omegaSquared]` (≈238–239, 249, 271) picks the positive root. §3 requires `ω²` as a **single symbol**
precisely so no branch is chosen. ⚠ Harmless to this run's algebra, but the sibling engine may resolve it
differently ⇒ **a cross-engine divergence for no physical reason.**
⭐ **Fix:** keep the ansatz in terms of `omegaSquared` throughout, or ⭐ emit the branch choice as an
explicit tag so the divergence is visible rather than silent.

## ⛔ R13 — the brane dimension in Q6 is unconstrained

`braneDimension` is never tied to the run's integer `D` outside `runRankStrata` (≈1065), so the Q6
coefficient dimensions are closed functions of a symbol **nothing else in the run constrains**.
⭐ **Fix:** emit the relation between `braneDimension` and the package's concrete `D` explicitly, and ⭐ emit
the Q6 result **both** symbolically in `D` **and** specialised to the run's `D`.

## ⛔ R14 — distinct strata are merged by their printed form

`GatherBy` on the printed form of `Allowed` (≈1034–1037) groups strata by **string equality** of `Reduce`'s
output and keeps only the first group's equations. ⇒ two genuinely distinct strata that print alike yield
**one** point. ⭐ **Fix:** group by a **structural** comparison, ⛔ never by printed form.

## ⛔ R15 — Q6's audit set omits several emitted expressions

Dimensions are not audited for `N1`, `N5`, `N6`, the Q8 minors, or Q7. ⭐ **Fix:** extend the audit set to
every emitted dimensionful expression, as Q6 requires.

## ⛔ R16 — a hand-typed list presented as a record

`SYMBOLIC_SIMPLIFIERS` (≈1124–1126) is a hand-typed list that would print identically if none of those
simplifiers were used. ⭐ Either record what was **actually invoked**, or ⛔ delete the tag.

---

## ⭐ What NOT to change

⛔ `N3` (the stacked-rank transverse count), `N7` (the dual-algorithm count residual), the separate
construction of `M_A` and `M_B`, per-package re-entry at the action, Q5's scaling method, and Q7's
independent gradient symbols. ⭐ All four survived one-sided corruption and form ablation on both legs.

## Report back — ⛔ under 25 lines

1. One line per `R1`–`R16`: fixed / partially fixed / not fixed, and the line numbers touched.
2. New tag count, wall-clock runtime, exit code.
3. ⛔ **Do not report what any value came out to be**, and ⛔ do not say whether anything "worked".
4. ⭐ Anything in this list you believe is **wrong**, or any fix that would break something the review legs
   confirmed is correct. ⭐ This is wanted.
