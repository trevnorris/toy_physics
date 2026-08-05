# HOW WE WORK — ⛔ read before touching anything

⭐ **Fifteen rules. Every one exists because ignoring it cost a session.**
⭐ Full process: `docs/development_pipeline.md` · what we're building: `docs/development_plan.md` ·
you-are-here: `STATUS.md`.

---

## ⭐⭐⭐ THE METHOD — ⛔ the CAS answers, ⛔ not me

1. ⭐⭐ **TWO ENGINES EXIST SO THEY CAN DISAGREE.** Independent *construction*, ⛔ not hidden answers.
   ⭐ **The disagreement IS the measurement.**
2. ⭐⭐ **A SCRIPT PRINTS COMPUTED OBJECTS. ⛔ IT NEVER STATES CONCLUSIONS.** Emit **both operands and the
   residual**, ⭐ *then* guard. ⛔ A residual asserted zero always prints `0` and carries no information.
   ⭐ Interpretation belongs to the **step record**.
3. ⭐⭐⭐ **NAME THE OBJECT. ⛔ DO NOT SPECIFY THE RECIPE.** ⚠ If a review is arguing about *how* to derive
   something — is this quotient well-defined, is this weight unique — ⛔ **the question was manufactured by
   specifying a derivation path.** ⭐ Ask for the object; ⭐ let the engine hand over what it built.
4. ⭐⭐⭐ **IF I AM DECIDING IN PROSE WHAT THE ENGINES SHOULD COMPUTE, I HAVE INVERTED THE METHOD.**
   ⚠ **The tell:** many turns reasoning toward an answer that a script would settle in one.
   ⚠ **The cause, measured:** the instrument was broken, so no measurement was available — ⇒ ⭐ **FIX THE
   INSTRUMENT**, ⛔ do not reason around it.
5. ⭐ **The spec says what to COMPUTE. ⛔ Never what anything EQUALS**, is expected, or was measured.
   ⛔ Withhold exactly one thing: **an acceptance criterion referencing an expected value** — a builder
   iterating to exit 0 converges on any target it can see.
6. ⭐ **A DISAGREEMENT IS A FINDING.** ⛔ Do not try to make divergence impossible with more careful prose
   — ⚠ that defeats the reason there are two engines.

## ⭐⭐ THE GATES — ⛔ non-negotiable

7. ⭐ **Whatever writes does not review.** ⭐ **TWO** legs on anything physics-bearing — ⚠ **and the shared
   spec is physics-bearing**, because an error there makes both engines agree on the same wrong thing.
   Orchestrator-written → **Codex + Grok**. Codex-written → **fresh Claude agent + Grok**.
8. ⭐ **LAUNCH LEGS ON SIGHT — ⛔ before I look at the result.** ⚠ A self-check **discharges the felt need**
   for an independent one, ⭐ and it is most convincing when it *finds things*.
9. ⛔⛔ **NO COMMIT BEFORE BOTH LEGS REPORT.** The commit is the **last** step. ⛔ Reviewing the *directive*
   does not pay the tax for the *build*. ⛔ **My own verification is not a leg.**
10. ⭐⭐ **STOP when nothing outstanding changes what is COMPUTED or what may be CLAIMED** — ⛔ **not** when
    both legs are green. ⚠ "A leg that finds nothing is weak evidence" is **my prior**; ⛔ put it in a
    leg's prompt and it becomes a quota.

## ⛔ THE TRAPS — ⚠ all of these were walked into

11. ⛔⛔ **CORRECTNESS IS KING; COST IS NEVER A REASON** to drop a control, narrow a check, or skip a leg.
    ⭐ Scaling work down is the **user's** call.
12. ⛔ **DO NOT BUILD BLINDNESS APPARATUS.** ⚠ The measured failure is **absence of computation**, ⛔ not
    anchoring. ⭐ Quarantine is **CUT**; the three clauses replaced it. ⛔ A do-not-read list is a denylist,
    and a denylist means the architecture is wrong.
13. ⭐ **A finding is not a mandate — ⭐ verify it myself.** ⚠ Legs have been wrong in **both** directions.
14. ⭐ **Ablate to test, ⛔ don't read.** ⭐ A **FORM** control tests physics; a **COEFFICIENT** control
    tests arithmetic. ⛔ Demand a script and its literal stdout — ⚠ a prose re-derivation is the same defect
    relocated into the review.
15. ⭐ **If successive revisions keep breeding defects in the material just changed → ⭐ CHANGE THE AUTHOR.**
    ⛔ Do not fold a fourth time.
16. ⭐ **Prior art is an ORACLE, ⛔ never a PREMISE.** ⭐ Check our computed result against it; ⛔ never
    assume its result for our object. ⚠ Its **conditions** may not be ours.
