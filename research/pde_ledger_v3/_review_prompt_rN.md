# REVIEW — v3 ledger plan set. Read-only, one pass.

This is an **iterative** review. Prior rounds returned blocking findings which have been folded.
**HEAD has advanced; run `git log --oneline -3` to see it.** Your job: decide whether this is now safe to execute, and find anything the
folding broke or missed.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.

**Read, in order:**
1. `research/pde_ledger_v3/V3_STEP_PLAN.md`   ← the main object
2. `research/pde_ledger_v3/CHARTER.md`
3. `research/pde_ledger_v3/DEFECT_REGISTER.md`
4. `research/pde_ledger_v3/SESSION_REASONING.md`

## Fixes folded since the last round

Round 6: both reviewers STILL BLOCKING, converging on the same first item. Folded:

- ⭐ **The GOVERNING DOCS were stale and contradicted HEAD.** `docs/derivation_walkthrough_plan.md` and
  `research/pde_ledger_v2/walkthrough/DECISIONS.md` both still said the `R2.a_pin` class is OPEN and
  must not be resolved, while commit `407eed94` removed the relation. Both are now current; DECISIONS'
  superseded consequences are marked history, not instruction.
- **S16's remainder was dimensionally inhomogeneous** (force on the left, acceleration in the
  remainder). Corrected, and the length relabelled `a_WT`.
- **S15's force law had the profile normalization silently set to 1** while `I_F` was listed as debt in
  the same plan. `I_F,12` is back in the equation.
- **S1.5 rewritten whole** — no longer leads with a uniqueness claim it retracts, no longer asserts S1's
  retracted three-field inventory; `Cρ` is a recorded choice, not an asserted consequence.
- **S22's table replaced** (not annotated) — partitioned by locus; `J`, `m_defect` and the geon are
  three separate rows with the missing bridges listed as debt.
- **S12 reclassified** as a non-variational source plus separate boundary data.
- **S20a's two classifications made non-optional**; A11 residuals folded.

⭐⭐ **Hunt for the dominant defect class**: a correction *appended* while the original sentence,
heading, table row, equation or summary still asserts the pre-correction claim. Three rounds running,
that has been the main finding — including in documents OUTSIDE `research/pde_ledger_v3/` that the plan
depends on.

## Operating constraints
- **READ ONLY.** Do not modify or commit. **One pass, no clarifying questions** (you will hang).
- You may run Python/SymPy; write scratch outside the repo.
- ⭐ Open the cited file rather than trusting the document's summary of it.
- ⛔ For any claim of the form "there is no X", read the WHOLE artifact first. A prior round caught
  exactly that error (a universal negative drawn from 4% of a 3239-line file).

## Output

```
# Review — <your name/model>

## VERDICT
<CLEAN — SAFE TO EXECUTE / STILL BLOCKING / REGRESSED> + 2-3 sentences
⭐ If CLEAN, say so plainly. Do not manufacture findings to seem rigorous.

## BLOCKING
<numbered; each: claim, file:line, why wrong, what to change. EMPTY IS THE EXPECTED ANSWER
 IF THE PLAN IS SOUND — an empty list is a real verdict, not a failure to look.>

## NON-BLOCKING
<things worth fixing that do not stop execution>

## S0.5 → S1 → S1.5 → S2
These four will be executed side by side with the user next. What concretely goes wrong?

## MATH FLAGS
<table: claim | file:line | your result | agree/disagree>
```

## Standards
A matching number is not evidence. Dimensional agreement is not physical agreement — ask whether both
sides are indexed by the same thing. Falsification is welcome and first-class. Apparatus above physics
has killed two efforts here. Absence of a denial is not evidence.

⛔ Do not be agreeable — but ⛔ do not invent blockers either. If it is ready, say it is ready.
