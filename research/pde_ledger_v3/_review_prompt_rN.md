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

Round 5: one reviewer CLEAN (third consecutive, empty blockers); the other found four, all of one
class — **an amendment appended below an original claim, leaving both voices readable**, so a step
could bank the retracted version. All four are now fixed **by rewriting the original in place**:

- **S1's field count is 2, corrected in S1 itself** (`ψ` + the A13-selected order field). `U(ρ)` is an
  EOS energy density, not a field.
- **A11 is a gate on S16 and S22**, with distinct symbols `a_WT` (worldtube profile width) vs
  `a_mouth`, and the missing bridge between them listed as debt.
- **S22 is gated on A4, C4, C9**, partitioned by locus, and records `J` / `m_defect` / the geon as
  three distinct unbridged objects.
- **The charter's scope section is rewritten**, not amended: the justification is stated as conditional
  and response-side, "will not touch" → "will not claim to discharge", and a clean close no longer
  claims far-field independence.

⭐⭐ **Hunt specifically for this class**: any place where a correction was *appended* while the original
sentence, heading, table row or summary still asserts the pre-correction claim. That has been the
dominant defect for three rounds.

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
