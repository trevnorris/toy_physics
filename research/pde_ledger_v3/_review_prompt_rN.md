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

## Fixes folded since the last round — verify each is REAL and ENFORCEABLE, not a warning

Round 4: one reviewer CLEAN (empty blockers), the other found three, all of the same shape — a finding
recorded in one place while the operative step entries said the opposite. All three are folded:

- **A13 is now a GATE attached to S1, S5 and S12.** S1 previously said "Register: none", so under the
  plan's own done-rule it could bank without confronting the real-vs-complex `χ_B` contradiction.
- **S8/S9 committed to one layout** — `{ρ_br, μ_R}` and the R10 debt both start at **S8**; S9 is pure
  consequence. The earlier "pick one" left the ambiguity live.
- **S14/S15/S20 are each marked `CONDITIONAL ON S14a`**, the charter's carry-forward list now lists the
  `1/r²` law, the attractive sign and the stage009 falsifier as conditional rather than carried, and the
  order is `S13 → S14a → S14`.

⭐ **Check specifically for the "warning vs rule" failure**: a finding stated in the register or in a
note, while the step that must act on it says something else.

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
