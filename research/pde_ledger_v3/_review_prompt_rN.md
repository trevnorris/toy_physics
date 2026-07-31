# REVIEW — v3 ledger plan set. Read-only, one pass.

This is an **iterative** review. Prior rounds returned blocking findings which have been folded.
**HEAD is `0e279e23`.** Your job: decide whether this is now safe to execute, and find anything the
folding broke or missed.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.

**Read, in order:**
1. `research/pde_ledger_v3/V3_STEP_PLAN.md`   ← the main object
2. `research/pde_ledger_v3/CHARTER.md`
3. `research/pde_ledger_v3/DEFECT_REGISTER.md`
4. `research/pde_ledger_v3/SESSION_REASONING.md`

## Fixes folded since the last round — verify each is REAL, not cosmetic

- **S1.5 added and numbered** — the substrate action / Madelung balances, previously only a
  parenthetical. Check it actually supplies what S4's "core balance" and S12's momentum/energy
  partners need.
- **S0.5 restated as a real step with its own review legs** — removing `c_γ`/`λγ` from the medium
  contract breaks `C-M1`, `able_to_fail.py` and `test_registry.py`, which are built around `R3`.
  Check the step now says enough to be executable.
- **S16 no longer claims to legitimise the scope**; recast as a conditional response-side
  approximation whose supplied inputs must be enumerated.
- **S15's correction corrected** — S16 does NOT license S15; S15 carries its own Noether /
  control-surface assumptions.
- **S20a separates** the derived ratio `λγ = c_γ/c_s` from the calibrated equality `λγ = 1`.
- **New register rows A11 (worldtube profile scale ≠ mouth radius), A12 (a third healing-length
  convention), A13 (`χ_B` real vs complex)**.

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

## S0.5 → S1.5 → S1 → S2
These four will be executed side by side with the user next. What concretely goes wrong?

## MATH FLAGS
<table: claim | file:line | your result | agree/disagree>
```

## Standards
A matching number is not evidence. Dimensional agreement is not physical agreement — ask whether both
sides are indexed by the same thing. Falsification is welcome and first-class. Apparatus above physics
has killed two efforts here. Absence of a denial is not evidence.

⛔ Do not be agreeable — but ⛔ do not invent blockers either. If it is ready, say it is ready.
