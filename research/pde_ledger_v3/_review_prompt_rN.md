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

## Fixes folded since the last round — verify each is REAL, not cosmetic

Round 3 returned a split verdict (one CLEAN, one REGRESSED with five blockers). All five are folded:

- **S0.5** now specifies an executable lifecycle (retire the two quantities, remove `R3`) and
  enumerates the touchpoints, including `registry_read.py`'s propagation smoke test and `README.md`'s
  canonical example, both hardwired to `lambda_gamma`.
- **S1.5 order corrected to S1 → S1.5** (it consumes S1's primitives), plus: `U(ρ)` is **not** uniquely
  fixed by the pressure identity (`U = Kρ⁵/4 + Cρ`; the reference setting `C = 0` must be stated), and
  the S1 field count is flagged for correction.
- **S8/S9 provenance order** — `{ρ_br, μ_R}` must enter at S8, not after it.
- **A13** upgraded from listed to must-resolve-before-S1/S5/S12.
- ⭐ **S14a added and BLOCKING** — correcting S12 to the dynamical `Γ_B` law severed the chain to the
  imported gravity results, so S14/S15/S20 are now marked CONDITIONAL until a bridge derives the
  projected order-loss source, `J`, the `J → Q` map and the `ω→0` return law.

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
