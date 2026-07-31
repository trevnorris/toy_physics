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

Round 7: both reviewers converged on the governing method doc again (§8 had been rewritten, §1.1 — the
classification test applied at **every step** — had not). Rather than patch a seventh instance, the
**whole documentary blast radius was enumerated** and closed:

- ⭐ **The structural finding:** step ①'s quarantine rule ("fix what computes, quarantine what only
  narrates") had a gap — **some v2 prose GOVERNS**. The method doc's §1.1 executes at every step, and
  the walkthrough step records are re-banked by v3. Those were quarantined on a false premise.
- **Governing docs brought current:** `docs/derivation_walkthrough_plan.md` §1.1 + the REOPENED
  paragraph · `STATUS.md` · `walkthrough/01_sound_speed.md` (v3 re-banks it as S2) · `stage005` (S2's
  source, which also carried older single-classification wording for `λγ`).
- **Four superseded `_scratch` docs** carry `HISTORY — NOT INSTRUCTION` banners rather than rewrites.
- **The D-01 clause is KEPT** (a relation from unit pins is not a defining equation) — it outlives the
  case that produced it. Only the "open, do not resolve" instruction is removed.
- **S12's closing "only a boundary condition" line replaced** with two separate inventories (source /
  controller functions vs boundary / domain data). **S15's heading** no longer offers S16 as its premise.

⭐⭐ **Hunt for**: (a) any remaining place where a correction was *appended* while the original sentence,
heading, table row, equation or summary still asserts the pre-correction claim — **including outside
`research/pde_ledger_v3/`**; and (b) any *governing* instruction that contradicts HEAD.

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
