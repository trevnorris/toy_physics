# REVIEW — v3 scope widened to charge and magnetism. Read-only, one pass.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.
**Run `git log --oneline -3` first.** ⛔ Do not trust a hash written in any document, including this one.

The changeset under review is the most recent commit. `git show --stat HEAD` and `git diff HEAD~1` will
show you exactly what moved.

## What changed, factually

1. **`CHARTER.md`** — v3's scope was `gravity · light · gravitomagnetism`. It is now
   `gravity · light · gravitomagnetism · charge · magnetism`. The EM/charge cluster had been listed under
   "Not carried — deferred"; that deferral was retired.
2. **`V3_STEP_PLAN.md`** — a new `PHASE 4b — charge and magnetism` with seven steps `Q1`–`Q7`, placed
   after PHASE 4 (gravity) and before PHASE 5 (the knit). `S22`, the deliverable, was re-scoped from
   gravity's interior debts to a listing organized around the throat solve, sectioned per sector.
   §1's step-shape counts moved `6/8/4` → `9/9/7`.
3. **`docs/derivation_walkthrough_plan.md`** — §3's phase-6 rule. The table cell asserted "should
   introduce **nothing** new"; a paragraph below it retracted that; the replacement forbade all new
   inputs at integration. Both the cell and the replacement were rewritten.
4. **`DEFECT_REGISTER.md`** — a new prose entry `A-CAND` under section A.

## What to check

⭐ These are the questions, not a list of known answers. **An empty BLOCKING list is a real verdict.**

1. **Is the scope widening supported by the loci it cites?** The charter's argument rests chiefly on
   `docs/model_map.md:178` and `research/pde_ledger_v2/notes/parameter_register.md:143`, `:228`, `:329`.
   ⛔ Open them. Do the cited lines carry the claim made of them, or has a sentence been stretched?
2. **Are `Q1`–`Q7` faithful to their sources?** Every class assignment (`EARNED` / `DERIVED` /
   `POSTULATED` / `R1_REQUIRED` / `FREE-UNREDUCED`) and every quotation carries a locus. Open them.
   ⭐ A promotion — something recorded as conditional or postulated in the source but reading as earned
   in the step — is the specific failure that matters most here.
3. **Does the rewritten phase-6 rule still have teeth?** It now distinguishes new *inputs* (forbidden if
   they revise an earlier sector) from new *consequences* (expected). Is that distinction operable, or
   can any new input be relabelled a consequence to evade it?
4. **`A-CAND`** claims to record a candidate pattern *without* making the identification. Read it and
   decide whether it in fact makes the identification while disclaiming it.
5. ⭐⭐ **Append residue.** Hunt for any place — **including outside `research/pde_ledger_v3/`** — where
   this changeset appended a correction while an original sentence, heading, table row, equation or
   summary still asserts the pre-correction claim. Two of the four files carry retraction paragraphs
   sitting near text they retract.
6. **Do the loci resolve?** Cited line ranges in this repo are not content-verified by any tool. Several
   files in this changeset changed length. ⭐ Spot-check citations into the four changed files and report
   any that no longer land on what they claim.
7. **Placement.** `PHASE 4b` runs after gravity. `docs/derivation_walkthrough_plan.md` §3's phase table
   orders defects/charge *before* flow/gravity. Is the divergence adequately recorded, and is the stated
   dependency (PHASE 4b needs PHASE 3, not the PN ladder) correct?

## Operating constraints

- **READ ONLY.** ⛔ Do not modify, stage, or commit anything. ⛔ One pass, no clarifying questions.
- You may run Python/SymPy. Write scratch **outside** the repository.
- ⭐ Open the cited file rather than trusting a document's summary of it.
- ⛔ For any claim of the form "there is no X", read the WHOLE artifact first. `wc -l` before any
  universal negative. A prior round produced exactly that error — a universal negative drawn from 4% of
  a 3239-line file.

## Output

```
# Review — <your name/model>

## VERDICT
<CLEAN — SAFE TO EXECUTE / STILL BLOCKING / REGRESSED> + 2-3 sentences
⭐ If CLEAN, say so plainly. ⛔ Do not manufacture findings to seem rigorous.

## BLOCKING
<numbered; each: claim, file:line, why wrong, what to change.
 EMPTY IS A REAL VERDICT, not a failure to look.>

## NON-BLOCKING

## LOCI THAT DO NOT RESOLVE
<table: citation | what it claims | what the lines actually say>

## MATH FLAGS
<table: claim | file:line | your result | agree/disagree>
```

## Standards

A matching number is not evidence. Dimensional agreement is not physical agreement — ask whether both
sides are indexed by the same thing. Falsification is welcome and first-class. Apparatus above physics
has killed two efforts on this project. Absence of a denial is not evidence.

⛔ Do not be agreeable — ⛔ but do not invent blockers either. If it is ready, say it is ready.
