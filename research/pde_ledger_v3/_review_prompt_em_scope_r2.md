# REVIEW — round 2. The round-1 findings are folded. Read-only, one pass.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.
**Run `git log --oneline -3` first.** ⛔ Do not trust a hash written in any document, including this one.

Round 1 (two independent reviewers) returned STILL BLOCKING on the scope-widening commit. Those findings
are folded in the most recent commit. `git show --stat HEAD` and `git diff HEAD~1` show exactly what moved;
`git show -s --format=%B HEAD` states what each change was for.

⭐ **Your job: decide whether this is now safe to execute, and find what the folding broke or missed.**
A half-fold — fixing one instance of a defect while leaving its siblings — is the specific failure this
project keeps having, and round 1 found four of them.

## What was folded

1. **Scope residue.** The widening had been appended to `CHARTER.md` §1 while operative scope statements
   elsewhere still declared the old scope — including acceptance criterion (4), the definition of done,
   and `STATUS.md`, the repository front door. Four sites were corrected.
2. **Q2** asserted the deflection "is HELD by the puncture." It now earns only the source identity and
   far-field form, carries the core holder as an R1 debt, records three candidate holder mechanisms
   without selecting one, and defers `C4`/`C1` explicitly.
3. **Q4** split into four recorded classes rather than one "interior quantities" bucket.
4. **Q6** split into four objects across two routes; edge `R71` added.
5. **Drain dependency** changed from a phase-wide claim to per-step dependencies.
6. **S21** rewritten against the corrected phase-6 rule (it carried the superseded one).
7. **A-CAND** recast from a stated shared shape to a question to be tested.
8. **Two stale counts** ("ten" → "thirteen"), one inside §0.

## What to check

⭐ These are questions, ⛔ not a list of known answers. **An empty BLOCKING list is a real verdict.**

1. **Is each fold actually correct**, or does it substitute a new error? ⛔ Open the cited sources; do not
   grade the change against the commit message's description of it.
2. ⭐⭐ **Are there SIBLINGS still unfixed?** For every defect above, search for other instances of the
   same defect elsewhere — **including outside `research/pde_ledger_v3/`**. Round 1's largest finding was
   that a correction had been applied in one place and not its four siblings.
3. **Q2's three candidate mechanisms.** Is recording competing mechanisms without selecting one an honest
   deferral, or does it let an unresolved question read as characterised? Are the three fairly stated
   against their sources?
4. **Q4 and Q6's class assignments** — verify each against `parameter_register.md` and the stage notes.
   ⭐ A promotion (something conditional or postulated in the source reading as earned in the step) is the
   failure that matters most.
5. **Per-step drain dependencies** — are they right for each of Q1–Q7?
6. **S21 against `docs/derivation_walkthrough_plan.md`'s corrected rule** — consistent now, and does the
   rule still have teeth as S21 applies it?
7. **A-CAND** — does the question form still assert what it disclaims?
8. ⛔ **LOCI.** Cited line ranges are not content-verified by any tool here, and four separate citations
   went stale during this work — one of them landing on a **partial** match that a careless check would
   pass. Several files changed length again. ⭐ Check citations into every changed file, and note that
   citation FORMAT varies: `FILE.md:120`, `` `:120-125` `` after a filename, and prose references all occur.
9. **Anything round 1 missed.** Two reviewers read this and two further defects were found afterwards by
   verifying their reports. Assume more remain.

## Operating constraints

- **READ ONLY.** ⛔ Do not modify, stage, or commit. ⛔ One pass, no clarifying questions.
- You may run Python/SymPy. Write scratch **outside** the repository.
- ⭐ Open the cited file rather than trusting a document's summary of it.
- ⛔ For any claim of the form "there is no X", read the WHOLE artifact first. `wc -l` before any universal
  negative.

## Output

```
# Review — <your name/model>

## VERDICT
<CLEAN — SAFE TO EXECUTE / STILL BLOCKING / REGRESSED> + 2-3 sentences
⭐ If CLEAN, say so plainly. ⛔ Do not manufacture findings to seem rigorous.

## BLOCKING
<numbered; each: claim, file:line, why wrong, what to change.
 EMPTY IS A REAL VERDICT, not a failure to look.>

## SIBLINGS STILL UNFIXED
<instances of a folded defect that remain elsewhere; empty is a real answer>

## NON-BLOCKING

## LOCI THAT DO NOT RESOLVE
<table: citation | what it claims | what the lines actually say>

## MATH FLAGS
<table: claim | file:line | your result | agree/disagree>
```

## Standards

A matching number is not evidence. Dimensional agreement is not physical agreement — ask whether both
sides are indexed by the same thing. Falsification is welcome and first-class. Apparatus above physics has
killed two efforts on this project. Absence of a denial is not evidence.

⛔ Do not be agreeable — ⛔ but do not invent blockers either. If it is ready, say it is ready.
