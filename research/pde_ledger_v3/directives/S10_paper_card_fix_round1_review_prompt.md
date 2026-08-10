# Independent review — a fix list, before any builder reads it

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_paper_card_fix_round1_decisions.md`

Orchestrator-written. It will be handed to a builder that repairs a LaTeX paper card. **No builder has
launched.** You are one of two independent legs; the other is not visible to you.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the fix list first

1. `/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step
   record. ⭐⭐ **The authority.** Closed and reviewed.
2. `/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — the
   card as built.
3. **Write down, before step 4:** which claims in that card are not licensed by that record, at that
   strength? Which conditions does the record attach to a result that the card states without them? Which
   citations did you check, and what did you find at the cited locus? ⭐ Keep this list.
4. **Only now** read the fix list.

⭐ You may also read `directives/S10_paper_card_claim_ledger.md`, `..._dropped.md`, `..._unresolved.md`
after step 3 — they are the builder's self-report and are themselves under review.

## What to check

1. ⭐⭐ **Is every item in the fix list real?** For each, quote the card text it names and the record text
   it fails against. ⛔ **An item that is not a genuine defect is a finding** — a fix round that repairs a
   non-defect breeds new defects in the material it touches.
2. ⭐⭐ **Is any item BROADER than the defect it names?** ⚠ A finding arrives with a scope, and the scope is
   part of the finding. An item that generalises a specific defect into a category can downgrade correct
   material elsewhere in the card. Name any item where this is happening and say what would be damaged.
3. ⭐⭐ **Compare your step-3 list against the fix list.** What is on your list and absent from the fix
   list? That is the round's most important failure mode and it is invisible from the list alone.
4. ⭐ **Is any acceptance criterion satisfiable without doing the work it names** — by transcription, by
   deletion, or by asserting the outcome? ⚠ Measured here before: an acceptance item that stated its own
   outcome was passed by copying it.
5. ⭐ **Does any item state an expected value or outcome** the builder could converge on rather than read
   out of the record? ⚠ A prohibition leaks as surely as an assertion.
6. ⭐ **Does any item commission something the record refuses, or restore something it retracted?**
7. ⭐ **Does the scope fence forbid the builder from touching a file an item requires changing?**
8. ⭐ **Is any item unsatisfiable, or satisfiable only by inventing content?**

## Physics filter

Report a finding only if it catches a way the **physics could be wrong** or the **card could misrepresent
the record**. ⛔ Not style, not wording, not LaTeX.

## Method

- ⭐ Quote **both sides** for every finding. A finding without both quotations is not usable.
- ⭐ For anything turning on a citation, **open it and say what you found**.
- ⛔ Do not rewrite the fix list. ⛔ Do not modify the working tree. Read-only.
- ⭐ End with: **which of your step-3 items the fix list handles, and which it misses.**

⚠ **A leg that returns "nothing survives the filter" is weak evidence.** If that is your conclusion, say so
plainly and state what you checked that could have failed and did not.
