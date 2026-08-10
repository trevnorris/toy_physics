# Build — fix round 1 on the S10 paper card

## Your task

Apply every item in `research/pde_ledger_v3/directives/S10_paper_card_fix_round1_decisions.md` to
`research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`.

**Read first, in this order:**

1. `research/pde_ledger_v3/directives/S10_paper_card_fix_round1_decisions.md` — the fix list. It governs.
   ⭐ It has been through two independent review legs and is folded and closed; ⛔ do not renegotiate it.
   ⚠ Note its **anti-broadening fences** and the item-5 half that was **cut** — ⛔ do not reinstate that.
2. `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step record. ⭐⭐ **The authority for
   every claim and every number.** ⛔ Do not re-derive; ⛔ do not open an engine to check the record.
3. `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — the card to repair.
4. `research/pde_ledger_v3/directives/S10_paper_card_decisions.md` — the governing property the card was
   built to, still in force.

⚠ For the **evidence loci** the fix list asks you to cite, the record already cites them for the same
objects. ⭐ Use what the record cites. ⭐ Where an item asks you to check a locus, **open it**.

## This is a REPAIR, not a rebuild

⭐ The card is largely correct: two independent legs verified **every number** against the raw transcripts
and the comparator and found none wrong. ⛔ **Do not rewrite it.** Change what the fix list names and what
those changes require for coherence, and ⛔ nothing else.

⛔ **Do not fix anything not on the list.** If you find something else, **report it** — it goes to the
orchestrator, not into this diff.

## What you may write

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`
- `research/pde_ledger_v3/directives/S10_paper_card_claim_ledger.md`
- `research/pde_ledger_v3/directives/S10_paper_card_dropped.md`
- `research/pde_ledger_v3/directives/S10_paper_card_unresolved.md`

⛔ Nothing else. ⛔ Do not edit the step record, the shared spec, `V3_STEP_PLAN.md`, any engine, any other
card, or any prior-art document. ⛔ Do not commit.

## Verification you run yourself

From `research/pde_ledger_v3/paper/`:

```
pdflatex -interaction=nonstopmode pde_ledger_v3.tex
pdflatex -interaction=nonstopmode pde_ledger_v3.tex
```

⭐ Twice, no error. ⚠ The pre-fix card builds clean at 35 pages, so a failure is yours. ⭐ Report the literal
final `Output written on` line.

## Report

- **Item 1: name the package you checked the headline's condition list against, and what separated it.**
  ⛔ A report that only asserts the list is now complete does not discharge it.
- For every citation you touched: **what you opened and what you found there.**
- The ledger's `CLAIM`/`EXPOSITION` tallies, and **every sentence you deleted, with why.**
- ⭐ Anything on the list you could **not** satisfy, and why. ⛔ Do not invent content to satisfy an item —
  an item you could not meet is a finding worth more than a build that pretends it met it.
