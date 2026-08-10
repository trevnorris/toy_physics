# Build — rebuild the S10 paper card against its step record

## Your task

Rebuild `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` so that it satisfies every
item in its decision list.

**Read first, in this order:**

1. `research/pde_ledger_v3/directives/S10_paper_card_decisions.md` — the decision list. It governs. It has
   been through two independent review legs and is closed; ⛔ do not renegotiate it.
2. `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step record. ⭐⭐ **This is the
   authority for every claim and every number.**
3. `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — the card as it stands.

⭐ For house style, `research/pde_ledger_v3/paper/steps/S9_light_requires_shear.tex` is the sibling card in
the same part. ⚠ **Style only** — ⛔ its numbers are from a deleted instrument and are not a model for
anything.

## What you may write

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`
- `research/pde_ledger_v3/directives/S10_paper_card_claim_ledger.md`
- `research/pde_ledger_v3/directives/S10_paper_card_dropped.md`
- `research/pde_ledger_v3/directives/S10_paper_card_unresolved.md`

⛔ **Nothing else in the working tree.** ⛔ Do not edit the step record, the decision list, any engine, any
transcript, or any other card. ⛔ Do not commit.

## Scale

The card is 299 lines and a large fraction of it does not survive the decision list. **A rebuild is
expected; a patch is unlikely to reach the property.** Write the card the record supports.

⭐ Keep the card readable as a paper section. The decision list constrains what may be **claimed**; it does
not ask for a terse document. A reader should finish it understanding what was computed, under what
conditions, and what it does not settle.

## Verification you run yourself

From `research/pde_ledger_v3/paper/`:

```
pdflatex -interaction=nonstopmode pde_ledger_v3.tex
pdflatex -interaction=nonstopmode pde_ledger_v3.tex
```

⭐ It must complete twice with no error. ⚠ Baseline at HEAD is a clean build, 35 pages — so a failure is
yours, ⛔ not pre-existing. ⭐ Report the literal final `Output written on` line.

## Report

State, in your final message:

- what you rebuilt versus what you kept;
- the count of card sentences classified `CLAIM` and `EXPOSITION`;
- **anything in the decision list you could not satisfy, and why** — ⭐ an item you could not meet is a
  finding worth more than a build that pretends it met it. ⛔ Do not invent content to satisfy an item.
