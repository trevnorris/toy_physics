# Independent review — a paper card against the record it summarises

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`

A rebuilt LaTeX section for a physics ledger. You are one of two independent legs; the other is not
visible to you.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not open the card first

For a document, blindness comes from **reading order**, and it is as load-bearing here as quarantine is
for a script. Reading the card first anchors you to its framing, which is the thing under test.

1. Read `/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step
   record. ⭐⭐ **This is the source of truth.** It is closed and reviewed; every number and cited locus in
   it was verified by three independent passes.
2. **Write down, before step 3:** what does this record establish, at what strength, and what does it
   explicitly refuse? What conditions does it attach to its headline result? Keep this list.
3. **Only now** read the card.

## ⛔ Do not read

- `research/pde_ledger_v3/directives/S10_paper_card_decisions.md` and
  `research/pde_ledger_v3/directives/S10_paper_card_build_prompt.md` — the directives this card was built
  from. ⭐⭐ **An artifact can satisfy its directive and still misrepresent its source, and that case is
  exactly what this leg exists to catch.**
- any git history of the card, and any prior review output.

⚠ You **may** read `research/pde_ledger_v3/directives/S10_paper_card_claim_ledger.md`,
`..._dropped.md` and `..._unresolved.md` — but ⛔ **only after step 2 is written down.** They are the
builder's own self-report and are themselves under review.

## What to check

**The governing question: can a reader quote this card without quoting something S10 did not establish?**

1. ⭐⭐ **Every claim, against the record.** For each assertive sentence in the card, does the record
   license it, at that strength? Report every sentence that claims more than the record does, states a
   condition more weakly than the record states it, or asserts something the record refuses.
2. ⭐⭐ **Misrepresentation by OMISSION.** A card can contain no false sentence and still mislead by
   dropping a qualification, a condition, or a limitation the record attaches to a result it quotes. ⚠ In
   particular: does every result appear with the conditions the record attaches to **it**, at the site
   where the result is stated — ⛔ not exiled to a later list?
3. ⭐ **Attribution.** Where the card explains *why* a measured result comes out as it does, does the
   record's own control evidence support that explanation? ⚠ A passage can be correct algebra about the
   supplied model and still attribute a measured outcome to the wrong feature of it.
4. ⭐ **Numbers and loci.** Every number in the card: does the record report it, and does the cited source
   actually contain it? ⛔ **A path that resolves is not a source.** Check the cited locus, not just that
   the file exists.
5. ⭐ **Dead artifacts.** Does the card cite, quote, or report output from any instrument, file, or
   registry that does not exist in the live tree? Verify by looking.
6. ⭐ **Scope.** Does the card contain a result belonging to a different step or a different physics
   sector? ⚠ Note the distinction: the record's own sentence handing an object onward to a later step is
   the record's boundary and is legitimate; a claim about what that later step will find is not.
7. ⭐ **The self-report.** Once you have your own list: does the claim ledger classify sentences honestly?
   ⛔ An `EXPOSITION` row that asserts a measured outcome is misclassified. Is any card sentence missing
   from the ledger entirely?
8. ⚠ **Visibility.** `research/pde_ledger_v3/paper/macros.tex` suppresses one field **by name** in the
   default build. Is any reader-critical content sitting inside it?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way the **card could
misrepresent what was measured**. ⛔ Do not report LaTeX style, wording preference, or *"a reader might
misunderstand"* absent a concrete misreading the text supports.

## Method

- ⭐⭐ **Quote both sides for every finding**: the card text, and the record text it fails against. A
  finding without both quotations is not usable and will be discarded.
- ⭐ For any finding about a number or a locus, **state what you looked at and what you found there.**
- ⛔ Do **not** edit the card or anything else in the working tree. This is read-only.
- ⛔ Do **not** rewrite the card or propose replacement prose beyond what is needed to state the finding.
- ⭐ End with: **which items on your step-2 list are absent from the card**, and whether each absence is
  a defect or a legitimate scope choice.

⚠ **A leg that returns "nothing survives the filter" is weak evidence.** If that is genuinely your
conclusion, say so plainly and state what you checked that could have failed and did not.
