# Independent review — a repaired paper card against the record it summarises

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`

A LaTeX section for a physics ledger, just repaired. You are one of two independent legs; the other is not
visible to you.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not open the card first

For a document, blindness comes from **reading order**. Reading the card first anchors you to its framing,
which is the thing under test.

1. Read `/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step
   record. ⭐⭐ **This is the source of truth.** Closed and reviewed.
2. **Write down, before step 3:** what does this record establish, at what strength, and what does it
   explicitly refuse? What conditions does it attach to its headline result, and where does it attach
   them? ⭐ Keep this list.
3. **Only now** read the card.

## ⛔ Do not read

- Any file under `research/pde_ledger_v3/directives/` matching `S10_paper_card_*decisions*` or
  `*build_prompt*` — the directives this card was built and repaired from. ⭐⭐ **An artifact can satisfy
  its directive and still misrepresent its source, and that is exactly what this leg exists to catch.**
  (They have been moved out of the tree; if you find one, do not read it.)
- Any prior review output.

⚠ You **may** read `directives/S10_paper_card_claim_ledger.md`, `..._dropped.md`, `..._unresolved.md` —
⛔ but only **after step 2 is written down**. They are the builder's self-report and are under review.

## What to check

**Governing question: can a reader quote this card without quoting something S10 did not establish?**

1. ⭐⭐ **Every claim, against the record.** For each assertive sentence: does the record license it, at
   that strength? Report anything that claims more than the record, states a condition more weakly than
   the record states it, or asserts something the record refuses.
2. ⭐⭐ **Misrepresentation by OMISSION.** A card can contain no false sentence and still mislead by
   dropping a qualification the record attaches to a result it quotes. ⚠ Does every result appear with the
   conditions the record attaches to **it**, **at the site where the result is stated** — ⛔ not exiled to
   a later list?
   ⭐⭐ **Test this concretely: take the card's headline sentence alone, as a reader would quote it. Is
   there any case the record measures that satisfies every condition that sentence lists and yet gives a
   different answer?** If so, name it.
3. ⭐ **Attribution.** Where the card explains *why* a measured result comes out as it does, does the
   record's control evidence support that explanation? ⚠ Correct algebra about a supplied model can still
   attribute a measured outcome to the wrong feature of it.
4. ⭐ **Strength of any general statement.** Does the card state at unrestricted generality something the
   record establishes only for the cases it measured?
5. ⭐ **Numbers and loci.** Every number: does the record report it, and **does the cited locus actually
   contain it?** ⛔ A path that resolves is not a source. Open the loci.
6. ⭐ **Dead artifacts.** Does the card cite, quote, or report output from any instrument, file, or
   registry absent from the live tree? Verify by looking.
7. ⭐ **Scope.** Any result belonging to a different step or sector? ⚠ The record's own sentence handing an
   object onward to a later step is legitimate; a claim about what that step will find is not.
8. ⭐ **The self-report.** Does the claim ledger classify honestly? ⛔ An `EXPOSITION` row asserting a
   measured outcome is misclassified. Is any card sentence missing from the ledger?
9. ⚠ **Visibility.** `research/pde_ledger_v3/paper/macros.tex` suppresses one field **by name** in the
   default build. Is any reader-critical content inside it?

## ⭐⭐ 10 — THE DIFF CHECK, and this leg is the only place it happens

This card was **repaired**, not rebuilt. Compare it against its previous committed version:

```
git -C /var/projects/toy_physics diff 2998029f -- research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex
```

⛔⛔ **A repair can evade a defect by DELETING the sentence rather than fixing it.** For every substantive
sentence that the diff **removes**, ask: was that content something the record requires the card to carry?
⚠ A card that silently drops a characterised departure, a limitation, or a measured result is **worse**
than one that stated it imperfectly, and it will not show up in any check of what the card now says.

⭐ Also check the converse: did the repair **introduce** anything new that the record does not license?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way the **card could
misrepresent what was measured**. ⛔ Not LaTeX, not style, not wording preference.

## Method

- ⭐⭐ **Quote both sides for every finding**: the card text, and the record text it fails against. A
  finding without both quotations is not usable.
- ⭐ For anything turning on a number or a locus, **state what you opened and what you found**.
- ⛔ Do not edit the card or anything else. Read-only.
- ⭐ End with: **which items on your step-2 list are absent from the card**, and whether each absence is a
  defect or a legitimate scope choice.

⚠ **A leg that returns "nothing survives the filter" is weak evidence.** If that is genuinely your
conclusion, say so plainly and state what you checked that could have failed and did not.
