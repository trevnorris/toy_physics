# Independent review — a decision list, before any builder reads it

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_paper_card_decisions.md`

It is orchestrator-written. It will be handed to a builder that rebuilds a LaTeX paper card. **No builder
has launched.** You are one of two independent legs; the other is not visible to you.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the decision list first

This is a **document** review, and for a document the blindness comes from **reading order**. Reading the
decision list first anchors you to its framing, which is the thing under test.

1. Read `/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md`.
   This is the **source of truth**: a closed, reviewed measurement record for step S10.
2. Read `/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex`.
   This is the paper card **as it stands at HEAD** — the artifact the builder must rebuild.
3. **Form your own view, in writing, before step 4:** which specific claims in that card are not
   supported by that record, and which are contradicted by it? Write that list down first.
4. **Only now** read the decision list.

## What to check

The question is **whether this decision list, handed to a competent builder, produces a card that says what
the record says and nothing more.** Concretely:

1. **Is any acceptance criterion satisfiable without doing the work it names?** In particular, could a
   builder satisfy it by transcription, by deletion, or by a mechanical text operation, while leaving a
   real defect in place? This exact failure has occurred here before: an acceptance item that stated the
   outcome could be passed by copying it, and a *"every cited path resolves"* item was passed by a live
   file cited at the **wrong lines**.
2. **Does the list state, or strongly imply, a measured value or an expected outcome** that the builder
   could converge on rather than read out of the record? A builder iterating to a green acceptance
   converges on any target it can see. ⚠ A **prohibition** leaks as surely as an assertion.
3. **Compare your step-3 list against what this decision list would actually catch.** Name every item on
   your list that the property, applied honestly by a builder, would **leave in place**. That is the
   list's most important failure mode and it is invisible from the list alone.
4. **The converse: what would it remove that should stay?** The card is a paper artifact with a reader.
   The list draws a distinction between a *claim* (must be licensed by the record) and an *exposition*
   (free). Does that distinction survive contact with the actual card text, or are there passages where it
   is genuinely ambiguous which one they are?
5. **Does any item specify a derivation path where it should name an object?** A list that argues *how* to
   arrive at something manufactures questions; the correct form asks for the object.
6. **Is any item unsatisfiable, or satisfiable only by inventing something?** A requirement whose only
   honest outcome is a fabricated value is worse than no requirement.
7. **Does the list's own scope fence forbid the builder from touching a file that one of its items
   requires changing?**
8. **Does anything the list commissions contradict the record, or restore something the record retracted?**
   ⚠ A superseded artifact is a source of dropped limits, never of settled physics — if the HEAD card
   states an inference the record does not, treat the record as governing and say so.

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way the **card could
misrepresent the record**. Do not report style, formatting, or *"a builder might misread this"* absent a
concrete reading that produces a wrong card.

## Method

- Quote **both sides** for every finding: the decision-list text, and the record or card text it fails
  against. A finding without both quotations is not usable.
- ⛔ Do **not** rewrite the decision list. Report findings; the orchestrator folds them.
- ⛔ Do **not** modify anything in the working tree. This is a read-only review.
- ⭐ State explicitly, at the end, **which of your step-3 items the list handles and which it misses.**

## Scope

Read-only. Do not run engines, do not run `pdflatex`, do not build the card. Everything you need is in the
three documents named above; you may read `research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md` and
`research/pde_ledger_v3/paper/macros.tex` if a finding turns on them.
