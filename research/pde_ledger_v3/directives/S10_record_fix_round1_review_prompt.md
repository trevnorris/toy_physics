# Independent review — S10 record FIX-ROUND decision list, before any builder sees it

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_record_fix_round1_decisions.md`

**You are one of two independent legs with this identical prompt.** ⛔ Read-only: do not edit the working
tree. Verify with `git status` before you report.

## What this is

`research/pde_ledger_v3/steps/S10_two_transverse_photons.md` is the permanent record for S10 — the step
establishing light's transverse mode count from the brane's own dimension. It was rewritten, then read by
two independent review legs.

⭐ **Those legs found no false number and no unsupported locus in it.** Every defect they found was an
**omission**: something true of the committed artifacts that the record does not say, or says too softly.

This list is the fix round. It will be handed to a builder. ⭐ **It is checked zero times while everything
downstream is checked twice** — a defect in it propagates into the permanent record.

⚠ **A step record's job is to state what was measured and what may NOT be claimed.** The failure that
matters is a record a reader can honestly quote to overstate the result.

## What to check

1. ⛔⛔ **Take each item, imagine the laziest edit that satisfies it, and say whether the record is still
   wrong.** ⭐ **This question produced both blockers in the previous round of this review** — spend your
   effort here.
2. ⛔ **Does any item state a measured answer, value, count, sign or the reason behind one?** The builder
   must reach every number by measuring. ⚠ An item that supplies its own answer can be satisfied by
   transcription.
3. ⛔⛔ **A RESTORATION ROUND CAN RESTORE A DEFECT.** This list reinstates limits dropped from an earlier
   version of the record whose instrument has since been deleted. For each: is the thing being restored
   still true of the **committed artifacts**, or only of the deleted instrument? ⭐ Check the ones you can.
4. ⛔ **Item 2 draws a line** — one half of an old passage is to be restored and the other half
   permanently excluded. **Is that line in the right place?** Is the excluded half really an over-claim,
   and is the retained half really measured?
5. ⛔ **Is any item already true of the current record?** Open it and check. An item that cannot fail is
   noise and dilutes the ones that can.
6. ⛔ **What is MISSING** — including the opposite failure: anything this round would cause to be dropped
   or degraded. ⚠ The list forbids reworking material the previous legs confirmed correct; check that
   fence is not itself blocking a needed change.
7. ⛔ **Level-above:** if the builder makes exactly each named thing true, is the property it serves still
   false?
8. ⛔ **Recipe vs object** — does any item prescribe *how* to derive something instead of naming what must
   become true?
9. ⛔ **Does anything commission a check, residual or count whose operands are produced by the thing it
   polices?**
10. **Internal consistency** — does any item require what another forbids? Is the scope fence contradicted
    by an item inside the same list?

## Required method

⛔ **Read the current record and the live artifacts BEFORE forming a view on the list.** Reading the list
first anchors you to its framing, which is the thing under test.

⭐ **For every claim you make about what an artifact contains, run the command and quote its literal
output.** Save any script and its stdout to named absolute paths and report them. ⛔ A prose assertion
about what a file contains will be discarded.

⛔ Do not state anywhere what you expect any count or value to be before you measure it.

## Where things are

- Repo root `/var/projects/toy_physics`, branch `ledger-v3-rebuild`. ⚠ The working tree is **dirty by
  design** — the record under discussion and the front door are both modified and uncommitted.
- The record: `research/pde_ledger_v3/steps/S10_two_transverse_photons.md`.
  Its previous version: `git show HEAD:research/pde_ledger_v3/steps/S10_two_transverse_photons.md`.
- Committed transcripts: `research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out`,
  `research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`,
  `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out`.
- Also live: `research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`,
  `research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md`,
  `research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md`,
  `research/pde_ledger_v3/scripts/S9_exports.py`, `research/pde_ledger_v3/scripts/S10_exports.py`.
- Governing process: `CLAUDE.md`. Front door: `research/pde_ledger_v3/REBUILD_HANDOFF.md`.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not run or modify either engine or the comparator in the tree.
- ⛔ Never modify the working tree. Work under `/tmp`.
- ⭐ A long run that is printing is fine; one that has gone **silent** is a failed measurement — kill it
  and report, ⛔ do not raise its budget.

## Physics filter

Report a finding only if it catches a way **the physics, or the ledger's record of it, could be wrong or
could be honestly quoted to overstate the result.** ⛔ Not style, not formatting, not length.

⚠ **"Nothing survives the filter" is a permitted answer and is weak evidence on its own** — if that is
your conclusion, say which items you attacked and how.
