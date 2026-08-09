# Independent review — the S10 step-record DECISION LIST, before any builder sees it

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_step_record_decisions.md`

**You are one of two independent legs with this identical prompt.** ⛔ Read-only: do not edit the working
tree. Verify with `git status` before you report.

## What this is, and why it is being reviewed

This list will be handed to a builder who rewrites
`research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the step record for S10.

⭐ **The list is the one artifact in this pipeline that is checked zero times while everything downstream
is checked twice.** The builder trusts it. A defect in it propagates into the record and both of the
record's own review legs will be reading an artifact built to a broken instruction.

⚠ **A step record is not code.** Its job is to state what was measured and what may **not** be claimed on
the strength of it. The failure that matters is a record a reader can honestly quote to overstate the
result — not a typo.

## What to check

1. ⛔⛔ **Can any item's acceptance be satisfied with the defect it targets still in place?** ⭐ **This is
   the most important question in this review.** Take each item, imagine the laziest rewrite that passes
   it, and say whether that rewrite still leaves the record wrong. ⚠ Measured on this project: an
   acceptance criterion was written that would have passed with the defect fully intact, because the
   substitution it prescribed collapsed two symbols into one and no discriminating case existed.
2. ⛔ **Does any item state a measured answer, an expected value, a count, a sign, or the reason behind
   one?** The builder must reach every number by reading an artifact. ⚠ An item that says what the answer
   will be is an item the builder can satisfy by transcription.
3. ⛔ **Level-above misses.** For each item: if the builder makes exactly the named thing true, is the
   level above it still broken? ⚠ Measured on this project three times in one session — a named object was
   made true while the property it was supposed to serve stayed false.
4. ⛔ **Is any item already true of the current record?** Open the record and check. An item that cannot
   fail is noise, and it dilutes the ones that can.
5. ⛔ **What is MISSING.** ⭐ **A list is easiest to break by omission, so spend real effort here.** Read
   the current record and the live artifacts, and name anything the rewritten record would need in order
   to be honest that no item requires. ⚠ Include the opposite failure: material the list would cause to be
   **dropped** that a reader needs.
6. ⛔ **Is the scope fence contradicted by any item inside the same list?** ⚠ Measured on this project: a
   list demanded a new emission while forbidding any change to the file that receives it, and the builder
   correctly obeyed the fence.
7. ⛔ **Does any item specify a recipe where it should name an object or a property?** If an item tells the
   builder *how* to derive something rather than what must become true, say so.
8. ⛔ **Does anything in the list commission a check, residual or count whose operands are produced by the
   thing it polices?** ⚠ Seven such were commissioned on this project in one session; most were provably
   unable to fail.
9. **Is the list internally consistent** — does any item require something another item forbids?

## Required method

⛔ **Read the current step record and the live artifacts BEFORE forming a view on the list.** Reading the
list first anchors you to its framing, which is the thing under test.

⭐ **For every claim you make about what an artifact contains, run the command and quote its literal
output.** Save any script you write and its stdout to named absolute paths and report them. ⛔ A prose
assertion about what a file contains will be discarded.

⛔ Do not state anywhere what you expect any count or value to be before you measure it.

## Where things are

- Repo root `/var/projects/toy_physics`, branch `ledger-v3-rebuild`, tree clean at `2e64e5c5`.
- The record under discussion: `research/pde_ledger_v3/steps/S10_two_transverse_photons.md`.
- The front door: `research/pde_ledger_v3/REBUILD_HANDOFF.md`. Governing process: `CLAUDE.md`.
- Committed transcripts: `research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out`,
  `research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`,
  `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out`.
- A comparable record already rebuilt under the current pattern:
  `research/pde_ledger_v3/steps/S9_light_requires_shear.md`.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not run or modify either engine or the comparator.
- ⛔ Never modify the working tree. Work under `/tmp`.
- Wrap any long run in `timeout 600`.

## Physics filter

Report a finding only if it catches a way **the physics, or the ledger's record of it, could be wrong or
could be honestly quoted to overstate the result.** ⛔ Not style, not formatting, not length.

⚠ **"Nothing survives the filter" is a permitted answer and is weak evidence on its own** — if that is your
conclusion, say which items you attacked and how, so the strength of the pass is visible.
