# Independent physics review — `F9` round 4, the publish invariant

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`,
**the `F9` section only** — from `## ⭐⭐ F9 — WHEN A STEP'S EXPORT WRITES A KEY THE IMPORTED LEDGER
ALREADY HAS` to the start of `## ⛔ What this list does not decide`.

⚠ It is **uncommitted in the working tree**. Read it there.

`F1`–`F8` are prior context. Read them to judge contradiction; they are not under review.

This is a **decision list**, ⛔ not a build directive. Codex will rewrite the S11 SymPy engine's export
path from it. It is the one artifact in the chain the builder trusts without a second check, so a defect
here reaches the engine unfiltered.

## Background you need

Each step's SymPy engine imports the previous step's `LEDGER` (flat dict, `key -> {value, class, step,
…}`), computes its own physics, and re-exports the merged whole. The accumulated ledger is meant to end up
as the model's single list of every quantity defined and every knob.

- `research/pde_ledger_v3/scripts/S9_exports.py`, `S10_exports.py` — the two existing ledgers.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — wrote S10's; the export path
  and its collision guard are around `:1740-2140`.
- `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree) — the shared spec both S11
  engines read. `§5` corollaries, `§7` packages, `§8` tag grammar.
- `research/pde_ledger_v3/DEFECT_REGISTER.md` entry `C20`.
- `CLAUDE.md` — binds.

S11's `MAIN` action is S10's `MAIN` plus one additive stiffness term carrying a new coefficient. S11
computes the same *kinds* of object from a *different* action, landing on keys S10 already owns.

## What changed, and what to check

An earlier revision ended `F9`'s rule with **two** requirements on the publish guard. Two independent
review legs each built wrong exports that satisfy both and still publish, in four classes: a **payload
that is not this step's derived object**; a **proved-equal pair written under the prefixed key**;
**undecided collapsed to equal** by comparing tokens rather than live objects; and a **check not closed
over everything published**.

Those two requirements have been replaced by a single **invariant of the merged export**. That paragraph
is the primary thing under review.

1. **Does the invariant actually close all four classes?** Take each in turn. For each, construct a
   concrete wrong export and show whether the invariant admits or rejects it. ⛔ Do not reason about it in
   prose — build it and show the literal stdout.

2. **Is there a fifth class it does not close?** The four above were found by review, not by
   construction from first principles. Derive the space of wrong exports yourself and say what else is in
   it. This is the highest-value item in this review.

3. **Is the invariant checkable at all?** It asserts a property of the *published artifact* relative to
   the *imported* ledger. Say what a consumer or a control would have to hold to decide it, and whether
   anything in it is undecidable from the artifact alone. ⚠ If a clause can only be verified by trusting
   the writer's own claim, that is the finding.

4. **Does it point at, or restate?** It cites `S11_SHARED_PHYSICS.md §5` corollary 5 for the payload and
   operand wiring. Read corollary 5. Does the citation carry the obligation faithfully, is it weaker than
   the original, or does it claim coverage corollary 5 does not give? ⚠ In this workstream every attempt
   to re-word a spec obligation produced a weaker rule; the citation form is the fix under test.

5. **Regression.** The rest of `F9` — the three outcomes, the `UNDECIDED` handling, the blocked record,
   the deferred items — was cleared by two legs before this edit. Verify the edit did not weaken any of
   it, and that no contradiction was introduced against `F1`–`F7` or the shared spec.

6. **Leakage.** `CLAUDE.md` rule 5: a list says what to compute, never what anything equals, is expected,
   or was measured. An earlier revision breached this twice. Check that nothing in the current text hands
   a builder a value, count, membership, baseline or worked example — including inside the blocked record.

## Do not read

- Every other file in `research/pde_ledger_v3/directives/_legs/` — earlier legs and briefs, including the
  brief that produced this edit. ⛔ Their findings and their framing are not yours to inherit.
- `/tmp/f9_leg_*`, `/tmp/f9r2_leg_*`, `/tmp/f9r3_leg_*`, `/tmp/s11_fold_leg/`, `/tmp/f9r4_codex_build.txt`
  — quarantined workspaces and build transcripts from earlier rounds.
- Any `research/pde_ledger_v3/scripts/S11_*` file — the S11 engine is being rewritten; the existing file
  is not a premise.

## Required method

Read the **sources of truth first** — the two ledgers, S10's export path and its guard, the shared spec —
form your own view of what the chain does today and what a publish invariant would have to prevent, and
**only then** read `F9`. Reading `F9` first anchors you to its framing, which is the thing under test.

⛔ **A prose derivation is worth nothing.** Where a claim is checkable, check it with a script and show the
literal stdout. Save every script and its output to named absolute paths under `/tmp/f9r4_leg_<yourname>/`
and report those paths. Claims without them will be discarded. Items 1, 2, 3 and 5 are all decidable by
running something against the real ledgers and the real guard. Run it.

⛔ Wrap any long-running command in `timeout 600`. A 600s hit is a failed probe — report it and move on;
⛔ never raise the timeout, and ⛔ never conclude infeasibility from a timeout.
⛔ Do not modify the working tree. Copy anything you need to ablate.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way a builder could satisfy
`F9` while producing an engine whose export does not measure what `F9` claims. No style findings.

⚠ **This is the fourth revision of this section, and the previous three each broke in the material the
one before it changed.** ⛔ Do not treat that history as a reason to find something here, and ⛔ do not
treat the fact that a different author wrote this revision as evidence it is sound. If nothing survives
the filter, say so plainly, and state what you ran and what would have had to be true for you to find
something.

## Quarantine rule

S11's physics results are not yet computed and must not be. If a check would require computing S11's
spectrum, determinant, root structure or mode count, **stop and report that instead** — that is itself a
finding about `F9`.
