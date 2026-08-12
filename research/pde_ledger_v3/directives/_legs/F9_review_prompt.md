# Independent physics review — the F9 export-override decision

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`,
**the `F9` section only** (everything from `## ⭐⭐ F9 — WHEN A STEP'S EXPORT OVERRIDES AN IMPORTED KEY`
to the end of that section). `F1`–`F8` are prior context: read them to judge whether `F9` contradicts
them, but they are not what you are reviewing.

This is an **orchestrator-written decision list**. It is not a build directive. It will be handed to a
builder (Codex) who rewrites the S11 SymPy engine's export path from it. It is the one artifact in the
chain the builder trusts without a second check, so a defect here reaches the engine unfiltered.

## Background you need

The v3 ledger chain: each step's SymPy engine imports the previous step's `LEDGER` (a flat dict,
`key -> {value, class, step, ...}`), computes its own physics, and re-exports the merged whole. The
accumulated ledger is meant to end up as the model's single list of every quantity defined and every knob.

- `research/pde_ledger_v3/scripts/S9_exports.py`, `S10_exports.py` — the two existing ledgers.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the engine that wrote S10's,
  including the export path around `:1740-2140`.
- `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree) — the shared spec both S11
  engines read. §7 defines the packages; §8 the tag grammar; §5 the corollaries.
- `research/pde_ledger_v3/DEFECT_REGISTER.md`, entry `C20` — four earlier designs for this same problem
  and the measurement that killed each.
- `CLAUDE.md` — binds.

S11's `MAIN` action is S10's `MAIN` action plus one additive stiffness term carrying a new coefficient.
S11 therefore computes the same *kinds* of object (Lagrangian, Euler–Lagrange system, determinant,
spectrum) from a *different* action, and lands on keys S10 already owns.

## What to check

**Whether `F9` can be built from, and whether building it would measure the physics it claims to.**

Specifically, and in this order:

1. **Is the three-way outcome actually COMPUTABLE, or does it hide a judgment call?** `F9a`/`F9b`/`F9c`
   claim to be decided by computation with nothing declared. Construct the decision procedure yourself and
   say where it becomes a choice rather than a computation. If two competent builders could route the same
   key to different outcomes, that is the finding.

2. **`F9b`'s specialising substitution.** `F9` says it must be *derived* from the two actions, with the
   imported action read out of the imported `LEDGER`. Verify that the imported action is actually there and
   is in a form from which a substitution can be derived — do not assume it. Then: is the substitution
   unique? What happens when several substitutions carry one action to the other, or none does, or one does
   only after simplification?

3. **Does `F9b` lose something the chain needs?** A supersession overwrites the imported value. S10's result
   is a correct statement about a real sub-model. Say what is lost from the merged ledger and whether any
   consumer needs it.

4. **The guard.** Item 5 requires the guard to admit a nonzero value residual exactly when the reduction
   residual is zero. Is that condition sufficient to keep the guard able to fail? Construct a wrong export
   that passes it.

5. **Leakage.** `CLAUDE.md` rule 5: a spec or list says what to compute, never what anything equals, is
   expected, or was measured. A builder iterating to exit 0 converges on any target it can see. Does `F9` —
   including its measured-state table — hand the builder a value, a count, a membership, or an expected
   outcome for any key? The table is there to justify the decision to a human; judge whether it also tells
   the builder what to produce.

6. **Contradiction.** Does `F9` contradict `F1`–`F7`, the shared spec, or the corollaries in `§5`? Does its
   PY-only claim hold — is there any joint property here that the WL engine would need and cannot have?

7. **The measured table.** Verify every number in it against the artifacts. Report any that is wrong.

## Do not read

- `research/pde_ledger_v3/directives/_legs/` — every other file in it. These are earlier review legs on
  other artifacts; their findings are not yours to inherit.
- `/tmp/s11_fold_leg/` — quarantined. It contains computed S11 physics.
- Any file under `research/pde_ledger_v3/scripts/` named `S11_*` — the S11 engine is being rewritten and
  the existing file is not a premise.
- `docs/method_prior_art_brief.md`.

## Required method

This is a **decision list**, so read the **sources of truth first** — the two ledgers, S10's export path,
and the shared spec — form your own view of what the chain does today and what S11 needs from it, and
**only then** read `F9`. Do not read them in the other order: reading `F9` first anchors you to its framing,
which is the thing under test.

⛔ **A prose derivation is worth nothing here.** Where a claim is checkable against the artifacts, check it
with a script and show the literal stdout. Save every script and its output to named absolute paths under
`/tmp/f9_leg_<yourname>/` and report those paths. Claims without them will be discarded.

In particular, items 1, 2, 4 and 7 are all decidable by running something against `S9_exports.py`,
`S10_exports.py` and the S10 engine. Run it.

⛔ Wrap any long-running command in `timeout 600`. A 600s hit is a failed probe — report it and move on;
⛔ never raise the timeout, and ⛔ never conclude infeasibility from a timeout.
⛔ Do not modify the working tree. Copy anything you need to ablate.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way a builder could satisfy
`F9` while producing an engine whose export does not measure what `F9` claims. Do not report style, and do
not report "the list would be wrong on a different problem".

If nothing survives that filter, say so plainly — but a leg that finds nothing is weak evidence, so state
what you ran and what would have had to be true for you to find something.

## Quarantine rule

The S11 physics results are not yet computed and must not be. If your check would require computing S11's
spectrum, determinant, or root structure, **stop and report that instead** — that is a finding about `F9`
(it would mean `F9` cannot be verified without the answers it is supposed to gate).
