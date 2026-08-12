# Independent physics review — F9 round 3, the export-override decision

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`,
**the `F9` section only** — from `## ⭐⭐ F9 — WHEN A STEP'S EXPORT WRITES A KEY THE IMPORTED LEDGER
ALREADY HAS` to the start of `## ⛔ What this list does not decide`. That span includes a **blocked
record** of what round 1 proposed; the blocked record is part of what you are reviewing.

`F1`–`F8` are prior context. Read them to judge contradiction; they are not under review.

This is an **orchestrator-written decision list** — not a build directive. Codex will rewrite the S11
SymPy engine's export path from it. It is the one artifact in the chain the builder trusts without a
second check, so a defect here reaches the engine unfiltered.

## Background you need

Each step's SymPy engine imports the previous step's `LEDGER` (flat dict, `key -> {value, class, step,
…}`), computes its own physics, and re-exports the merged whole. The accumulated ledger is meant to end
up as the model's single list of every quantity defined and every knob.

- `research/pde_ledger_v3/scripts/S9_exports.py`, `S10_exports.py` — the two existing ledgers.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — wrote S10's; export path
  around `:1740-2140`.
- `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree) — the shared spec both S11
  engines read. §5 corollaries, §7 packages, §8 tag grammar.
- `research/pde_ledger_v3/DEFECT_REGISTER.md` entry `C20` — four earlier designs for this problem and the
  measurement that killed each.
- `CLAUDE.md` — binds.

S11's `MAIN` action is S10's `MAIN` plus one additive stiffness term carrying a new coefficient. S11
computes the same *kinds* of object from a *different* action, landing on keys S10 already owns.

## What to check

**Whether `F9` as it now stands can be built from, and whether the subtraction that produced it went too
far or not far enough.** In this order:

1. **The equality decision is the only remaining judgment. Is it now constrained enough?** `F9` declines to
   prescribe a routine and instead fixes four properties of it — total over the row shapes the imported
   ledger carries, computed under the governing premises, elementwise into containers, three-valued.
   Build a pair of objects on which two routines *satisfying all four* still disagree, and say what that
   costs. Then judge whether declining to prescribe the algorithm is defensible under rule 3, or whether it
   leaves a fork that changes what the ledger says.

2. **Does `F9c` break a consumer?** A prefixed key means this step's object is not addressable under the
   bare name. Scan for anything that would now read the predecessor's value believing it to be the current
   one — including the next step's import, `Q6r`'s lookup, and `dimension_key` references. Report each by
   name, or report that you found none and say what you scanned.

3. **Does the merged ledger become self-contradictory?** After `F9c`, the ledger holds both steps' answers
   for the same question under different keys, with nothing marking which reflects the current model.
   Is that a defect of `F9`, or correctly deferred to the step record? Argue it, don't assert it.

4. **The blocked record.** (a) Is the stated reason correct — verify `B_comp > 0` is a declared premise and
   that this does defeat the branch. (b) Round 2's record was found to leak measured S11 payload behaviour;
   it has been rewritten to state the three defects as *properties* with the examples removed. Judge
   whether the properties as now written still hand a builder anything about what S11 should compute —
   and, separately, whether stripping the examples has made any of the three unfalsifiable.

5. **Leakage, generally.** `CLAUDE.md` rule 5: a list says what to compute, never what anything equals, is
   expected, or was measured. Round 1 failed this — its measured-state table handed the builder ledger
   counts and the equality outcome for the whole existing baseline. That table is gone. Check that nothing
   equivalent survives anywhere in `F9`, in any form: count, membership, baseline, line number, or example.

6. **Contradiction.** Against `F1`–`F7`, the shared spec, and the `§5` corollaries. `F9` claims to
   supersede `F2`'s second branch and to amend `F1`; check those are complete and that no third
   contradiction is left standing. Check the PY-only claim.

7. **The guard.** `F9` now retracts its earlier "no guard change" claim and states two requirements: the
   guard sees every write under its **final** key, and it admits an `F9b` whose predecessor is **any**
   imported row. Execute the current guard at `S10_brane_mode_spectrum_sympy_audit.py:2111` against all
   three outcomes and say whether those two requirements are **sufficient** — construct a wrong export that
   satisfies both and still publishes.

## Do not read

- Every other file in `research/pde_ledger_v3/directives/_legs/` — earlier legs on other artifacts,
  including round 1 of this one. Their findings are not yours to inherit.
- `/tmp/f9_leg_grok/`, `/tmp/f9_leg_codex/`, `/tmp/f9r2_leg_grok/`, `/tmp/f9r2_leg_codex/`,
  `/tmp/s11_fold_leg/` — quarantined leg workspaces from earlier rounds.
- Any `research/pde_ledger_v3/scripts/S11_*` file — the S11 engine is being rewritten; the existing file
  is not a premise.

## Required method

Read the **sources of truth first** — the two ledgers, S10's export path, the shared spec — form your own
view of what the chain does today and what S11 needs from it, and **only then** read `F9`. Reading `F9`
first anchors you to its framing, which is the thing under test.

⛔ **A prose derivation is worth nothing.** Where a claim is checkable against the artifacts, check it with
a script and show the literal stdout. Save every script and its output to named absolute paths under
`/tmp/f9r3_leg_<yourname>/` and report those paths. Claims without them will be discarded. Items 1, 2, 4a
and 7 are all decidable by running something. Run it.

⛔ Wrap any long-running command in `timeout 600`. A 600s hit is a failed probe — report it and move on;
⛔ never raise the timeout, and ⛔ never conclude infeasibility from a timeout.
⛔ Do not modify the working tree. Copy anything you need to ablate.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way a builder could satisfy
`F9` while producing an engine whose export does not measure what `F9` claims. No style findings; no
"would be wrong on a different problem".

⭐ **This is the THIRD round.** Round 1 was withdrawn as unbuildable; round 2 was corrected on three
measured defects — an equality test with no branch for an undecided result, a leak of measured payload
behaviour inside a *blocked* record, and a false claim that the publish guard needed no change. ⛔ Do not
treat any of that as a reason to find something here, and ⛔ do not treat the list's brevity or its
amendment history as evidence it is now sound.

⚠ **Round 2's two legs disagreed on two points** — whether the equality routine must be pinned, and whether
the blocked record leaked. ⭐ Form your own view on both; ⛔ do not assume the stricter reading won because
it was acted on.

If nothing survives the filter, say so plainly, and state what you ran and what would have had to be true
for you to find something.

## Quarantine rule

S11's physics results are not yet computed and must not be. If a check would require computing S11's
spectrum, determinant, or root structure, **stop and report that instead** — that is itself a finding
about `F9`.
