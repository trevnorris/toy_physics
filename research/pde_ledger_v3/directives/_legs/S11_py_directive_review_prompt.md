# Independent review — the S11 SymPy build directive (round 2)

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_build_directive.md`
(untracked in the working tree; 47 lines).

⚠ **Round 1 of this review, by two independent legs, blocked it on the EXPORT half only.** Its findings, so
you do not re-find them: the export's **row population** was undetermined; the **tag→key transform** was
unspecified in a way that let a builder make `F9` fire on nothing while looking identical; `F6`'s
**publication choice** was unmade; `class` assignment was undetermined; and the engine was no longer
forbidden from writing `S10_exports.py`. ⭐ All five have been edited. ⛔ Round 1 found the **stdout half
clean and buildable**, its citations faithful, no leak, and — on a full 81→38 diff — **nothing dropped**.
⛔ Do not treat any of that as settled by assertion: ⭐ re-derive what you can, ⛔ but do not spend the
review re-litigating the stdout half.

It replaces `S11_sympy_no_ledger_build_directive.md`, which is deleted in the working tree and recoverable
with `git show HEAD:research/pde_ledger_v3/directives/S11_sympy_no_ledger_build_directive.md`. ⚠ That
predecessor was **blocked by two legs** and its scope has since been reversed — ⛔ do not treat it as a
standard the new one must meet, ⭐ but do use it for item 5 below.

## What this artifact is

The build directive Codex will rewrite the S11 SymPy engine from. It is **not** the physics: the physics is
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`, which binds and wins every conflict. The
directive's job is to carry only what the spec does not — scope, the export path, and boundaries.

The engine's products are a flushed stdout tag stream and `scripts/S11_exports.py`, an accumulating ledger
the next step imports.

## Background

- `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — the shared spec. `§4` structural computation,
  `§5` clauses/corollaries/locus protocol, `§6` Q1–Q11, `§7` packages and run record, `§8` tag grammar,
  `§9` premise inventory, `§10` reporting.
- `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — `F1`–`F9` and, after `F9`, an
  **OWED TO THE BUILD REVIEW** checklist. `F9` is the key-write rule.
- `research/pde_ledger_v3/scripts/S9_exports.py`, `S10_exports.py` — the two existing ledgers.
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the previous step's engine and
  its export path (`:1740-2140`). Mechanical precedent only.
- `CLAUDE.md` — binds.

## What to check

1. **Is it buildable without diverging?** Read it as the builder. Name every place two competent builders
   would produce materially different engines from it. ⭐ This is the highest-value item. In particular
   decide, from the directive plus what it points at, exactly **which rows the engine writes** to
   `S11_exports.py` — if you cannot determine that, say so, and say what a builder would most likely do.

2. **Does it point, or restate?** It claims to cite spec obligations rather than re-word them. For each
   citation, open the cited section and judge: does the pointer carry the obligation, is it weaker, or does
   it claim coverage the section does not give? ⚠ In this workstream every re-wording produced a weaker
   rule, and the previous directive's three blocking findings were all this defect.

3. **Does anything it forbids need to exist somewhere?** It bans a completion registry, a directive-local
   exit policy, and encoding the review checklist as engine guards. For each ban, verify the obligation is
   genuinely held elsewhere and name where. ⛔ A ban on machinery whose obligation nothing else carries is
   a hole.

4. **Leakage.** `CLAUDE.md` rule 5: a directive says what to compute, never what anything equals, is
   expected, or was measured. Check nothing hands the builder a value, count, membership, baseline or
   worked example for any S11 object.

5. **What the five edits cost.** Diff `git show HEAD:...S11_sympy_no_ledger_build_directive.md` against
   the current file. ⭐ The question is no longer "what was lost" — round 1 settled that — ⛔ it is whether
   the **five additions** introduced anything: a restatement of an obligation the spec or `F9` already
   carries, a control the engine would apply to itself, a value or expectation, or a new divergence.
   ⚠ In particular the key/imported-key emission (item 2 of the fix) is new machinery in the engine's
   output — judge whether it is genuinely operands, or a claim wearing operands' clothes.

6. **Contradiction.** Against the shared spec, `F1`–`F9`, and `CLAUDE.md`.

## Do not read

- Every other file in `research/pde_ledger_v3/directives/_legs/` — earlier legs and briefs, including the
  brief that produced this rewrite. ⛔ Their framing is not yours to inherit.
- `/tmp/f9*_leg_*/`, `/tmp/s11_fold_leg/`, `/tmp/s11_directive_rewrite.txt` — quarantined workspaces and
  build transcripts.
- `research/pde_ledger_v3/scripts/S11_*` and `research/pde_ledger_v3/_asbuilt/` — the S11 engine is being
  rewritten and the existing file is not a premise.

## Required method

Read the **spec and `F9` first**, form your own view of what a directive must add on top of them, and
**only then** read the directive. Reading it first anchors you to its framing, which is the thing under
test.

⛔ **A prose derivation is worth nothing.** Where a claim is checkable, check it and show literal output —
item 5's diff and item 1's "which rows" question both are. Save anything you run, with its stdout, under
`/tmp/s11dir_leg_<yourname>/` and report the paths.

⛔ Wrap any long-running command in `timeout 600`; a 600s hit is a failed probe, ⛔ never a measurement of
infeasibility. ⛔ Do not modify the working tree. ⛔ Do not run the S11 engine or any engine.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way a builder could satisfy the
directive and produce an engine that does not measure what the spec asks for. No style findings.

⚠ **Brevity is not a defect and is not evidence of one.** The predecessor was longer and was blocked.
⛔ Do not report that something is missing without completing item 5's test on it — most of what is gone
is gone because the spec already carries it.

If nothing survives the filter, say so plainly, and state what you ran and what would have had to be true
for you to find something.

## Quarantine

S11's physics results are not computed and must not be. If a check would require computing S11's spectrum,
determinant, root structure or mode count, **stop and report that instead** — that is itself a finding.
