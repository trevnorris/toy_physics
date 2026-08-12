# Independent review — the S11 SymPy build directive (round 5 — FINAL, NARROW)

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_build_directive.md`
(uncommitted in the working tree; 55 lines).

⚠⚠ **NARROW REGRESSION CHECK, ⛔ NOT A REVIEW.** ⭐ Four rounds of two legs each preceded it.
⛔ **Do not re-review anything except the four questions below.** Everything else in this file has been
cleared by at least two independent legs, and re-litigating it is the cost this is trying to stop paying.

⭐⭐ **The context: one FACTUAL error was found and fixed.** The population sentence said the export carries
*"every tag emitted for the primary package."* ⛔ That was false about the predecessor's own export, and
⭐ three findings across two rounds were symptoms of it. ⭐ The census is committed at
`research/pde_ledger_v3/directives/_measurements/S11_sympy_build_directive.md` — ⭐ **re-run its commands
rather than trusting the transcription.**

⭐⭐ **Answer only these:**

1. ⭐⭐ **Does the new population sentence match the measurement?** *"the computed objects of the package §7
   identifies as primary and the free symbols they contain."* ⛔ Re-run the census yourself against the
   committed predecessor export and its exporter. ⚠ Does that sentence describe what the predecessor
   actually does — ⛔ over-broad, ⛔ under-broad, or right?
2. ⭐ **Do the three symptoms now have a referent?** (a) the `class` rule's *"a declared symbol carries its
   declared class"*; (b) a primary-package coefficient absent from the import getting a row; (c) `§Q6r`'s
   two-step lookup having a coefficient row to resolve into. ⛔ For each: closed by the population change,
   or still needs its own sentence?
3. ⭐ **Is *"§Q6r's lookup shape, applied to the export this engine writes"* now unambiguous**, and does it
   still ask for a computed object rather than an end-state?
4. ⭐ ***"Apply `F3` to the rows this step writes."*** ⛔ Does pointing rather than transcribing still oblige
   what `F3` exists for, and does bounding it to this step's rows remove the collision with the
   non-regeneration sentence — ⚠ without dropping an obligation `F3` places on a row this step **does**
   write?

⭐ Then one general question, ⛔ and nothing else: did any of the four introduce a **new** divergence, a
value or expectation, or a control the engine would apply to itself?

⚠ ⭐ **If your answer to all five is clean, say so plainly** — ⛔ a leg that manufactures a finding to avoid
returning empty is worse than one that returns empty. ⭐ State what you ran and what would have had to be
true for you to find something.

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
`/tmp/s11dir4_leg_<yourname>/` and report the paths.

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
