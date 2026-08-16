# Independent review: a Wolfram engine repair (script branch)

## Artifact
The WORKING-TREE file
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
— it carries an uncommitted repair on top of commit `19b607ab` (diff it yourself against
`git show 19b607ab:research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
to see exactly what the repair touched).

## What to check
The repair claims to replace unbounded quantifier-elimination calls with one uniform bounded route
(primary QE attempt under a time budget, then an exact-rational fallback decision, then an honest
undecided outcome) at all three QE sites, without changing any emitted object's meaning, weakening
any record, or gating on cell identity. The decision list it must satisfy is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round1_brief.md`
(seven obligations; the spec it cites is `directives/S11_SHARED_PHYSICS.md`, locus protocol §5).
Review the repair's soundness: can it emit a wrong token, a dead payload, or a silently weaker
object; is the route genuinely uniform; can its fallback actually decide anything (able-to-fail
both ways)?

## What you are handed
The artifact, the brief and its `_measurements/` twin, the spec, the committed baseline output
`research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`, and the
build's run records under `/home/trevnorris/.s11_build/fix1_build/` (regression logs, guarded cell
records, emission timelines).

## Required method — this is a SCRIPT artifact
Derive independently. Ablate every load-bearing check and report its literal output; code-reading
alone has repeatedly missed real defects here. Write your own probe scripts BEFORE settling any
contested question, and save both the script and its literal stdout to named absolute paths —
without these, derivation claims are discarded.

⛔⛔ A FORM ABLATION IS MANDATORY. Change the STRUCTURE of a load-bearing object — flip a sign AND
an off-diagonal in the action's stiffness structure, or collapse two independent symbols into one —
re-run an affected cell, and report the literal diff of emitted tags. A coefficient rescale tests
arithmetic; only a form change tests physics. Ask of every emitted claim: WHICH LINE COMPUTED
THIS? Report any payload that survives a form ablation byte-identical, and any `assert`-like
guard that precedes the value it protects.

One-sided corruptions of the repair itself (each: corrupt, re-run one small cell, report the
literal tag diff, restore the copy):
1. Corrupt the fallback decision helper so it returns a wrong-but-well-formed answer — do the
   status tokens and downstream records move? If nothing moves, the fallback is dead code.
2. Corrupt the primary route's returned outcome the same way — same question.
3. Shrink the primary budget to near-zero — the route must visibly fall back (operands must show
   the attempt record), and tokens must remain honest, never flip to a decided value the fallback
   did not compute.
4. Check emission conditionality: no tag's PRESENCE may depend on a payload's value — only on
   package and quantity. Report any conditional emission the repair introduced.
5. Check uniformity in the diff: any route selection on package, dimension, root, stratum, or tag
   identity is a BLOCKING finding.

## Operational constraints (identical for both legs)
- ⛔ Wrap EVERY kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it, move on.
- ⛔ NEVER raise the timeout, never run more than one kernel at a time (two-seat licence).
- ⛔ Do NOT run `XKIN_ANISO 4` or `XKIN_ANISO 2` — both exceed the machine's memory in the strata
  (measured; the guard records are in the build scratch). Use `MAIN 2` (~15 s), `XFORM_*` cells
  (~15–30 s), or `XKIN_ANISO 3` (~255 s) for every run and ablation.
- ⛔ Copy the artifact to /tmp and ablate the COPY. Never modify the working tree.
- ⭐ Save every ablation script AND its literal stdout to named absolute paths and report them.

## Physics filter
Report a finding only if it catches a way the physics record could be wrong; do not report "the
script would be wrong on a different input."

## Report format
Verdict line; BLOCKING findings (each with the corrupted line, the command run, and the literal
output that exposes it); non-blocking findings; what is solid — with the ablation evidence table
(probe → cell run → literal outcome path).
