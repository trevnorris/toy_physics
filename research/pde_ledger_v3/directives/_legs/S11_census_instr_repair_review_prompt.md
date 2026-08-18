# Independent review — S11 census-instrument REPAIR diff and its round-2 census

## Artifact

The repair diff `git diff 90ab5e2d..89ed80c9 -- research/pde_ledger_v3/reduction/` (five modules;
working tree matches `89ed80c9`), the five repaired modules themselves under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`, and the round-2 census outputs they
produced under `/home/trevnorris/.s11_build/census_build2/` (sha256s in `artifact_sha256.stdout`
there; build report `build_report.md` there).

## Scope — this is a SCOPED review of the diff, not a re-review of the instruments

A full two-leg review of the round-1 instruments (commit `90ab5e2d`) already ran; its verified
defect list became the repair directive
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_repair_directive.md`
(commit `33babf8d`) — read it; its seven numbered classes are the acceptance surface for this
diff, together with the folded brief
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_locus_census_instrument_brief.md`
(commit `fa8c58b3`) whose eight obligations still govern. Do not re-litigate logic outside the
diff that the full legs cleared. Your questions:

1. **Does each of the seven repairs actually implement the directive's required semantics** — in
   the code path the census executes, not in a comment or an unreachable branch? For each class,
   find the diff hunk, then demonstrate the behavior with your own script against real record
   payloads.
2. **Did any repair weaken verdict logic the directive says must stay?** The directive's closing
   constraint names them: reducer arithmetic, population reconciliation, the NONZERO refutation
   direction, the PROVED_CONSISTENT decisions, the two-route residual comparisons. A repair that
   silences a false-verdict class by also silencing the true-verdict class is the worst outcome
   this review can miss.
3. **Is the round-2 calibration real?** Verify the planted cases traverse the production
   entrypoint/parser/solver (code path, not comment), are byte-shaped like the committed records,
   and that each planted case was demonstrated ABLE TO FAIL (the preserved
   `calibration.failed1.stdout` and `calibration.pass_before_*.stdout` files in
   `census_build2/` are the builder's own failure evidence — check they show what a failing
   plant looks like under the FINAL code, or that the final calibration would still catch a
   regression).
4. **Are the round-2 verdicts trustworthy where they CHANGED from round 1?** The build report's
   delta table claims: all 1,019 round-1 parse failures eliminated; large drops in spurious and
   omitted verdicts; EmptySet/zero-sample decisions demoted to undecided-confirmed. Sample at
   least five changed verdicts per census (containment WL, containment PY, probe WL, probe PY)
   and re-decide each with your own independent script from the record's own operands. A verdict
   that changed WITHOUT a corresponding repair class is a finding of the first order.
5. **Reducer**: recompute its counts (`889 failures / 358 findings / 885 limitations /
   ROUND_FAIL`) from the round-2 census stdout files with your own script; any discrepancy is a
   finding. Verify the taxonomy is now CLOSED: every verdict token the four census stdouts
   actually contain maps to exactly one bucket, and an unrecognized token fails acceptance.

## What you are handed

The repair directive and brief (paths above); the spec
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (§5, Q8a/Q8b);
both committed records (WL `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` at
`a4cf6539`, PY `scripts/out/S11_stray_longitudinal_sympy_audit.out` at `19591194`); round-1 census
outputs under `/home/trevnorris/.s11_build/census_build/` (for delta verification); round-2
outputs under `/home/trevnorris/.s11_build/census_build2/`; `DEFECT_REGISTER.md` (7 entries,
deliberately untouched — a census rediscovering one is corroboration, not a finding).

## Required method

These artifacts are SCRIPTS. Write your own re-decision scripts BEFORE reading the repaired code,
from the brief, the spec's locus protocol, and the committed records alone; save every script and
its literal stdout to named absolute paths under `/tmp/` and report those paths — a prose
derivation is worth nothing and will be discarded.

⛔⛔ **ABLATION IS MANDATORY.** For each repair class you verify, ablate it in a /tmp copy of the
module: revert or corrupt the repaired predicate ONE class at a time — e.g. restore the flipped
verdict token, break the new sheet rule, drop the premise conjunction — re-run against a small
slice of the real record, and report the literal verdict diff. A corruption that moves nothing is
a finding of the first order: it means the repair is dead code and the verdict changed for some
other reason. Also run one-sided corruption on any two-route comparison the diff touches.

Probe for: a verdict verified using the object that produced it; an expected value living inside
the instrument or its calibration; an `assert` before the print it guards; a per-record line
whose verdict token does not follow from the printed operands; a calibration plant that cannot
fail under the final code.

## Physics filter

Report a finding only if it catches a way a round-2 census verdict could be wrong, a real record
class could be silently unmeasured, or a certified-correct verdict direction was weakened. Style,
performance, and hypothetical-input robustness are out.

## Sandbox

Copy anything you modify to /tmp and ablate the copy; never modify the working tree; never
commit. No Wolfram kernel is needed or permitted — everything is Python + the record files. Cap
any single script run at `timeout 600`; a timeout is a failed probe to report, not a reason to
retry harder. The full censuses take hours — run your re-decisions on slices, never the full
population.

## Report format

Numbered findings, most severe first: claim, script path, literal stdout excerpt, instrument
file:line (at `89ed80c9`). If nothing survives the filter, list the ablations, samples, and
literal outputs you ran instead.
