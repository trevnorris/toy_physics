# Independent review — S11 census-instrument round-3 repair diff and its census

## Artifact

The repair diff `git diff 89ed80c9..fd9a5835 -- research/pde_ledger_v3/reduction/` (four modules;
working tree matches `fd9a5835`), the modules themselves under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`, and the round-3 census outputs
under `/home/trevnorris/.s11_build/census_build3/` (sha256s in `artifact_sha256.stdout`; build
report `build_report.md`).

## Scope — SCOPED review of the round-3 diff

Full reviews already ran on the round-1 instruments and the round-2 diff; the round-3 acceptance
surface is the folded directive
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_repair3_directive.md`
(commit `ef9085c6` — read it and its `../_measurements/` twin; its six numbered defect classes
were verified by two directive legs plus orchestrator computation), under the brief
(`S11_locus_census_instrument_brief.md`, `fa8c58b3`). Do not re-litigate what earlier legs
cleared outside this diff. Your questions:

1. **Does each of the six repairs implement the directive's Required semantics** in the executed
   census path? Ablate each one class at a time in a /tmp copy against a real-record slice and
   report the literal verdict diff. In particular: is the union test now the defined-branch
   union product (piecewise coverage decided, undefined siblings excluded, `COVERAGE_UNDECIDED`
   never collapsed)? Are premises substituted FIRST, split into conjuncts, classified
   TRUE/FALSE/contingent, with concrete-number assumption atoms deciding and no
   pre-substitution re/im expansion deciding truth? Is the sampler pairwise-distinct within
   each assumption pool?
2. **Did any repair weaken a certified direction?** The directive's closing constraint lists
   them (seven round-2 classes, reducer failure/finding arithmetic, population reconciliation,
   generic NONZERO refutation, `PROVED_CONSISTENT`, two-route residual comparisons).
3. **Is the round-3 calibration real?** Plants through production entrypoints, byte-shaped, and
   each of the nine directive-ordered plants demonstrated ABLE TO FAIL — several are specified
   to fail at `89ed80c9`; verify at least the piecewise-union plant and the radical-realness
   plant actually fail there (run the round-3 calibration case against a `89ed80c9` copy of the
   instruments in /tmp and show the miss).
4. **Are the changed verdicts trustworthy?** The build report claims: WL omitted 13 → 3, WL
   witness failures 172 → 104, WL `WITNESS_VALIDATED` reachable again, probe decided
   561/328 → 589/328, reducer `failures=917 findings=348 limitations=815`. Sample at least five
   changed verdicts per census and re-decide each with your own independent script from the
   record's own operands. A verdict changed WITHOUT a corresponding round-3 repair class is a
   finding of the first order. Check specifically: do the 3 surviving WL omitted records carry
   defined, genuinely-uncovered candidates? Are the 104 surviving witness failures free of the
   pre-substitution-expansion artifact (re-decide the two STRATUM5/6 coincidence witnesses named
   in the directive twin — they must NOT be failures now)? Do the demoted-to-validated and
   demoted-to-undecided witnesses follow from their printed conjunct truths?
5. **Reducer**: recompute failures/findings/limitations from the four round-3 census stdouts
   with your own script (object-level counting per directive defect 5 — sheet lines are
   evidence, not objects); any discrepancy is a finding. Verify the taxonomy stays closed over
   every token actually emitted.

## What you are handed

The round-3 directive + twin; the brief; the spec
(`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` §5,
Q8a/Q8b); both committed records (WL `a4cf6539`, PY `19591194`); round-2 outputs under
`/home/trevnorris/.s11_build/census_build2/` (delta baseline); round-3 outputs under
`/home/trevnorris/.s11_build/census_build3/`; `DEFECT_REGISTER.md` (7 entries, untouched —
census rediscovery is corroboration).

## Required method

These artifacts are SCRIPTS. Write your own re-decision scripts BEFORE reading the repaired
code, from the directive's Required sentences, the brief, the spec, and the committed records
alone; save every script and its literal stdout to named absolute paths under `/tmp/` and
report those paths — a prose derivation is worth nothing and will be discarded.

⛔⛔ **ABLATION IS MANDATORY.** One class at a time, /tmp copy, real-record slice, literal
verdict diff. A corruption that moves nothing is a finding of the first order. Run one-sided
corruption on any two-route comparison the diff touches.

Probe for: a verdict verified using the object that produced it; an expected value living
inside the instrument or its calibration; an `assert` before the print it guards; a per-record
line whose verdict token does not follow from its printed operands; a calibration plant that
cannot fail under the final code.

## Physics filter

Report a finding only if it catches a way a round-3 census verdict could be wrong, a real
record class could be silently unmeasured, or a certified-correct direction was weakened.
Style, performance, and hypothetical-input robustness are out.

## Sandbox

Copy anything you modify to /tmp and ablate the copy; never modify the working tree; never
commit. No Wolfram kernel — everything is Python + the record files. Cap any single script run
at `timeout 600`; a timeout is a failed probe to report, not a reason to retry harder. The full
censuses take hours — run slices, never the full population.

## Report format

Numbered findings, most severe first: claim, script path, literal stdout excerpt, instrument
file:line (at `fd9a5835`). If nothing survives the filter, list the ablations, samples, and
literal outputs you ran instead.
