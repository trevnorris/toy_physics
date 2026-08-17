# Independent review — S11 census instruments (scripts) and their first real census

## Artifact

The six instrument modules at commit `90ab5e2d` (working tree matches):
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/s11_census_common.py`,
`s11_census_math.py`, `s11_containment_census.py`, `s11_undecided_probe_census.py`,
`s11_calibrate_censuses.py`, `s11_acceptance_reducer.py` — plus the census outputs they produced
under `/home/trevnorris/.s11_build/census_build/` (sha256s pinned in
`/home/trevnorris/.s11_build/census_build/build_report.md`).

## What to check

These instruments adjudicate physics: their verdicts (spurious branch, omitted branch, witness
failure, decided-undecided) will drive the engine-repair round. The governing decision list is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_locus_census_instrument_brief.md`
(folded, commit `fa8c58b3`) with its generated twin in `../_measurements/`. The first real census
reported: 1,019 in-population parse failures, 971 decided-undecided records, 217 spurious
branches, 355 omitted-branch records, 46 witness failures. Every one of those numbers is only as
good as the code that computed it. Check whether a verdict can be wrong:

1. **Spurious-branch logic**: does it substitute the branch as written (never sign-flipping
   radicals inside the branch), evaluate every equation under every global sheet of the remaining
   radicands, and only then emit the verdict? Sample at least five `SPURIOUS` per-record lines
   from each engine's census output and re-decide them with your own independent script.
2. **Omitted-branch logic**: is the independent solve genuinely independent of the producer's
   route (the PY producer uses `sp.solve`; an "independent" census calling the same route
   reproduces the producer's omissions and cannot see them)? Is containment decided up to
   algebraic equivalence, with sampling clearly tokened when used? Sample at least five `OMITTED`
   verdicts and re-decide them independently.
3. **Decided-undecided logic**: 971 is a large claim. Sample at least five from each engine and
   verify the probe's decision by your own computation on the record's own operands.
4. **The 1,019 parse failures**: classify what actually fails to parse (the per-record lines name
   each). Is any "parse failure" hiding a payload the brief's population section requires the
   instrument to handle (attempt-provenance lists, `ConditionSet`, multi-megabyte lines)? Is any
   in-population record silently absent from the per-record lines entirely (count them yourself)?
5. **The reducer**: recompute its counts from the census stdout files with your own script; any
   discrepancy is a finding. Check that findings are carried as per-record lines, not summary
   bits, and that nothing in the failure/finding split lets a real defect hide.
6. **Calibration**: could it pass while the census path diverges? Verify the planted records
   traverse the same entrypoint/parser/solver as the real census (code path, not comment).

## Required method

These artifacts are SCRIPTS. Write your own derivation/verification scripts BEFORE reading the
instrument code, from the brief, the spec's locus protocol
(`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` §5,
Q8a/Q8b), and the committed records alone; save every script and its literal stdout to named
absolute paths and report those paths — a prose derivation is worth nothing and will be
discarded.

⛔⛔ **ABLATION IS MANDATORY.** For an instrument the form ablation is verdict-logic ablation:
in a /tmp copy, corrupt ONE load-bearing test at a time — flip the membership predicate, break
the sheet consistency, disable the branch-as-written rule — re-run against a small slice of the
real record, and report the literal verdict diff. A corruption that moves nothing is a finding
of the first order. Also run one-sided corruption on any two-route comparison you find.

Probe for: a verdict verified using the object that produced it; an expected value living inside
the instrument; an `assert` before the print it guards; a per-record line whose verdict token
does not follow from the printed operands.

## What you are handed

The brief + twin; the spec; both committed records (WL at `a4cf6539`, PY at `19591194`); the
census stdout files and build report under `/home/trevnorris/.s11_build/census_build/`;
`DEFECT_REGISTER.md` (7 entries, deliberately untouched — the census may rediscover them; that
is corroboration, not a finding against the instrument).

## Physics filter

Report a finding only if it catches a way a census verdict could be wrong or a real record class
could be silently unmeasured. Style, performance, and hypothetical-input robustness are out.

## Sandbox

Copy anything you modify to /tmp and ablate the copy; never modify the working tree; never
commit. No Wolfram kernel is needed or permitted — everything here is Python + the record files.
Cap any single script run at `timeout 600`; a timeout is a failed probe to report, not a reason
to retry harder.

## Report format

Numbered findings, most severe first: claim, script path, literal stdout excerpt, instrument
file:line. If nothing survives the filter, list the ablations and samples you ran with their
literal outputs.
