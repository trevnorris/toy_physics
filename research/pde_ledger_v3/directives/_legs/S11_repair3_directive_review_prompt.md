# Independent review — round-3 census-instrument repair directive (decision list)

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_repair3_directive.md`
at commit `61eec069`, with its measurements twin in `../_measurements/` (same name).

## What to check

This is an orchestrator-written decision list. A builder will implement it verbatim against the
instrument modules at `89ed80c9`; every ambiguity or wrong requirement becomes a build round. The
five defect claims were verified by two prior review legs and orchestrator re-computation (the
twin carries the commands); your job is the DIRECTIVE, not re-litigating the defects:

1. **Is each required semantics CORRECT?** For each of the five defects: does the "Required"
   sentence actually eliminate the measured failure without introducing a new one? Check each
   against the actual code at `89ed80c9`
   (`research/pde_ledger_v3/reduction/s11_census_math.py`, `s11_containment_census.py`,
   `s11_acceptance_reducer.py`) and against real record payloads — e.g. is per-branch
   membership with undefined-branch-as-non-covering sound on a candidate that is undefined on
   EVERY emitted branch? Is "pairwise-distinct within each assumption pool" sufficient for the
   refutation direction, or can a coincidence residual still evade a distinct-value sampler?
   Does the witness partition (bound premises decide, unbound printed verbatim) mis-handle any
   real witness in either committed record — e.g. one binding SOME mode variables?
2. **Is anything under-specified** — a requirement the builder can satisfy two materially
   different ways with different census verdicts?
3. **Does any requirement weaken a certified direction?** The closing constraint lists what must
   stay. Check no "Required" sentence quietly contradicts it.
4. **Level-above misses**: a defect interaction the directive does not order correctly (e.g. do
   defects 1 and 2 interact — a pole candidate is also the poisoning sibling; which rule wins
   and is the directive unambiguous there?); a calibration plant that cannot fail as specified;
   a plant whose specified behavior encodes an expected census outcome the builder could
   iterate toward.
5. **Count statements**: the directive states measured round-2 facts as evidence. Verify each
   against `/home/trevnorris/.s11_build/census_build2/` outputs (commands in the twin). A wrong
   count in a directive breeds a wrong acceptance surface.

## What you are handed

The directive + twin; the round-2 repair directive
(`S11_census_instr_repair_directive.md`, `33babf8d`) it supersedes in part; the folded brief
(`S11_locus_census_instrument_brief.md`, `fa8c58b3`) whose obligations it cites; the instruments
at `89ed80c9`; the round-2 census outputs under `/home/trevnorris/.s11_build/census_build2/`;
both committed records (WL `a4cf6539`, PY `19591194`); `DEFECT_REGISTER.md` (untouched, 7
entries).

## Required method

Where a directive claim is checkable by computation (a count, a code-line citation, a defect
reproduction, a semantics counterexample), check it by computation and show the literal command
and output; save scripts and stdout to named absolute paths under `/tmp/` and report the paths.
A prose objection to a computable claim will be discarded. For semantics questions (1) and (4),
construct a concrete counterexample payload (byte-shaped like the records) if you claim a
required semantics is wrong — a hypothetical is not a finding.

## Physics filter

Report a finding only if it changes what the builder would build, what the census may verdict,
or what the round may claim. Style, tone, and formatting are out.

## Sandbox

Read-only on the working tree; write only under /tmp; never commit. No Wolfram kernel. Cap any
single script run at `timeout 600`.

## Report format

Numbered findings, most severe first: claim, the directive line it targets, evidence (command +
literal output or concrete counterexample payload). If nothing survives the filter, list what
you checked with literal outputs.
