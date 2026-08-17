# Independent review — S11 locus-census instrument decision list

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_locus_census_instrument_brief.md`
at commit `68f458e6`, with its generated measurements twin
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11_locus_census_instrument_brief.md`.

## What to check

This decision list governs the next build: two census instruments (an independent
containment/soundness census over emitted `_EQUATIONS`/`_SOLUTION` locus pairs, and an
undecided-record probe census) that must run on the real committed records of both engines. The
brief is orchestrator-written and the builder will trust it unreviewed — an error here makes the
build converge on the wrong thing. Check it the way you would check a spec: every obligation must
be satisfiable, testable, and sufficient; a repair satisfying the letter of every numbered item
while leaving the measured wall in place is a defect in the brief, and you should construct such a
defective repair explicitly for any obligation where one exists.

Specifically probe:
- Can an instrument pass obligation 1 (parse 100%) by classifying hard lines into a benign
  category instead of parsing them? Is "the record classes it audits" tight enough?
- Does obligation 2's calibration actually bind the census that follows — could the planted
  records exercise a code path the real census never takes?
- Is obligation 3 consistent with an instrument whose verdicts must feed a later round's
  acceptance? Does "exit code reflects execution, never verdict content" create a
  green-while-failing hazard anywhere downstream?
- Obligation 4 names sheet-consistency for radical systems. Is a global-sheet sweep sufficient
  for the record classes in the twin, or does some emitted object class need per-branch sheet
  data the brief fails to demand?
- Obligation 5's no-starvation floor cites the engine budgets in the twin. Are there record
  classes whose emitting route had NO budget (so the floor is undefined)?
- Is anything in the brief an acceptance criterion referencing an expected value of the real
  census (it must not be)? Is anything missing that the round-2 brief's obligation 4 demanded
  (residual recompute, live operands) that this round silently drops?

## What you are handed

The brief and twin (paths above); the round-2 brief
`S11_wl_engine_fix_round2_brief.md` (its obligation 4 is the parent obligation);
`/var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md` (7 entries — this round
builds the instrument, repairs nothing);
the committed records named in the twin (WL at `a4cf6539`, SymPy at `19591194`);
the spec `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`
(locus protocol §5, Q8a/Q8b);
the old broken instruments `~/.s11_build/fix2_build/sympy_completeness.py` and
`~/.s11_build/fix2_build/probe_census.wl` (as the measured failure, not as a base to repair).

## Required method

This artifact is a DOCUMENT. Read the sources of truth first — the spec's locus protocol, the
register, the twin, and sample lines of both committed records — form your own view of what a
working census must do, and only then read the brief. Quote both sides for every finding. Where a
claim in the brief or twin is checkable by a command (record counts, pairing claims, budget
constants), run the command and report its literal output; a prose disagreement without a command
behind it will be discarded.

## Physics filter

Report a finding only if it catches a way the built instrument could pass this brief while
mismeasuring the records — a wrong verdict, a silently skipped class, an uncatchable defect. Do
not report style, redundancy, or "a different structure would be cleaner."

## Sandbox

Read-only review. Copy anything you want to experiment on to /tmp; never modify the working tree;
never spawn a Wolfram kernel (the records are text; nothing here needs one).
