# Independent review: the decision list for WL engine fix round 2 (strata memory wall)

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round2_brief.md`
at commit `af6be9e5`, with its generated measurements twin
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11_wl_engine_fix_round2_brief.md`.

## What to check

This decision list will govern a builder repairing two measured, deterministic memory deaths in
the Wolfram engine. The list is the one artifact the builder trusts unreviewed — an error here
makes the builder and every downstream review agree on the same wrong thing. Review the LIST, not
the engine. Specifically:

1. **Is every obligation well-defined, decidable from artifacts, and able to fail?** An
   acceptance that cannot fail for a bad repair is the defect class that has cost this project
   the most.
2. **MANDATORY defective-repair construction.** Construct, precisely (pseudocode or exact edit
   description), at least one repair that satisfies EVERY obligation as written while leaving
   the emitted record wrong or dishonest, or the wall effectively alive. Rounds 1–6 of this
   program each had a greenable acceptance found only by this probe. If you cannot construct
   one, report what stopped each attempt — the attempts are the evidence.
3. **Contradictions**: within this brief; between it and the round-1 brief sections it binds
   verbatim (`S11_wl_engine_fix_round1_brief.md`, obligations 1, 5, 6, 7 and the Post-legs
   fold); between it and the spec's locus protocol (`directives/S11_SHARED_PHYSICS.md` §5); and
   between it and the committed guarded death artifacts it cites.
4. **Level-above misses**: given the measured wall (read the death records yourself), what
   obligation is MISSING? In particular: does the list actually force the D2 accumulation
   mechanism to be fixed, or only the terminal call? Does any obligation quietly permit an
   undecided flood, an identity-gated route, or a weakened record?
5. **Twin audit**: does every measured claim in the brief have a backing command block in the
   twin; is any twin block hand-typed rather than command-backed; does any block leak a record
   value for the two open cells (their locus solutions, ranks, nullities, counts, or status
   tokens are deliberately withheld)?

## What you are handed

The brief and twin above; the round-1 brief and its twin (same directories, `round1` names); the
engine `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`; the spec
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`; the defect register
`research/pde_ledger_v3/DEFECT_REGISTER.md`; the death artifacts under
`/home/trevnorris/.s11_build/fix1_build/` (`guarded_cells_record.log`,
`guarded_XKIN_ANISO_D{2,3,4}.out`, `.emission_times.tsv`, `.progress.log`, and the `.prerepair2`
/ `.killed_run1` variants); full read access to the repository.

## Required method

Document branch: read the sources of truth FIRST — the engine at the cited lines, the death
records, the round-1 brief, the spec — and form your own view of what a correct round-2 decision
list must require; only then read the artifact, and quote both sides for every finding. Where a
brief claim is checkable by computation (emission counts, gaps, tag order, source lines), check
it with a command and report the command and its literal output — a prose re-derivation is worth
nothing.

## Operational constraints (identical for both legs)

- ⛔ Do NOT launch any Mathematica kernel or `wolframscript` process for any purpose. This
  review needs none; the licence has two seats and they are reserved.
- ⛔ Never modify the working tree. Write any scripts/outputs under a directory you create in
  /tmp (`mktemp -d /tmp/s11_fix2rev_XXXX`) and cite absolute paths.
- Wrap every script run in `timeout 600`.

## Physics filter

BLOCKING: a defect that would let a build satisfy every obligation while the emitted record is
wrong or dishonest, the wall survives in practice, or a genuine disagreement gets silently
absorbed. Everything else is non-blocking.

## Report format

Verdict line first; BLOCKING findings with literal evidence (quote the obligation text and the
defective-repair construction or contradiction); non-blocking; the defective-repair
construction(s) in full; commands + literal output for every checked claim.
