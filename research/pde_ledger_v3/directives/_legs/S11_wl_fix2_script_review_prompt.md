# Independent physics review — S11 WL engine after fix round 2

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
at commit `d6316978` (the working tree matches this commit byte-for-byte; verify with
`git -C /var/projects/toy_physics diff --stat d6316978 -- research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`).

## What to check

This engine is the Wolfram half of a dual-engine audit of a linearized brane action: it computes,
per control package and dimension, the mode spectrum, polarization structure, and the stratified
degeneracy loci of the dispersion variety. Fix round 2 rebuilt every CAS call class on the live
cell path into uniformly bounded routes (time AND memory per attempt, measured-expiry fallbacks,
honest undecided tokens carrying attempt provenance). The decision list that governed the repair is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round2_brief.md`
with its generated measurements twin in `../_measurements/`; round-1 obligations 1, 5, 6, 7 bind it
(`S11_wl_engine_fix_round1_brief.md`). The emitted-object spec is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (locus protocol §5;
solution sets emitted "exactly as your CAS returns it", spec:245; opaque non-answers stay opaque,
spec:271–272).

Check whether the physics this engine emits can be wrong: a bounded route that silently narrows a
solution set, an undecided token where the object is decidable, a fallback selected by cell identity
instead of measured expiry, an emitted object that is typed rather than computed.

## Three open measured findings you must attack with computation

1. **Joint-admissibility arity change (obligation-6 HALT).** Between the round-1 partial record
   (committed canonical `.out`, commit `dfc2568c`) and the new complete run, the token lists of
   `WL_S11_XKIN_ANISO_D4_ROOT2_RANK_DROP_JOINT_REAL_ADMISSIBLE` and `..._ROOT3_...` changed arity:
   44 candidate points (mixed ADMISSIBLE/EXCLUDED/UNDECIDED) → 16 and 12 points (all decided). The
   underlying `..._RANK_DROP_JOINT_..._SOLUTION` branch lists changed shape between the old unbounded
   `Solve` and the new bounded route. Determine by computation whether the new branch list and the
   old one describe the same solution set — an omitted branch is the failure mode that matters.
   Both record versions are in the files listed below.
2. **Neither obligation-4 instrument has ever run on real records.** The build's SymPy containment
   harness (`/home/trevnorris/.s11_build/fix2_build/sympy_completeness.py`) parse-fails on all 484
   equations/solution locus pairs in the real records (log:
   `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/sympy_completeness_final.log`),
   and the WL probe census (`/home/trevnorris/.s11_build/fix2_build/probe_census.wl`) dies on
   `ToExpression::sntxi` before emitting a single probe line while still exiting 0
   (log: `/home/trevnorris/.s11_build/fix2_final/probe_census_final.log`). Both calibrated only
   against planted synthetic records, never the real file format.
   Write your own independent containment check for at least the two HALT loci of finding 1 and any
   three further loci of your choosing: solve the emitted `_EQUATIONS` under the record's premises in
   an engine that is not Wolfram, and test membership of every branch you find in the emitted
   `_SOLUTION`. Save the script and its literal stdout.
3. **`"NOT_APPLICABLE"` solution slots.** Many `_LOCUS_SOLUTION` records carry the string
   `"NOT_APPLICABLE"` while their paired `_LOCUS_EQUATIONS` record carries live equations. Decide by
   computation on sampled cases whether those equation systems have solutions under the record's
   premises. Read the spec's locus protocol first and report whether this emission is
   spec-sanctioned emptiness or a coerced non-answer (spec:245, 271–272).

## What you are handed

- The artifact (path above) and the spec `S11_SHARED_PHYSICS.md`.
- The round-2 brief + measurements twin, and the round-1 brief it binds.
- Final-engine complete records (produced by the committed engine bytes):
  `/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out` (1,898 emissions),
  `.../guarded_XKIN_ANISO_D3.out` (2,501), `.../guarded_XKIN_ANISO_D4.out` (5,412).
- The committed canonical WL record (19 cells + round-1 partials):
  `git show dfc2568c:research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`.
- The SymPy sibling's committed record for cross-engine reference:
  `/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`
  (note: its `STRATUM_ORDERING` is empty on every cell — the strata surface is WL-only by
  measurement; do not report that as a defect).
- `/var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md`: the three entries dated
  2026-08-16 are registered and deliberately out of scope this round — do not re-report them; the
  register's TEXT, not its line numbers, identifies them.
- Measured and already explained, do not re-report: (a) between repeat runs, a near-budget attempt
  may expire on its time budget in one run and its memory budget in another — the decided values
  agree and only the recorded `Failure[...]` head differs; (b) `Q7_CURL*` families exist only in
  D=3 in every package (3D cross product), so their absence from D2/D4 is structural.

## Required method

This artifact is a SCRIPT. Write your own derivation script **before** opening the artifact, and
save **both the script and its literal stdout** to named absolute paths under your scratch
directory. Without these, your derivation claims will be discarded. Derive independently; a prose
derivation is worth nothing.

Ablate every load-bearing check and report its literal output. Probe for: a value verified using
the predicate or definition that produced it; a conclusion emitted as an unconditional literal; a
check whose expected value lives inside the artifact it checks; an `assert`/`Abort` that precedes
the value it guards. Ask of every claim: WHICH LINE COMPUTED THIS — give the line number or report
it as uncomputed.

⛔⛔ **A FORM ABLATION IS MANDATORY.** Change the **structure** of a load-bearing object — flip a
sign *and* an off-diagonal in the action's kinetic block, collapse two independent symbols into
one — re-run a cheap cell, and report the **literal** diff of the emitted tags. A coefficient
rescale tests arithmetic; only a form change tests physics. Additionally ablate the ROUTE layer:
set one primary budget to 0 in the copy and verify by literal output that the fallback fires
**uniformly at every call site of that helper** and that the emitted record carries the primary's
actual expiry — a fallback reachable only for some cells' identity is the defect class this round
existed to remove.

## Physics filter

Report a finding only if it catches a way the physics could be wrong; do not report "the script
would be wrong on a different input."

## Ablation sandbox and kernel discipline

⛔ Copy the artifact to /tmp and ablate the COPY. Never modify the working tree, never commit.
⛔ Wrap EVERY kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it and move on.
⛔ NEVER raise the timeout, and never run more than one kernel at a time (the licence has TWO seats
and one must stay free).
⛔ Do NOT run `XKIN_ANISO` dimension 2 or 4 at all — they need 16–77 minutes and up to 25 GB. Use
`MAIN 2` (~8 s), `XFORM_CURLONLY 2` (~7 s), or `XKIN_ANISO 3` only if you must touch the package
(~260 s), or synthetic operands.
Invocation contract: `wolframscript -file <abs .wl> <PACKAGE> <D>`.
Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Report format

Numbered findings, most severe first. For every finding: the claim, the script path, the literal
stdout excerpt that shows it, and the artifact line number(s). If nothing survives the physics
filter, say so explicitly — and list the ablations you ran with their outputs.
