# Independent review: a decision list for a Wolfram engine fix round

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round1_brief.md`
with its generated measurements twin
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11_wl_engine_fix_round1_brief.md`
(both at commit `9190d569`).

## What to check
The artifact is an orchestrator-written decision list that will govern a builder's repair of the
Wolfram engine `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`.
The engine's XKIN_ANISO D4 cell was killed after 15+ hours of CPU-pegged silence inside a
quantifier-elimination call, and D2 previously died on memory inside its strata; the brief's
obligations are supposed to make a repair that (a) terminates, (b) computes exactly the objects the
shared spec names, and (c) cannot be greened by a defective repair. The decision list is the ONE
artifact the builder trusts unchecked — a defect here multiplies into a full build round.

## What you are handed
- The artifact and its twin (above).
- The spec the objects come from: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`
  (the locus protocol is §5, lines ~231–295).
- The engine source named above, and its committed output
  `research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`.
- The SymPy sibling engine `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
  and its committed output beside it in `scripts/out/` — the independent construction of the same
  spec, useful as a comparison point for what the records look like when the protocol completes.

## Required method (document branch)
Read the SOURCES OF TRUTH first — the spec's locus protocol, the engine's `emitLocus`, the
committed `.out` record — and form your own view of what a sound fix-round decision list for this
wall must require. Only then read the artifact, and quote both sides for every finding.

Then attack the artifact along these lines, with computation wherever a claim is checkable:

1. **Construct a defective repair that satisfies every numbered obligation.** This is the highest-
   value finding class: the orchestrator's acceptance bars have been greened by defective repairs
   six consecutive times in this program. Try at least: an ops-gated route that silently emits the
   undecided token everywhere it used to hang; a repair that special-cases exactly the cells the
   bars measure; a repair that satisfies the byte-identity bar by never changing anything and the
   completion bar by weakening a stage the bars do not census. If any such repair passes every
   obligation as written, report it as BLOCKING with the obligation text it slips through.
2. **Check every measured claim in the brief against the twin, and the twin's commands against the
   artifacts they cite.** Re-run any command you doubt (they are all read-only). A claim in the
   brief with no command behind it in the twin is itself a finding.
3. **Check the obligations against the spec.** Does any obligation permit an emitted object weaker
   than the spec's definition (tokens, payload field sequences, branch completeness, the
   computed-from-_EQUATIONS requirement at spec line ~251)? Does any obligation accidentally
   prescribe a recipe where the spec names only an object?
4. **Check internal consistency.** Can obligations 1, 4 and 5 be satisfied simultaneously by at
   least one honest repair? If every honest repair must violate one of them, that is BLOCKING.
5. **Check the budgets and bars for fabrication forcing.** Is there any cell state whose only
   honest outcome fails a bar as written (so a builder iterating to green must fabricate)?

## Physics filter
Report a finding only if it catches a way the physics record could be wrong or a way the build
round could be wasted; do not report "the document could be prettier."

## Constraints — these bind you
- ⛔ Do NOT launch any Mathematica kernel or `wolframscript` process — the licence seats are
  reserved. All computation you do runs in Python/SymPy.
- ⛔ Never modify the working tree. Scratch work goes under a directory you create in /tmp.
- Wrap every script run in `timeout 600`; a timeout hit is a data point, not a licence to retry
  bigger.
- ⭐ For every computed claim in your report, save the script AND its literal stdout to named
  absolute paths and cite them; unbacked derivation claims will be discarded.

## Report format
Verdict line, then BLOCKING findings (each: obligation text quoted, defect, evidence), then
non-blocking findings, then what is solid. Attribute every finding to a numbered obligation or to
a specific brief/twin line.
