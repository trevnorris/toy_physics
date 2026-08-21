# Directive-design review — S11b cross-engine comparator build directive

You are one of two independent legs reviewing an ORCHESTRATOR-written build directive BEFORE any builder
runs against it. Your job: find defects in the DIRECTIVE — leaks, under/over-specification, a contract
mis-transcription, an acceptance test that can't fail, or a mechanism that would let the comparator hide the
very cross-engine disagreements it exists to surface. A leg that finds nothing is weak evidence; if you
find nothing, show the specific checks you ran.

## The artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_comparator_build_directive.md`

## Sources of truth — read these and form your own view FIRST, then read the directive
- The frozen T7 comparator contract: `research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md`
  (search "T7") — especially "reject a native boolean as a residual operand" and three-valued undecided.
- The §10 tag grammar: `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` §10 (both engines emit
  parallel `<ENGINE>_S11B_<QUANTITY>` tag sets; `_LOCAL_` infix; one tag per named object).
- The mechanical precedent the directive says to reuse:
  `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py` (957 lines). ⚠ Note its `residual()` uses
  `py_value == wl_value` — check whether that is safe once native booleans are in the stream (S10 had none;
  S11b has 108 PY + 26 WL boolean payloads).
- What the comparator must be capable of surfacing (do NOT let the directive suppress these):
  `research/pde_ledger_v3/steps/S11b_wl_engine_review_disposition.md` — a genuine cross-engine SIGN
  disagreement on one object, and a genuine energy-basis COUNT disagreement (`ENERGY_BASIS_COUNT`).
- The two real engine transcripts exist but are the ORCHESTRATOR's frozen run, not the builder's input:
  `scripts/out/S11b_interface_coupling_law_sympy_audit.out`,
  `mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`.

## Required checks — report a finding only if it catches a way the COMPARATOR could be wrong
1. **Rule-5 leak.** Does the directive state, anywhere, what the real engines equal, whether they agree, an
   expected residual, the contested sign, or the basis count? Quote any leak.
2. **Acceptance decisive & able-to-fail.** For each of the 7 synthetic acceptance tests, is there a concrete
   BAD implementation that would PASS it? An acceptance that a defective comparator also passes is worthless.
   In particular: does the native-boolean test actually catch the `False`-vs-`0`-scores-AGREE bug, and does
   the repoint-ablation test actually force a DISAGREE?
3. **T7 contract fidelity.** Is native-boolean rejection specified so it covers ALL forms the parsers
   produce (Python `bool`, SymPy `S.true`/`S.false`, a Wolfram `True`/`False` token)? Is three-valued
   undecided preserved (never forced binary)? Is the join by object name, residual of paired payloads?
4. **Does any mechanism let the comparator DESIGN AWAY a disagreement?** Scrutinize the `W_0` rule: could
   "derive `n` from the `DIM_` tags and normalize" ever turn a genuine sign/structure disagreement into
   AGREE? Is it correctly bounded to dimension-licensed cases and flagged (never silent)? Could the
   Association key-matching or symbol transliteration mask a real difference?
5. **Blindness vs buildability.** The directive says the builder must not tune residual/normalization to a
   real payload, yet must parse real payload SHAPES (`<|…|>`, radicals). Is that boundary specified well
   enough to be followed, or is it ambiguous in a way that invites tuning toward agreement?
6. **Under/over-specification.** Anything the builder would trip on, or any parallel machinery the directive
   demands that the frozen contract does not need.

## Method
This is a DOCUMENT review. Quote both the source of truth and the directive for every finding. You may run
small scripts to check a claim (e.g. does `parse_mathematica` ingest `<|…|>`? does `sympify("False")==0`?);
if you do, save the script + its literal stdout to named absolute paths and report them — a prose claim
about code behaviour is discarded. ⛔ Do not edit the working tree.

## Report — under ~25 lines
Numbered findings, each: the directive line, the source-of-truth it violates or the failure it permits, and
the concrete fix. Then a one-line verdict: is the directive safe to build from, or must it be folded first?
