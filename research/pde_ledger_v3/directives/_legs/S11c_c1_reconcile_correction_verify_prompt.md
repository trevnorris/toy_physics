# Verify — S11c-c1 reconcile/step-record CORRECTIONS (did the orchestrator fold the 2 legs faithfully?)

Two independent legs (Grok + Codex sol) reviewed an S11c-c1 step record that claimed the two blind engines
"AGREE on the curved-interface bulk closure." They found the KERNEL agreement real but the CLOSURE-WIDE verdict
overstated. The orchestrator then FOLDED their findings and committed the corrections (`git show 072d0b75`). Your
job: confirm the fold is **faithful** (both legs' findings incorporated), **correct** (the new verdicts match what
was actually measured), **internally consistent** (no contradiction left across the artifacts), and **not
over-corrected** (nothing sound was needlessly downgraded, no NEW error introduced). ⛔ A prose judgement is worth
nothing where a claim rests on computation — RE-RUN it and report the LITERAL stdout with the command. Save scripts
+ outputs to named `/tmp` absolute paths and report them. ⛔ Copy anything you ablate to `/tmp`; never touch the
working tree.

Working dir `/var/projects/toy_physics` (branch `ledger-v3-rebuild`). Paths below under `research/pde_ledger_v3/`.

## The corrected artifacts (read these)
- `steps/S11c_c1_curved_bulk_closure.md` (the step record)
- `directives/_measurements/S11c_c1_comparator_reconcile.md` (the reconcile adjudication)
- `DEFERRED_HEAVY_RUNS.md` (the S11c-c1 section near the end)

## The two legs' reports the fold must honor (read these — the source of the required corrections)
- `directives/_measurements/S11c_c1_step_record_review_grok.txt`
- `directives/_measurements/S11c_c1_step_record_review_codex_sol.txt`

## The measurements behind the corrected verdicts (re-run to check them)
- Reproducible bridge (reads the committed `.out`; `datalad get` them first — `scripts/out/S11c_c1_bulk_closure_sympy_audit.out`,
  `mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out`):
  `python3 directives/_measurements/S11c_c1_reconcile_reproduce.py --py <PY.out> --wl <WL.out> --adversarial`
  (or `--runlog <a comparator run over DTN_KERNEL + FACE_RESPONSE_COEFFS>`).
- The engine sources: `scripts/S11c_c1_bulk_closure_sympy_audit.py`, `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl`.

## Checks (CLEAR / DEFECT, each with the command + literal stdout)
1. **Faithful fold.** Every MUST-FIX / defect the two legs raised is addressed in the corrected artifacts, and NOT
   silently dropped. Grok's DEFECT: the ESTABLISHED bullet claimed the raw DtN OPERATOR AGREE (reconcile left it
   UNDECIDED) — is the kernel now split from the operator, and is c2 told not to inherit an operator-AGREE? Codex
   sol's 4 findings: (1) rule-17 seal-5 mis-adjudication → is seal 5 now UNDECIDED (surfaced freeze), not AGREE?
   (2) deferred/uncovered physics pre-adjudicated AGREE → are the 4 giants now UNMEASURED (not "adjudicated
   AGREE") in `DEFERRED_HEAVY_RUNS.md`? (3) "coefficients proven EXACT ZERO" overstated → is the coefficient claim
   now scoped (LAB_HELD off-diagonal + MATERIAL on-diagonal; flat-resolvent leg-labeling surfaced)? (4) bridge
   scripts not self-contained → does `S11c_c1_reconcile_reproduce.py` run from the committed `.out` with no
   scratchpad dependency? Re-run it and report.
2. **Correct, not over-corrected.** Is the KERNEL AGREE still asserted at full strength (all 4 cases exact zero,
   corruptions bite)? Re-run the reproduction script and confirm. ⛔ The fold must NOT have downgraded the kernel
   result (which the legs AFFIRMED). Conversely, seal 5 / operator / giants / ENERGY / flat-leg must be UNDECIDED
   or DEFERRED, not re-asserted as AGREE anywhere.
3. **Seal-5 resolution.** The two legs DISAGREED on seal 5 (Grok: AGREE-with-caveat; Codex: keep surfaced). The
   fold adopted "UNDECIDED — surfaced rule-17 freeze; c2 re-adjudication mandatory." Is that the correct resolution
   of a GOVERNING-rule disagreement (rule 17: a required freeze is the finding; ⛔ do not fold a live field to a
   constant to declare AGREE)? Independently confirm neither engine differentiates the density (grep both).
4. **Internal consistency.** No artifact still says "the engines AGREE on the closure" / "closure-wide AGREE" /
   "already adjudicated AGREE" for the deferred families. The step record's STATUS box, seal list, ESTABLISHED/OWED
   split, method notes, and What's-next must all agree with the reconcile record and `DEFERRED_HEAVY_RUNS.md`.
5. **c2 blast radius.** Does the corrected record clearly tell c2 (a) the kernel + flat symbol + response are
   cross-engine AGREE and safe to consume, but (b) `dtn_operator` whole-form and the density seal are NOT
   cross-engine-closed and MUST be re-adjudicated in c2? A future c2 author must not be able to read an
   operator/density AGREE out of these files.

## Output
Per check CLEAR/DEFECT with the command, literal stdout, and your script/output paths. Rank real defects
most-severe first; a clean pass with citations is equally useful. ⛔ Do NOT propose an expected measured residual
for a deferred/UNDECIDED family. If the fold is faithful and correct, say so plainly — we are deciding whether
S11c-c1 is in the clear to move to c2.
