# Independent review — S11c-c2 N6 ablation-harness directive, FOLD round 2

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_ablation_harness_directive.md`
(a Codex-authored **build directive** = a spec). Review by reading + cheap mechanical checks, NOT execution:
the four harnesses do not exist yet (executable ablation is a later build review). Do not treat their absence
as a finding.

## Context — what this fold fixes
A prior revision (compact output contract + WL single-case budget pin) was two-leg reviewed and found
NOT-SOUND on three findings. THIS fold (round 2) claims to resolve exactly those three, touching only the
output-contract wording, the WL budget case-pin, and the CHANGE-LOG. The **ten-knife design is frozen**
(scope/site, Old, FORM, NO-OP, COEFFICIENT ×2, DEAD sets — byte-identical). The three findings resolved:

- **F1 — WL K_junk was inert.** The prior pin `{"LAB_HELD","RHOBR_CONSTANT"}` sits outside the engine's
  `actualJunkCase = {"MATERIAL_ADVECTED","RHO4_CONSTANT"}` (`audit.wl:21`), and
  `activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0]` (`:845`) is the only `case===` gate,
  so K_junk (`actualJunkCoefficient 0→1`) did nothing under that pin. The fold repins the WL harness to
  `{"MATERIAL_ADVECTED","RHO4_CONSTANT"}`.
- **F2 — `R_N6`/`R_cov`/`SPLIT_CHECK`/guard residuals are heavy PIT objects, not "small scalar direct
  payloads"** (a leak path). The fold removes that carve-out and routes them through the compact PIT
  fingerprint + digest.
- **F3 — the per-object digest is harness-computed, not an engine field.** The fold restates it as a
  harness SHA-256 over the full emitted payload before compaction (WL: full `WL_S11CC2_* = payload` text;
  SymPy: `n.sha` over the parsed emitted PIT payload).

## What you are handed
- The folded directive (path above).
- The fold diff (REV1 → REV2): `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c2_N6_fold_round2.diff`
- The four live engines (read to verify against source):
  - `research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl` (`actualJunkCase`:21, `activeJunk`:845,
    `materialNormalKnife`:860, `sourceChannel`:892, driver `Do`:1011-1014, `drawCount`:760, `ARITHMETIC`:804)
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py` (`ACTUAL_JUNK`:35, `actual_amplitudes` call:187-188)
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py` (`pit`:668, `sha`:91)
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py`
  - the failed-build WL harness (defines the digest mechanism the rebuild reproduces):
    `research/pde_ledger_v3/scripts/S11c_c2_N6_ablation_harness_wl.py` (`emitted_hashes`:190)

## Required method
Verify each claim against the engine source, not the directive's prose. Cite line numbers. Read-only.

## Settle these — reasoning + mechanical check each
**Q1 — F1 fully resolved, nothing else broken by the repin.** Confirm `{"MATERIAL_ADVECTED","RHO4_CONSTANT"}`
= `actualJunkCase` and that under it all three WL knives bite: K_junk (now `activeJunk = actualJunkCoefficient`),
K_carrier (`materialNormalKnife`, not case-gated), K_split_route (`sourceChannel`, not case-gated). Confirm the
three SymPy harness pins are left `LAB_HELD/RHOBR_CONSTANT` AND that every SymPy knife (incl. covariance
`ACTUAL_JUNK`, and diagnostic/reconcile knives) is live on that pin (i.e. no SymPy knife is case-gated to a
different case the way WL K_junk was). Name any knife — WL or SymPy — that is NOT live on its harness's pinned
case.

**Q2 — F2 fully resolved, no residual leak.** Confirm the "direct payload" carve-out for
`R_N6`/`R_cov`/`SPLIT_CHECK`/guard residuals is gone and they now route through compact PIT + digest. Confirm
the surviving "direct payload" allowance (a genuine scalar integer/rational with no ARITHMETIC/table payload)
cannot be satisfied by any heavy engine object. Is there any remaining path — in clause 1, 4, 5, 8, the
Deliverables sentence, or the per-harness certified/DEAD lists — by which a full symbolic object, raw residual
table, arithmetic DAG, PIT sample matrix, or raw engine transcript could reach a committed transcript
(including error paths)?

**Q3 — F3 faithful + implementable.** Is the harness-computed full-payload digest well-defined on both
engines? For WL confirm the driver can SHA the full emitted `WL_S11CC2_* = payload` text (see
`ablation_harness_wl.py:190`); for SymPy confirm `n.sha` applies to the parsed emitted PIT payload and NOT to a
raw DAG `Node`. Critically: does this full-payload digest actually witness a FORM change that does NOT move the
nonzero tally (same-support bite)? If a certified/DEAD object's only witness is this digest, is there any
FORM/COEFFICIENT change in the frozen knife list that would leave BOTH the fingerprint and the full-payload
digest unchanged (a true blind spot)?

**Q4 — no new defect; knives frozen.** From the diff, confirm no frozen knife block was altered. Did the fold
introduce any new ambiguity, inconsistency (e.g. clause 1 vs the coverage note vs the CHANGE-LOG), or
budget/output regression? Is the recorded Opus-Q3 sampler caveat (a DEAD control's digest can move under a
singularity-changing knife → possible false positive, never false negative) correctly scoped to the build
review and not a directive defect?

## Physics filter
Report a finding only if it catches a way the ablation cert could be wrong or misleading — an inert/invisible
knife, an uncovered object, a witness blind-spot, a residual output leak, a disturbed frozen knife, or a
self-contradiction a builder could implement wrongly. No style preferences.

## Output
Per-question verdict (SOUND / finding) with directive line + engine evidence, and an overall SOUND / NOT-SOUND
on whether the fold resolves F1/F2/F3, introduces nothing new, and leaves the frozen knife design intact.
