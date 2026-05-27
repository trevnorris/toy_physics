---
unit_id: 113
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 113

## Per-finding outcomes

No findings were raised by the auditor. The original report (`redteam/reports/stage_113.md`) front-matter records `verdict: clean`, `findings_count: 0`, and `paper_alignment: aligned`. The "Findings" section explicitly reads `(None.)`. No directive file was generated (`redteam/directives/stage_113.md` does not exist), and there is consequently no Codex iteration to verify against.

The auditor's verdict justification documents three adversarial attacks against this status-only card — (a) a check that the card or notes do not assert a claim that should have an in-unit assertion, (b) a check that no numeric constant in the appendix's `Robin and mixed-pole tests` / `Compensated Robin–mixed outlet` subsections (lines 399–456) is silently relied on without provenance, and (c) a check that the card's `Downstream use` language does not overpromise the compensated branch as an unconditional actual-branch theorem. All three attacks failed, consistent with the `clean` verdict.

## Status-only consolidation cross-check

This unit is registered as `is_status_only_candidate: true` in `redteam/MANIFEST.yaml:3894`, with `is_checkpoint: false` (line 3893) and both `files.sympy.path` / `files.mathematica.path` set to `null` (lines 3899–3904) with `exists: false`. Per the v2 audit rules, `missing_sympy` / `missing_mathematica` is not by itself a finding for status-only candidates. The card's `\stagefield{Verification}` line (`paper/stages/stage_113.tex`) openly disclaims any standalone audit ("SymPy audit: none yet. Mathematica audit: none yet."), so the absent scripts are honest bookkeeping rather than a covered-up gap.

The auditor's paper↔script cross-check table records four carry-forward deliverables relying on the upstream 107–112 outlet block:

- `\chi_Q^R = 3/(3-\rho_R)` (pure Robin) — upstream in 107–112 block.
- Standalone mixed pole forces `\sigma_W = 0` (no-go) — upstream in 107–112 block.
- Compensated branch `(\rho_R, \kappa_W) = (4\sigma_W, 1/3)`, `\chi_Q^{hyb} = (1 - 9\sigma_W\gamma_W)/(1 - \sigma_W)`, canonical iff `\gamma_W = 1/9` — upstream in 107–112 block.
- Audit-path summary "pure scale/argument or compensated Robin–mixed are the only viable routes" — consolidation claim of this card.

Each row is marked `n/a (status carry-forward)`, consistent with the carve-out. The auditor correctly declines to file a cross-unit coverage finding from within this unit's reading scope and notes that such a sweep belongs at the batch-tracker level.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for this unit; the orchestrator-named log path `scripts/output/moving_throat_pde_stage113_outlet_model_status_sympy_audit.txt` is correctly absent. Per the status-only carve-out, the absent log is not a regression.

**Mathematica:** exit=n/a. No `.wl` mirror exists; the orchestrator-named log path `mathematica/output/moving_throat_pde_stage113_outlet_model_status_mathematica_audit.txt` is correctly absent.

**Output freshness:** not applicable (no scripts → no outputs to refresh).

## Material-change assessment

`material_change`: false.

Nothing was edited in this verification cycle. There is no directive, no Codex `## Applied:` block, no script edit, no assertion change, and no numeric or symbolic result moved. Downstream units cannot inherit a stale result from a no-op consolidation card.

## Side observations (non-blocking)

1. The auditor's "Self-test notes" enumerate four adversarial traps (variable independence, path-specification, paper round-trip, carve-out gate) and find each negative. The paper round-trip in particular re-checked appendix lines 402–456 against the notes and confirmed equation-level agreement on the Robin, mixed-pole, and compensated identities. This is the standard adversarial pattern for status-only consolidation cards.
2. The auditor explicitly flags that a cross-unit coverage sweep (whether 107–112 actually verifies the compensated branch identities `(\rho_R, \kappa_W) = (4\sigma_W, 1/3)`, `\chi_Q^{hyb}`, and `\gamma_W = 1/9`) would be a batch-tracker concern rather than a stage-113 finding. Recorded so the orchestrator can decide whether to surface a batch-level coverage item for IV.2; this does not block verification of unit 113.

## Verdict justification

Audit unit 113 is a status-only consolidation card with `is_status_only_candidate: true` and both script paths `null` in the manifest. The auditor rated it `clean` with zero findings; no directive exists; no Codex iteration was triggered. The verifier confirmed (a) the manifest flag, (b) the absence of expected exec logs (consistent with no scripts), (c) the absence of any directive file, and (d) the auditor's carry-forward bookkeeping that defers the three outlet identities and the standalone-mixed no-go to the 107–112 upstream block while the card itself transparently disclaims an in-unit audit. Nothing changed downstream of the audit. Verdict: verified.
