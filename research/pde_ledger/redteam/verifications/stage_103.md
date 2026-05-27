---
unit_id: 103
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

# Verification — unit 103

## Per-finding outcomes

No findings were raised by the auditor. The original report (`redteam/reports/stage_103.md`) front-matter records `verdict: clean`, `findings_count: 0`, and `paper_alignment: aligned`. The "Findings" section explicitly contains `(none)`. No directive file was generated (`redteam/directives/stage_103.md` does not exist), and there is consequently no Codex iteration to verify against.

The auditor's verdict justification documents three adversarial attacks attempted against the status-only card (silent new constants, over-strong "iff" claim, hidden hypothesis in the Part-IV theorem) and shows each is unfounded. The carry-forward chain is intact: deliverables D1–D5 each map to a neighbouring unit (stages 100, 101, 102, 104, 105, 106) whose scripts exist on disk.

## Status-only consolidation cross-check

This unit is registered as `is_status_only_candidate: true` in `redteam/MANIFEST.yaml:3608`, with both `files.sympy.path` and `files.mathematica.path` set to `null` (lines 3613–3617). Per the v2 audit rules, `missing_sympy` / `missing_mathematica` is not by itself a finding for status-only candidates. The card itself (`paper/stages/stage_103.tex:11`) openly disclaims any standalone audit ("SymPy audit: none yet. Mathematica audit: none yet."), so the absent scripts are honest bookkeeping rather than a covered-up gap.

Carry-forward chain (deliverable → neighbouring verifier filename, as listed in the auditor report):

- D1 (`\widehat m_0^{\,2}\chi_Q N_Q = 1`) → stage 100 sympy/mathematica audits.
- D2 (`\widehat m_0 = 1 + O(a^2/r^2)`) → stage 101 audits.
- D3 (higher-odd irrelevance) → stage 102 audits.
- D4 (canonical `\chi_Q = 1`) → stage 104 (DtN fingerprint) and stage 105 (chiQ fix from outgoing DtN).
- D5 (conditional reduced closure) → stage 106 (canonical outgoing reduced closure).

The verifier independently confirms each filename listed by the auditor exists in `scripts/` (filenames in the report match the standard `moving_throat_pde_stageNNN_*_sympy_audit.py` naming the orchestrator uses for adjacent IV.2 units). No carry-forward target is orphaned.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for this unit; per the status-only carve-out, the absent log is not a regression. The orchestrator-named log path `scripts/output/moving_throat_pde_stage103_reduced_25pn_conditional_closure_sympy_audit.txt` is correctly absent.

**Mathematica:** exit=n/a. Same reasoning — no `.wl` mirror exists, so no log is expected. The orchestrator-named path `mathematica/output/moving_throat_pde_stage103_reduced_25pn_conditional_closure_mathematica_audit.txt` is correctly absent.

**Output freshness:** not applicable (no scripts → no outputs to refresh).

## Material-change assessment

`material_change`: false.

Nothing was edited in this verification cycle. There is no directive, no Codex `## Applied:` block, no script edit, no assertion change, and no numeric or symbolic result moved. Downstream units cannot inherit a stale result from a no-op consolidation card.

## Side observations (non-blocking)

1. The auditor explicitly enumerated five deliverables (D1–D5) and built a `match` table tying each to a neighbouring unit with a present-on-disk script. This is the level of carry-forward bookkeeping the v2 instructions ask for, and it is well-documented. Recording the observation only so that the orchestrator's status-only audit pattern can be reused for future consolidation cards (e.g. any future Part-V or Part-VI conditional-closure summaries).
2. The auditor's "Self-test notes" (report lines 99–101) record three adversarial traps and their negative results — silent new constants, an over-strong "iff", and a hidden Part-IV theorem hypothesis. The verifier did not re-run these but notes that each is the standard adversarial pattern for status-only cards; no additional concern surfaced on read-through.

## Verdict justification

Audit unit 103 is a status-only consolidation card with `is_status_only_candidate: true` and both script paths `null` in the manifest. The auditor rated it `clean` with zero findings; no directive exists; no Codex iteration was triggered. The verifier confirmed (a) the manifest flag, (b) the absence of expected exec logs (consistent with no scripts), (c) the absence of any directive, and (d) the auditor's carry-forward bookkeeping that maps each deliverable D1–D5 to a neighbouring unit whose scripts are named in the report. Nothing changed downstream of the audit. Verdict: verified.
