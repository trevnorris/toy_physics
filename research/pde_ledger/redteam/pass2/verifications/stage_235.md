---
unit_id: 235
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_patch_empty: true
relgate:
  sympy: pass        # exit 0, ends "All Stage 235 symbolic checks passed."
  mathematica: pass  # exit 0, ends "All Stage 235 Mathematica checks passed."
  byte_identical_to_committed: true
findings_count: 1
findings_resolved: 0
findings_deferred: 1
---

# Verify unit 235 — pass-2 (zero-script-correction)

## Disposition

ZERO-SCRIPT-CORRECTION unit. Codex applied no edits; the captured diff
`redteam/pass2/exec_logs/stage_235_diff.patch` is empty (0 bytes). No `.py`,
`.wl`, `.tex`, or notes file was touched.

## Reliability-gate re-run

Both fresh relgate outputs confirmed:
- `redteam/pass2/exec_logs/relgate/sympy_235.txt` — exit 0, ends
  "All Stage 235 symbolic checks passed." Byte-identical to committed
  `scripts/output/...stage235..._sympy_audit.txt`.
- `redteam/pass2/exec_logs/relgate/mma_235.txt` — exit 0, ends
  "All Stage 235 Mathematica checks passed." Byte-identical to committed
  `mathematica/output/...stage235..._mathematica_audit.txt`.

All M1–M6 / §1–§5 checks PASS; every residual is the zero scalar/vector/matrix.

## Finding disposition

### F1 — paper_misalignment / paper_missing_script_claim — DEFERRED (P4-51)

The report's lone finding is the CARD-TEXT-LAG class: card
`paper/stages/stage_235.tex:11` reads "Mathematica audit: none yet" and notes
§8 is titled "SymPy-backed status," despite a present, clean, fresh-output
`...stage235..._mathematica_audit.wl`. This is a stale verification-coverage
STATUS annotation, not a value mismatch and not a script defect.

Per the standing user decision, this card-text-lag is DEFERRED to
PAPER_CLEANUP P4-51 (paper-side documentation sync; Codex applies any resulting
card/notes edit only under a follow-up directive). The directive itself routes
this to user resolution and applies nothing in-unit. No action taken here.

## Value Reconciliation cross-check

Report's pass-2 reconciliation: 9 deliverable values checked, all MATCH
(M_rm, M_rm^2=I/det=-1, P_nt, P_eta, x_nt, x_eta, codim-two orbit-lock
equivalences, ||x_eta||^2=(1+c_eta^2)q_eta^2, Delta_x_static/Delta_x_orbit).
The sole MISMATCH is the coverage statement, folded into F1 (coverage, not a
value). No value-level discrepancy. Spot-confirmed against the fresh relgate
outputs: SymPy `||x_eta||^2 = qeta_line**2*(c_eta**2 + 1)` and MMA
`x_blind.x_blind = (1 + cEta^2)*qEta^2` agree; both engines' orbit-lock
`Solve` → `{R1->0, E1->0}` unique origin.

## Verdict

overall_verdict: verified — math holds in both engines, outputs fresh and
byte-identical to committed, no script change, sole finding is a deferred
card-text-lag (→ P4-51), not a value or script defect.
