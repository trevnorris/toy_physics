---
unit_id: 241
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_empty: true
relgate:
  sympy: pass        # exit 0; "All Stage 241 symbolic checks passed."
  mathematica: pass  # exit 0; "Stage 241 Mathematica audit passed."
  byte_identical_to_committed: true
findings_count: 1
findings_resolved: 1
---

# Verification — unit 241 (pass-2, zero-script-correction)

## Summary

Stage 241 ranks the primitive carrier weights along the selected one-parameter
twin-support curve (`eps_* = 1 - 3 varrho/2`, `sigma = 4/(3 varrho) - 2`,
`0 < varrho < 2/3`, `0 < beta < 2/11`), proving `w_chi = 2 w_Z`, deriving the two
phase-diagram thresholds `varrho_WLambda = 2(1+beta^2)/[3(2+beta^2)]` and
`varrho_ULambda = 2(1+beta^2)/[3(1+beta+beta^2)]`, and ordering the three ranking
regions.

This is a ZERO-SCRIPT-CORRECTION unit: Codex made no edits, the captured diff
`redteam/pass2/exec_logs/stage_241_diff.patch` is EMPTY (0 bytes), and the
reliability-gate re-run is byte-identical to the committed outputs.

## Finding disposition

### F1 — paper_misalignment (CARD-TEXT-LAG)
- **Subtype:** `paper_missing_script_claim` — the stage card
  `paper/stages/stage_241.tex:11` still says "Mathematica audit: none yet" while a
  genuinely independent (Solve-based, not a port) `.wl` exists and passes all
  M1–M8 checks.
- **Class:** card-text-lag (paper-side, not a script defect).
- **Disposition:** DEFERRED to P4-51 (published-paper card-text reconciliation
  pass). No Codex/script edit prescribed; the directive holds for user resolution
  (recommended: cite the `.wl` audit on card line 11). Not a blocker for stage
  verdict.

## Notes-only paper_misalignment re-check (HOLDS)
- The pass-1 notes fix ϱ_WΛ `193/369 → 125/369` STILL HOLDS in the notes
  (`...stage241..._sympy_audit.md` line 577, `\frac{125}{369}\approx 0.338753`).
- All four ϱ windows reconcile across `.py`, `.wl` outputs, notes, and appendix:
  `1/3`, `125/369` (~0.338753), `2/3`, `250/441` (~0.566893) — MATCH everywhere.

## Independent-derivation confirmation
- The `.wl` is genuinely INDEPENDENT, not a transliteration: it SOLVES the
  physical crossover equations `wW = wLambda` / `wUmag = wLambda` for the nonzero
  `eps_*` root and re-derives `varrho` via the support law, reaching the same
  closed-form thresholds. The `.py` instead asserts the threshold identities
  directly. Different route, same result = dual-engine standard met.
- Variable-independence self-test: both `d/dbeta` checks are non-vacuous (both
  thresholds depend on `beta`); region positivity uses three distinct in-window
  `varrho` samples. No tautology / absent-variable trap.

## Relgate sanity
- `exec_logs/relgate/sympy_241.txt`: exit 0, ends "All Stage 241 symbolic checks passed."
- `exec_logs/relgate/mma_241.txt`: exit 0, ends "Stage 241 Mathematica audit passed."
- Both byte-identical to committed outputs; diff patch empty.

## Verdict

`overall_verdict: verified`, `material_change: false`. Every mathematical
deliverable matches across both engines; the lone finding is the stale card line
(card-text-lag → deferred P4-51), and the 125/369 notes fix holds.
