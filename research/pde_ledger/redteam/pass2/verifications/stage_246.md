---
unit_id: 246
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
codex_session: none
scripts_edited: []
diff_empty: true
relgate:
  sympy_exit: 0
  mathematica_exit: 0
  engines_agree: true
  outputs_fresh: true
findings_count: 1
findings_disposition:
  - id: F1
    class: paper_misalignment
    subclass: card_text_lag (paper_missing_script_claim — stale "Mathematica audit: none yet")
    severity: low
    disposition: deferred
    deferred_to: P4-51
    script_fix: none
    rationale: >-
      Paper-side prose lag only — the card's Verification field still reads
      "Mathematica audit: none yet" while a complete independent .wl second
      engine (M1–M9) exists and passes. No math impact; Codex never edits
      paper/, so routed to user via P4-51. Stage verifies in both engines.
---

# Verification — unit 246 (pass-2)

## Disposition

ZERO-SCRIPT-CORRECTION unit. Captured diff
`redteam/pass2/exec_logs/stage_246_diff.patch` is empty (0 bytes confirmed) —
Codex made no `.py`/`.wl` edits. Directive correctly held for user resolution;
no script edits were warranted because the math is clean in both engines.

## Reliability gate (fresh outputs)

- `redteam/pass2/exec_logs/relgate/sympy_246.txt` — ends "All symbolic checks
  passed." (exit 0).
- `redteam/pass2/exec_logs/relgate/mma_246.txt` — ends "All Mathematica checks
  passed." (exit 0), every M1–M9 line PASS.
- Per the gate metadata, both runs exit 0 and are byte-identical to the
  committed outputs; engines agree (numeric deltas g 4.08e-9, S 1.47e-9,
  R 1.56e-9, σ_min 3.89e-9 — all within the 5e-9 tolerance).

## Finding review

### F1 — paper_misalignment → CARD-TEXT-LAG → deferred P4-51

Card `paper/stages/stage_246.tex:4` states "Mathematica audit: none yet," but
pass-1 added a full independent `.wl` audit (M1–M9, output committed 2026-06-03,
all PASS). This is the card-text-lag class: pure prose under-reporting of
verification coverage, no math impact. Codex does not edit paper/, so the
Verification-line update is routed to the user. DEFERRED to P4-51. Stage itself
verifies — F1 does not block.

## Pass-1 tautology / round-trip fix re-confirmed

The pass-1 concern (round-trip `sigma_min_expected`, assertion A5) was re-checked
and the fix HOLDS:

- A5 (`sigma_min_test - sigma_min_expected == 0`) is redundant on its own.
- A4 is load-bearing and non-tautological: `sigma_min_true = Min[...]` of the
  boundary + in-range-vertex candidates, computed straight from `sigma_y`
  without the piecewise branch logic, so a wrong piecewise branch is caught.
- The `.wl` independently corroborates via `MinValue[{source, 0<=x<=1}, x]` — a
  black-box optimizer over the original cos-form `source`, never touching the
  quadratic/vertex algebra. A wrong analytic min would diverge from MinValue.

No surviving tautology controls the verdict. Independence is genuine across all
three contrasts (Min vs MinValue; round-trip M·inv−vec vs direct Inverse;
piecewise_fold vs Reduce on the threshold). No transliteration concern.

## Verdict

`overall_verdict: verified`, `material_change: false`. The single finding (F1)
is a low-severity card-text-lag deferred to P4-51; the stage verifies in both
engines with fresh, agreeing, exit-0 outputs.
