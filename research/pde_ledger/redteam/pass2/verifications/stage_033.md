---
unit_id: 033
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 033

## Per-finding outcomes

### F1 — stale_output (residual self-labels + stale transcripts)

**Classification:** resolved

**What changed:**
The diff (`exec_logs/stage_033_diff.patch`) contains exactly three label-string edits, nothing else:
- `scripts/...stage033...sympy_audit.py:3` — docstring `Moving-throat PDE Stage 16 SymPy audit.` → `Stage 33`.
- `scripts/...stage033...sympy_audit.py:133` — `print("All Stage 16 checks passed.")` → `Stage 33`.
- `mathematica/...stage033...audit.wl:134` — comment `the Stage 16.1 monotonicity identity and the Stage 16.6 gate identity at` → `Stage 33.1` / `Stage 33.6` (base number only; `.1`/`.6` sub-stage suffixes preserved).

The two committed `.txt` transcripts were also regenerated (orchestrator re-run).

**Assessment:**
The edit matches the directive's "required change" exactly and is purely label/comment-string text. ZERO change to any equation, value, variable, assertion, or `expect_zero`/`expectZero` target — confirmed by the unified diff: the only non-context lines are the three string substitutions, all single-token `16`→`33`. No collateral edit. A grep for residual `Stage 16` / `STAGE 016` / `STAGE 16` / `16.1` / `16.6` across both scripts and both refreshed outputs returns empty. The finding (cosmetic banner-label drift) is fully closed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 33.1 — EXACT MICROSCOPIC NORMALIZATION PRODUCT` (canonical banner)
- `dN/dalpha - monotonicity formula = 0`; `alpha_crit - closed finite-throat form = 0`; `N_-(0) - beta0 kappa0^2 / A = 0`; both weak-loading coefficient residuals `= 0`; `K0_onset - [...] = 0`
- `denominator ratio (must be parameter-free) = 9*pi**2` (load-bearing gate check intact)
- closing `All Stage 33 checks passed.`

**Mathematica:** exit=0. Notable lines:
- `STAGE 033 — MICROSCOPIC NORMALIZATION EQUATION` (canonical banner)
- `PASS:` on all checks: monotonicity, alpha_crit, N_-(0), both weak-loading coefficients, `N_-(0) at K0_onset - NQ`, `K0_onset - [...]`, gate denominator parameter-free (`denominator ratio = -9*Pi^2`, `PASS`), and the two numeric monotonicity residuals (`0``78.83...` / `0``78.99...`, high-precision zero, `PASS`).
- closing `Stage 033 Mathematica audit passed.`
- Pre-existing `N::meprec` precision-limit notices appear in the two numeric monotonicity sub-checks; these are unchanged high-precision-arithmetic warnings (the residual still resolves to high-precision 0 with `PASS`) and are unrelated to the label-only edit.

A grep for `fail`/`error`/`traceback` across both refreshed transcripts returns empty.

**Output freshness:** confirmed. Both `.txt` outputs (mtime 1780636408) are newer than their corresponding script sources (`.py` and `.wl`, both mtime 1780635821), so the transcripts were regenerated post-fix.

## Material-change assessment

`material_change`: false.

The edit touched only docstring text, one print-summary string, and one `.wl` comment. No derived result, closed form, or assertion target changed. Every value-reconciliation entry from the auditor's report still holds (the refreshed transcripts reproduce the same closed forms: `alpha_crit`, `N_-(0)=8 beta0/(pi^2 A)`, weak-loading coef `64 beta0(8A+9DK)/(9 pi^4 A^2 DK)`, `K0_onset`, gate `den_ratio = ±9 pi^2`). No downstream unit can depend on a banner string.

## Side observations (non-blocking)

The `N::meprec: $MaxExtraPrecision = 50 reached` notices in the Mathematica monotonicity numeric cross-check are informational and pre-date this fix; the check still resolves to high-precision zero and reports `PASS`. Not part of any finding; flagged only for awareness.

## Verdict justification

The sole finding (F1, low-severity cosmetic stale_output) is fully resolved. The diff is purely the three prescribed label/comment-string substitutions with zero math/equation/assertion/value change; both engines re-run clean to exit 0 with canonical `STAGE 33.x` / `STAGE 033` banners and the closing `All Stage 33 checks passed.` / `Stage 033 Mathematica audit passed.`; every prior `PASS` / `= 0` residual is preserved; no FAIL; refreshed outputs are newer than their scripts; no residual `16` labels remain. No regressions, no new findings. Verdict: verified, material_change: false.
