---
unit_id: 169
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T17:30:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 169

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
In `mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl`, the "Stage 168 weighted transport" block: the transcribed decimal radius `rNum = SetPrecision[1.77799353547498, 30]` (old L101) is replaced by the canonical closed form `rExact = Sqrt[4107 - 100*Pi^2]/(10*Pi)` followed by `rNum = N[rExact, 30]` (new L101-102). The diff is exactly two-for-one line replacement; nothing else touched.

**Assessment:**
The change matches the directive's prescribed before/after verbatim. Verified by grep that the bare decimal `1.77799353547498` no longer appears anywhere in the `.wl` — the radius is now computed, not transcribed, so the Mathematica numeric evaluation of the three Family-1 coefficients is independent of the SymPy decimal. The three paper-comparison TARGET literals (`0.758035078944663`, `1.00314310113848`, `1.88373219118005`) remain hand-typed in the `checks` table at L111-113 — correctly preserved, since those are the paper-side values compared AGAINST. `gNum = SetPrecision[0.758035078944663, 30]` (L103) legitimately remains a decimal: no closed form available for `g_*`, as the directive notes. The sphere-average (integral-based `sphereAvg`) and matrix-invariant (`iXY = (dX . gMat . dY)/5`) sections were already structural/independent and were correctly left untouched. No collateral edits. The coefficient checks remain non-tautological: `coeffT/V/L` are extracted via `Coefficient[...]` from the symbolically-built `xiPerp` evaluated at the now-derived radius, then compared against the independent paper literals at the `1e-12` tolerance.

## Exec log assessment

**SymPy:** exit=n/a. The directive touched only the `.wl`; the SymPy script and output are unchanged (no sympy log was captured, expected). SymPy output `.txt` mtime is 2026-05-28, consistent with no change.

**Mathematica:** exit=0. Notable lines:
- `PASS: eps_perp - Xi_perp Igrp`
- `Xi_perp coeff on xiT = 0.7580...995391011`20.` (paper 0.758035078944663) → `PASS`
- `Xi_perp coeff on xiv = 1.0031...071024522529`20.` (paper 1.00314310113848) → `PASS`
- `Xi_perp coeff on xiL = 1.8837...212034479298`20.` (paper 1.88373219118005) → `PASS`
- `Stage 169 Mathematica audit passed.`
All eight checks PASS; coefficients agree with paper values to ~15 digits.

**Output freshness:** Mathematica `.txt` mtime 2026-06-08 17:12:52 is newer than the `.wl` mtime 2026-06-08 16:32:46 — output was re-generated post-fix.

## Material-change assessment

`material_change`: false. The edit only de-transcribes the radius: `N[Sqrt[4107-100*Pi^2]/(10*Pi), 30]` is numerically identical to the former hardcoded `1.77799353547498` to full precision, and the exec log confirms the three coefficients are unchanged (still pass the `1e-12` checks). No derived RESULT value changed; downstream units are unaffected by the substance, only the second-engine independence improved.

## Side observations (non-blocking)

None.

## Verdict justification

The single finding (F1, mathematica_transliteration) is fully resolved: the load-bearing Family-1 radius is now derived from its canonical closed form rather than copied from the SymPy decimal, the paper-comparison targets are correctly preserved as hand-typed values, `gNum` legitimately stays decimal, the script exits 0 with all three coefficient checks and the structural checks passing, and the output `.txt` was refreshed post-fix. The change is numerically identical (de-transcription only), so `material_change: false`. Verdict: verified.
