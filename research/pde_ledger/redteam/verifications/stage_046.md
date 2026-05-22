---
unit_id: 046
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 046

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl:55-129` was rewritten. The hand-typed `pR`, `p1`, `p2`, `dGExpected`, `dFExpected`, `gDiffExpected`, `fDiffExpected` literals (formerly lines 55-89) are removed entirely (confirmed by reading the file and by `stage_046_diff.patch` lines 9-43 showing deletion). The replacement (wl:55-106) derives `dGdR = Together[D[gTr, r]]`, `dFdR = Together[D[fTr, r]]`, `deltaG = Together[gTr - gFlat]`, `deltaF = Together[fFlat - fTr]` directly from `gTr`/`fTr`, then uses `Reduce[ForAll[...] dGdR < 0]`, etc., to confirm sign claims on the open domain, and `PolynomialQuotientRemainder` to verify `(1 - r^2)` divides the numerator of `deltaG` and `(1 - r)` divides the numerator of `deltaF`.

**Assessment:**
The edit matches the directive verbatim — the diff is the directive block placed at wl:55-106. Grep for `4*r^4*xi^3`, `18*r^2*delta^2*xi`, `18*r^3*delta^2*xi^2`, `dGExpected`, `dFExpected`, `gDiffExpected`, `fDiffExpected` returns only the comment line "no hand-typed p1/p2/gDiffExpected/fDiffExpected" (wl:75) — no remaining literal source. The Mathematica transcript (mathematica_audit.txt:13-32) shows the four `Reduce[...] = True` results and `PASS` lines for `dG/dR < 0`, `dF/dR > 0`, `G_tr > G_flat`, `F_flat > F_tr`, plus zero remainders for both polynomial divisions. The sign checks are non-tautological: `Reduce` actually solves the ForAll-quantified inequality over the rationals — a sign typo in `gTr` or `fTr` would yield a non-True reduction and `fail` would trigger. The two engines now perform genuinely different work: SymPy compares to hand-typed factored expected forms (preserved at sympy:67-141) while Mathematica derives factorisations directly and runs `Reduce`. No collateral edits beyond the directive block.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:143-186` adds the `3b. Sign verification of branch-difference factors` banner with four boundary-value `expect_zero` calls and three interior rational sample points (R=1/4, 1/2, 3/4 with varied xi/delta), with `raise AssertionError` on non-positive samples.
- `mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl:108-129` adds the analogous `3b. Boundary-value sign checks` banner with four `expectZero` boundary checks and the same three sample points via `Do[...]` with `fail` triggers on `!TrueQ[gs > 0]` / `!TrueQ[fs > 0]`.

**Assessment:**
Both insertions match the directive block-for-block. The SymPy transcript (sympy_audit.txt:38-44) shows all four boundary `= 0` lines and three sample lines with the rational values `225/8869`, `38617837960/99381932001`, `81/1736`, `759648230/1473329763`, `91/21935`, `5842146019415/70196178995856` — all positive and identical between the two engines (Mathematica mathematica_audit.txt:37-47), confirming engine cross-agreement at the sampled points. The boundary checks are non-tautological: a sign error in `(1 - R)` inside `F_diff_expected` (the regression scenario the auditor flagged) would not cancel here because these new assertions act on `F_flat - F_tr` and `G_tr - G_flat` constructed directly from `F_tr`, `G_tr`, `F_flat`, `G_flat`, not from the hand-typed `*_diff_expected` polynomials. The numerical sample assertions use strict `<= 0` triggers — a regression yielding negative `delta_G` or `delta_F` would raise; SymPy passing requires all three samples positive.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- sympy_audit.txt:13-14: `strong-split endpoint for G = 0`, `strong-split endpoint for F = 0`
- sympy_audit.txt:38-41: `G_tr - G_flat vanishes at R=1 = 0`, `F_flat - F_tr vanishes at R=1 = 0`, `G_tr at R=0 equals xi = 0`, `F_tr at R=0 equals 1/(1-xi) = 0`
- sympy_audit.txt:42-44: three sample lines with strictly positive rational values for both `G_tr - G_flat` and `F_flat - F_tr`
- sympy_audit.txt:53: `All Stage-29 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- mathematica_audit.txt:15-18: `Reduce[dG/dR < 0 on (0,1)^3] = True` and `Reduce[dF/dR > 0 on (0,1)^3] = True` with corresponding `PASS` lines
- mathematica_audit.txt:23, 27: `PASS: (1 - r^2) divides numerator of G_tr - G_flat`, `PASS: (1 - r) divides numerator of F_flat - F_tr`
- mathematica_audit.txt:29-32: `Reduce[G_tr - G_flat > 0 on (0,1)^3] = True` and `Reduce[F_flat - F_tr > 0 on (0,1)^3] = True` with PASS lines
- mathematica_audit.txt:37-44: four boundary-value PASS lines and three sample lines with positive values matching the SymPy values exactly
- mathematica_audit.txt:49: `Stage 046 Mathematica audit passed.`

**Output freshness:** Confirmed post-fix. mtimes:
- sympy_audit.py = 1779475951; sympy_audit.txt = 1779476060 (txt newer by 109s)
- mathematica_audit.wl = 1779475951; mathematica_audit.txt = 1779476070 (txt newer by 119s)
Both `.txt` outputs were regenerated after the script edits.

## Material-change assessment

`material_change`: false.

The edits added new verification assertions (F2) and replaced a transliterated check with independent CAS reasoning (F1). The closed-form expressions for `G_tr`, `F_tr`, `G_flat`, `F_flat` themselves are unchanged. No derived constant, no symbolic identity, and no downstream-consumed result is altered — the audit now proves more about the same expressions. Downstream units that depend on the tracking-branch bound claim see strengthened (not modified) evidence.

## Side observations (non-blocking)

- The SymPy script header docstring still says "Stage 29 SymPy audit" (sympy_audit.py:3) and section banner says "STAGE 29" (sympy_audit.py:39); the file is `stage046`. Pre-existing labelling inconsistency unrelated to either finding.
- The new SymPy sample-loop uses `if g_sample <= 0` / `if f_sample <= 0` directly on a SymPy expression rather than coercing to a numeric/Rational comparison. For the rational sample points chosen, SymPy returns concrete `Rational` results so the comparison is well-defined, but this would break if non-rational samples were introduced. Non-blocking.

## Verdict justification

Both findings are fully addressed by edits that match the directive verbatim, with no collateral changes beyond the inserted blocks (diff confirms only `+` lines in the specified regions). The Mathematica file no longer contains any hand-typed `pR`/`p1`/`p2`/`*Expected` literals; `Reduce`-based sign claims and `PolynomialQuotientRemainder`-based factor checks now drive Mathematica's independent verification. The new SymPy/Mathematica boundary and sample assertions are non-tautological because they operate directly on `G_tr - G_flat` / `F_flat - F_tr` (not on the hand-typed `*_diff_expected` polynomials), so a sign typo in those polynomials would now fail these new checks while still passing the legacy `expect_zero(... - *_diff_expected)` lines. Both engines exit 0; outputs are fresh; no regressions appear in the diff.
