---
unit_id: 154
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 154

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl`:
- Line 26: banner now reads `"STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP"`.
- Lines 37–39: `rShiftSeries = Normal[Series[rFun /. g -> gStar + dg, {dg, 0, 2}]]`; `rShiftExpected` is now the paper-side closed form only.
- Lines 47–49: manual substitution dictionary replaced with single-eps parameterization `piExprEps = piExpr /. {dSigma0 -> epsLin*dSigma0, dR -> epsLin*dR, dS -> epsLin*dS}` and `piLin = Expand[Normal[Series[piExprEps, {epsLin, 0, 1}]] /. epsLin -> 1]`.

**Assessment:**
The directive asked for `Series`-driven derivation; the orchestrator's single-eps parameterization is the correct Mathematica idiom for joint linear truncation (multivariate `Series` would have retained cross-terms, as the orchestrator notes). The `expectZero` targets `rShiftExpected` and `dPiExpected` remain the paper claims, so the checks now compare engine-derived series against paper-side identities rather than mirroring a SymPy substitution table. No forbidden substitution rules (`dSigma0*dR -> 0`, etc.) remain. Two `Series[...]` calls present. Banner correct in saved output (line 3). Non-tautological: a sign error in either expected form would still fail.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `R(g_star) - 1/4 = 0`
- `exact shifted R formula = 0`
- `dPi identity = 0`

**Mathematica:** exit=0. Notable lines:
- `STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP`
- `PASS: R(g_star) - 1/4`
- `PASS: exact shifted R formula`
- `PASS: dPi identity`

**Output freshness:** mathematica .txt mtime 1779989691 > .wl mtime 1779989585; sympy .txt 1779989421 > .py 1779945093. Both regenerated post-fix.

## Material-change assessment

`material_change`: false. No derived results changed; only the Mathematica derivation path. The `expectZero` targets are identical to the pre-fix versions, and both engines still print residual `0`. No downstream propagation.

## Side observations (non-blocking)

None.

## Verdict justification

F1 is fully resolved: the .wl no longer mirrors the .py — `rShift` and `piLin` are produced by Mathematica's `Series` machinery (the single-eps trick being the correct workaround for joint linear truncation), the forbidden substitution rules are gone, the banner is fixed, both `expectZero` checks PASS with non-tautological targets, and the saved outputs are fresh. Engine cross-check is now meaningful.
