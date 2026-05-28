---
unit_id: 166
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 166

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py:52-56`: a new `banner("General inversion forms (paper Sec. 2)")` block inserted after the four solved-form prints (line 50) and before `banner("Forward verification")`, with four assertions:
  - `expect_zero("drho general", sol[drho] - sp.Rational(1,2)*dTheta)`
  - `expect_zero("da general", sol[da] - (sp.Rational(1,2)*dKs - sp.Rational(1,4)*dTheta))`
  - `expect_zero("dcs general", sol[dcs] - (sp.Rational(1,2)*dKs - sp.Rational(1,4)*dTheta + sp.Rational(1,5)*dP))`
  - `expect_zero("dZ general", sol[dZ] - (dKq - sp.Rational(2,5)*dP))`
- Mathematica `mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:47-51`: matching `expectZero` block (`drho general`, `da general`, `dcs general`, `dZ general`) inserted after the four `Print` lines (45) and before `banner["Forward verification"]`.

**Assessment:**
Correct and matches the directive's "required change" verbatim. The two load-bearing additions (`drho general` and `da general`) retain the full `dTheta` term, so they directly pin the headline §2 slopes `½ dTheta` (drho) and `−¼ dTheta` (da) that were previously only printed. These are non-tautological: they compare the `solve`/`Solve` output against the independently-stated paper §2 boxed forms, not against the equations the solver was fed; a wrong slope would make `sol[drho] - dTheta/2 ≠ 0` and trip the assertion. The `dcs`/`dZ` general checks are the (already-pinned-via-bundle) completeness additions the directive flagged as optional. Exec logs confirm `drho general = 0`, `da general = 0`, `dcs general = 0`, `dZ general = 0` (sympy lines 13-16) and matching `PASS:` lines (mathematica 14,16,18,20). No collateral edits in this block.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
1. Banner fix (both engines): `.py:33` and `.wl:26` changed from `"STAGE 149 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS"` to `"STAGE 166 — …"`. Confirmed in the git diff and in both refreshed transcripts (sympy line 3, mathematica line 3 both read "STAGE 166").
2. Independent matrix-inverse cross-check (mathematica only, `.wl:59-76`): a new `banner["Independent matrix-inverse cross-check"]` block inserted after `Forward verification` and before the full-bundle-form banner. It builds the forward map `Mmat`, computes `inv = Inverse[Mmat]`, `solVec = inv . {dTheta,dKs,dKq,dP}`, asserts each component against the paper §2 forms (`matrix drho/da/dcs/dZ`), and a `matrix round-trip` check `Total[(Mmat . solVec - {dTheta,dKs,dKq,dP})^2]`.

**Assessment:**
Correct. I verified `Mmat`'s rows transcribe eq1..eq4 of the forward map `(drho,da,dcs,dZ) → (dTheta,dKs,dKq,dP)`:
- row1 `{2,0,0,0}` = eq1 `dTheta = 2 drho` ✓
- row2 `{1,2,0,0}` = eq2 `dKs = drho + 2 da` ✓
- row3 `{0,-2,2,1}` = eq3 `dKq = -2 da + 2 dcs + dZ` ✓
- row4 `{0,-5,5,0}` = eq4 `dP = -5 da + 5 dcs` ✓

This is a genuinely distinct derivation path (matrix inverse) from the SymPy/Mathematica `Solve` path, and the four `matrix drho/da/dcs/dZ` checks pin M1–M4 against the paper §2 forms independently. The `matrix round-trip` line uses the orchestrator's documented sum-of-squares scalarization `Total[(…)^2]` (replacing the directive's raw length-4 vector residual, which `expectZero`'s `res === 0` test rejects for a list); per the verifier instructions this is an intentional, documented fix and not a deviation to flag. The round-trip is `0` for any invertible `Mmat` (it pins invertibility/consistency, not the matrix entries — those are pinned by the `matrix drho/da/dcs/dZ` comparisons), so the block as a whole non-tautologically establishes M1–M4 by a second route. Exec log confirms all five lines `= 0` / `PASS` (mathematica 37-46). No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (no FAIL/traceback; all `… = 0` lines present, all expect_zero passed). Notable lines:
- `drho general = 0` / `da general = 0` / `dcs general = 0` / `dZ general = 0` (the F1 additions, lines 13-16)
- `STAGE 166 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS` (line 3, F2 banner fix)
- `rho_w^(chi) = 0.403417022451042341 * lambda_mu^(-1)` (line 49, matches paper numeric)

**Mathematica:** exit=0. Notable lines:
- `PASS: drho general` / `PASS: da general` (lines 14-15, F1)
- `PASS: matrix drho` … `PASS: matrix round-trip` (lines 38-46, F2 independent cross-check)
- `Stage 166 Mathematica audit passed.` (line 86; script ends `Exit[0]`)
- `STAGE 166 — …` (line 3, F2 banner fix)

19 expectZero checks all PASS; consistent with the orchestrator note ("Re-run: 19 PASS, 0 FAIL, exit 0").

**Output freshness:** confirmed re-generated post-fix.
- SymPy: script mtime 2026-05-28 15:54:48 < output mtime 16:10:13.
- Mathematica: script mtime 2026-05-28 16:12:20 < output mtime 16:13:24.
Both outputs are newer than their scripts.

## Material-change assessment

`material_change`: false.

No derived result changed. F1 added assertions over the *already-computed* `sol` (no new quantity, no altered formula); the solved forms `drho=½dTheta`, `da=½dKs−¼dTheta`, `dcs=½dKs−¼dTheta+⅕dP`, `dZ=dKq−⅖dP` and the numeric `ρ_w^(χ)=0.403417…` are byte-for-byte unchanged from the pre-fix transcript. F2 changed a banner string and added a redundant second-route confirmation of the same coefficients. Nothing downstream units consume is altered.

## Side observations (non-blocking)

- The `matrix round-trip` check `Total[(Mmat . solVec - obs)^2]` is identically zero for any invertible `Mmat` (it composes `M · M⁻¹`), so it verifies internal consistency/invertibility rather than the correctness of `Mmat`'s entries. The entries themselves are correctly pinned by the four `matrix drho/da/dcs/dZ` comparisons against the §2 forms, so the block's coverage is sound. Not a defect — noted only for completeness.
- The original auditor labelled A1–A4 (forward re-substitution) "near-tautological"; these remain in both scripts unchanged. They were not the subject of either finding's required change, and F1's new general-form assertions now provide the non-tautological coverage that motivated the finding. No action needed.

## Verdict justification

Both findings are fully resolved. F1's required general-form assertions are present in both engines verbatim, retain the `dTheta` slopes they are meant to pin, and pass non-tautologically (comparing solver output to the independently-stated paper §2 forms). F2's banner fix landed in both files and the Mathematica script now carries a genuinely independent matrix-inverse derivation of M1–M4 plus a transcribed-correctly round-trip; I verified `Mmat` against eq1..eq4. The only deviation from the directive (round-trip scalarization to `Total[(…)^2]`) is the orchestrator's documented, intentional fix. Both exec logs are fresh, show all assertions passing, and exit 0. No regressions in the diff, no collateral edits. Verdict: verified.
