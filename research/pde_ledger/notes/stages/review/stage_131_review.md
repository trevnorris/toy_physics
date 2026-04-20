# Review: Stage 131 — Representative positive families

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage131_representative_positive_families.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage131_representative_positive_families_sympy_audit.py`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Script runs without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** MINOR
**Notes Derivation Review:** Uniform family: g_u = 2/pi, S_u = 2 tanh(pi/2)/pi verified. Derivative family: g_d = pi/4, S_d correct. Retuning formulas from Stages 129-130. Convex interpolation affine in lambda. Zero-crossing lambda_(Pi,0) ~ 0.816, lambda_(T,0) ~ 0.813. Stage 109 consistency 1-lambda ~ xi_* confirmed to ~1.7e-15. Broadening raises canonical point, sharpening lowers it.
**Script Review:** Moments from exact formulas. A_T, B_T from Stage 130. Interpolation solved symbolically. All genuine. All pass (exit code 0).
**Issues Found:**
1. **(MINOR)** dPi_d/eps displayed as ...921 in notes but script gives ...928 (7 ULP). Several other last-digit discrepancies.
2. **(MINOR)** "Stage-110 broadening fraction" should be "Stage-109" — xi_* originates in Stage 109.
3. **(MINOR)** xi_star hardcoded as 15-digit string rather than exact algebraic formula.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The representative uniform and self-matched derivative families reproduce the stated first-order bias and traction shifts, and the convex interpolation behaves affinely as claimed.
2. The zero-crossing interpolation points are in the expected interval and agree with the earlier exact compensation broadening fraction to numerical precision.
3. The stage remains mathematically sound, but the last-digit mismatches and the Stage-110/109 cross-reference typo are real enough to keep the verdict at `MINOR`.

**Script Review:**

The audit script evaluates the exact moments, computes the retuning coefficients, and solves the interpolation zero points. I independently rechecked the interpolation coefficients; the formulas are right, but the notes still disagree with the script in the last digits and the stage reference.

**Issues Found:**

1. The `dPi_d/eps` value in the notes still differs from the script output at the last digits.
2. The Stage-110 cross-reference should still be Stage-109.
3. `xi_*` is still hardcoded rather than given as an exact expression.

---
