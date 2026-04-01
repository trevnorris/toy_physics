# Review: Stage 135 — Family1 actual correction

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage135_family1_actual_correction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage135_family1_actual_correction_sympy_audit.py`

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
**Verdict:** PASS
**Notes Derivation Review:** Covariances Cov_*(c,R_*) ~ 0.0648, Cov_*(K_q,R_*) ~ 0.0389. Moment shifts delta g ~ -0.0648, delta S ~ -0.0389. Retuned Pi_corr ~ 2.416, T_corr ~ 1.173. One-step Picard: Pi_1 ~ 2.539, T_1 ~ 1.210 (same direction, stronger). Effective broadening lambda_eff ~ 0.38 from Stage 131 interpolation — Pi and T components agree to 3 sig figs (nontrivial cross-check).
**Script Review:** R_*(x) from Stage 133 formula. Covariances by mpmath numerical integration. All downstream quantities derived (not hardcoded). One-step Picard by exp(-Pi*x-R_*) normalized. Stage 131 constants verified. All 16 outputs match notes. All pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The full-profile correction is consistent with the earlier mouth-layer construction. The residual is tangent-matched at the mouth, the covariances have the stated signs, and the retuned `Pi_corr` and `T_corr` move in the same direction as the one-step nonlinear check, which is the right behavior for a broader source profile.

**Script Review:**

The script computes the residual numerically, evaluates the two covariances, applies the rigidity projection, and checks the one-step nonlinear iterate. The saved output matches the quoted numbers and the direction of the correction is internally consistent.

**Issues Found:**

None.

---
