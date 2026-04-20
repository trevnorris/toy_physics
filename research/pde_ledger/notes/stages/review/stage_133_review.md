# Review: Stage 133 — Full profile residual

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage133_full_profile_residual.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage133_full_profile_residual_sympy_audit.py`

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
**Notes Derivation Review:** T_s(x) shell profile satisfies T_s(0)=0, T_s'(0)=1. T_q(x) D/N profile with C_q, A_q coefficients correct. Fixed-point Pi=Sigma_m(4-S_q). Full potential Phi_* = 4 Sigma_m T_s - Sigma_m T_q. Residual R_*(x) = Phi_* - Pi x. Tangent matching: R_*(0)=0, R_*'(0)=0 verified by explicit cancellation (4-S_q-(4-S_q)=0). Curvature R_*''(0) = -3 Sigma_m Pi/(1-e^{-Pi}) verified — factor 3 from shell/mixed ratio 4-1. Negative curvature → source broader than exponential tangent.
**Script Review:** All profiles, derivatives, evaluations. Six expect_zero checks (T_s(0), T_q(0), T_s'(0)-1, R(0), R'(0), R''(0)-target). All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The full-profile shell and mixed channel profiles satisfy the mouth boundary conditions and tangent-matching identities exactly.
2. The curvature residual `R_*''(0)` has the correct negative sign and matches the stated `-3 Sigma_m Pi_*/(1-e^{-Pi_*})` formula.
3. The stage correctly upgrades the problem from branch choice to a finite profile-correction issue around the selected lower compensated branch.

**Script Review:**

The audit script checks the boundary values, first derivatives, and second derivative target directly. I independently verified the curvature target algebra, and it matches the saved output.

**Issues Found:** None.

---
