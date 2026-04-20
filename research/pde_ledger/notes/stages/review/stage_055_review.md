# Review: Stage 055 — Explicit branch thresholds

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage055_explicit_branch_thresholds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage055_explicit_branch_thresholds_sympy_audit.py`

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

**Notes Derivation Review:**

1. **Equation-level correctness.** Threshold surfaces Upsilon_fail and Upsilon_suff correctly substituted from Stage-54/49. Two asymptotic regimes verified by hand: shell-gradient dominated (kappa ~ (4/5)Lambda_ell^2) and compression dominated (kappa ~ 4 chi_s^2). Large-alpha approximations (cosh ~ sinh ~ exp/2) correctly applied. All coefficients verified.

2. **Logical flow.** Clean: exact surfaces → physical amplitude thresholds → asymptotic regimes.

3. **Physical interpretation.** Compression-dominated branches face steeper thresholds. Sound.

**Script Review:** Delta_0, Delta_inf from kernel formulas. Threshold surfaces checked. Both asymptotic regimes verified via large-alpha approximations with expect_zero. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The explicit branch substitution is correct: `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`, `eta = Lambda_ell`, and the wall-loading thresholds `Upsilon_fail` and `Upsilon_suff` are the direct `Stage 49 / 52` surfaces rewritten on the explicit branch.
2. The physical amplitude thresholds `V0_fail^2` and `V0_suff^2` are consistent rewrites of the same bounds, and the two asymptotic regimes match the large-`alpha` reductions stated in the notes.
3. The branch interpretation is sound: shell-gradient dominated branches are easier to cross than compression-dominated ones, with the compression regime producing the steeper sufficiency scaling.

**Script Review:**

The audit script faithfully computes `Delta_0`, `Delta_inf`, the explicit `Upsilon` thresholds, the `V0^2` forms, and both asymptotic limits. The saved output agrees with the notes, and the independent SymPy check reproduces the published coefficients.

**Issues Found:** None.

---
