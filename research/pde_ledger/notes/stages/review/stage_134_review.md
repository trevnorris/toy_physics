# Review: Stage 134 — First order selected correction

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage134_first_order_selected_correction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage134_first_order_selected_correction_sympy_audit.py`

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
**Notes Derivation Review:** First-order source: Sigma_*(1-Rtilde_*) with Rtilde_* = R_*-<R_*>_* correct. Moment shifts delta g = -Cov_*(c,R_*), delta S = -Cov_*(K_q,R_*) correct. Bias retuning delta Pi = Cov_*(c,R_*)/g'_*. Traction shift delta T = A_T delta g + B_T delta S from Stage 130. R_*(x) from Stage 133 is sole remaining unknown.
**Script Review:** Covariance identities for delta_g, delta_S checked via expect_zero. deltaPi and deltaT printed but not asserted (minor coverage gap). All pass (exit code 0).
**Issues Found:** None (minor: deltaPi/deltaT not asserted, only printed).

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The centered-residual form of the first-order selected correction is correct, and the moment shifts are the expected covariance projections of `R_*(x)`.
2. The bias and traction retuning formulas are consistent with the linearized expansion and with the Stage 130 rigidity coefficients.
3. The stage correctly identifies the remaining unknown as the actual numerical covariances of `R_*(x)`, not a new branch-selection ambiguity.

**Script Review:**

The script verifies the covariance identities directly. The printed `deltaPi` and `deltaT` forms are consistent with the script’s sign conventions and the note’s covariance interpretation.

**Issues Found:** None.

---
