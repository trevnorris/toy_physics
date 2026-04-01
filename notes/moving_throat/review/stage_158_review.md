# Review: Stage 158 — Wall normalized load shape

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage158_wall_normalized_load_shape.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage158_wall_normalized_load_shape_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** Wall-normalized variables with exact homogeneity degrees: Delta~K^2, Q~K^3, P~K^2. Bundle factorizations Z0=K*Upsilon, N0=Lambda^2. Defect rewrites Sigma_Z=delta ln Upsilon, Sigma_N=delta ln(Lambda^2/K). Conservative-shape theorem. No-go for common self-similarity (Xi_load=-delta_K). Outgoing-load theorem: 2 delta ln Lambda_r = delta_K. All verified.
**Script Review:** B0 factorization, homogeneity laws (genuine polynomial expansion), bundle factorizations, differential identities via dlog helper (chain rule through composite expressions). 13 expect_zero checks. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The homogeneity bookkeeping is correct. `B0`, `Delta`, `Q`, `P`, `Z0`, and `N0` factor with the stated powers of `K`, the defect identities reduce to logarithmic derivatives exactly, and the conservative-shape and outgoing-load corollaries follow cleanly. The no-go statement for common self-similarity also matches the algebra.

**Script Review:**

The audit script checks the factorization laws, the composite logarithmic identities, and the downstream corollaries. The saved output shows every expected zero check succeeding, and I do not see any hidden assumption that would weaken the result.

**Issues Found:**

None.

---
