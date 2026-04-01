# Review: Stage 042 — Operator branch residual bounds

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage042_operator_branch_residual_bounds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage042_operator_branch_residual_bounds_sympy_audit.py`

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

1. **Equation-level correctness.** zeta_phys from substituting Pe_* into Stage 40 family. Bracket zeta_- <= zeta_phys <= zeta_+ from Stage 41 Pe bracket + Stage 40 monotonicity. Residual brackets with correct sign convention. Xi_fail = Pe_req/Delta_inf and Xi_suff = Pe_req/Delta_0 from bracket inversion. Ordering Xi_fail <= Xi_suff from Delta_inf >= Delta_0. Weak-coupling expansion zeta_phys = A_K[1 + ((4-pi)/pi) Xi Delta_0 + O(Xi^2)] correct.

2. **Logical flow.** Clean three-zone structure (fail/intermediate/success) from bracket.

**Script Review:** Abstract Omega function for bracket checks, explicit Omega_Pe for series coefficient. Residual center identity, Xi threshold identities all verified. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. Evaluating the Stage-40 physical ratio on the Stage-41 fixed-point branch is done correctly: the monotonicity in `Pe` transfers the branch bracket into the exact `zeta_- <= zeta_phys <= zeta_+` support window.
2. The residual envelope `R_- <= R_phys <= R_+` and the two theorem gates follow with the correct sign convention. The threshold definitions `Xi_fail = Pe_req / Delta_inf` and `Xi_suff = Pe_req / Delta_0` are the right inversions of the branch bracket.
3. The weak-coupling expansion uses the same Stage-40 small-`Pe` coefficient `((4-pi)/pi)` and is consistent with substituting `Pe_* = Xi Delta_0 + O(Xi^2)`.

**Script Review:**

The script mirrors the note structure cleanly: it checks the abstract branch brackets with `Omega`, verifies the residual-center identity and both threshold identities, and then uses the explicit `Omega_Pe` series for the weak-coupling coefficient. I reran the Stage 42 audit locally and it passed unchanged.

**Issues Found:** None.

---
