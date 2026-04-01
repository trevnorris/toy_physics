# Review: Stage 170 — Orbit quotient closure

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage170_orbit_quotient_closure.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage170_orbit_quotient_closure_sympy_audit.py`

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
**Notes Derivation Review:** Finite upgrade of Stage 169 tangent-space result. Key insight: monomial invariants make fibre condition exactly linear in log-ratios, so same rank-3 matrix M_* governs both infinitesimal and finite problems. All three rows of M_* verified against monomial exponents. 3×3 minor determinant = 1+chi_{0,*} > 0 confirmed. Solution formulas for (Delta_eta, Delta_T, Delta_mu) match Stage 169 orbit laws identically. Quotient M_+/G_* = (R_{>0})^3 correct (free action on positive states). Observable formulas from Stage 166 correctly cited. Double-equivalence chain for final theorem is clean. Claims properly qualified — remaining open question (branch selection) explicitly identified as outside scope.
**Script Review:** M_* matrix and row equations cross-checked. Minor determinant verified. sp.solve compared against hand-derived formulas. Back-substitution into all three equations. No hardcoded values, no tautologies. All pass (exit code 0).
**Issues Found:** None. Clean capstone stage.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The quotient-closure argument is a clean finite version of Stage 169: the fibre condition is exactly the similarity orbit, and the quotient by that orbit is `(R_{>0})^3`. I independently checked the finite log-ratio equations against the orbit laws from the previous stage.

**Script Review:** The audit script reproduces the matrix form, solves the finite equations, and back-substitutes the result into all three constraints. The saved output is clean, with no missing cases or branch issues.

**Issues Found:** None.

---
