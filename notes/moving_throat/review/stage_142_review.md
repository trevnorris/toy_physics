# Review: Stage 142 — Hybrid outlet projection

**Batch:** 18 — Linear Defect Transport & Final
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage142_hybrid_outlet_projection.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage142_hybrid_outlet_projection_sympy_audit.py`

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
**Notes Derivation Review:** delta Sigma_0 cancels in delta C (exact). delta C = -16 sigma_* delta R. 2×2 det = -27s forces delta C = delta kappa_W = 0 → delta g = 0. Final Delta_Q = -9 sigma_* delta gamma_W/(1-sigma_*). All linearizations verified.
**Script Review:** Hybrid outlet linearization, mouth-gain transport, canonical-even solve. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The hybrid-outlet linearization is correct. Expanding around the compensated canonical branch gives the stated formulas for `delta E2`, `delta E4`, and `Delta_Q` in terms of `delta C`, `delta kappa_W`, and `delta gamma_W`.
2. The key cancellation `delta C = delta rho_R - 4 delta sigma_W` correctly removes the direct `delta Sigma_0` contribution from the mouth-gain transport, leaving only the transported ratio defect.
3. Solving the canonical-even preservation conditions forces `delta C = delta kappa_W = 0`, so the last linear 2.5PN defect really does collapse to pure `delta gamma_W`.

**Script Review:**

The script faithfully linearizes the hybrid outlet, checks the mouth-gain transport cancellation, solves the canonical-even constraints, and verifies the final reduced defect formula. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
