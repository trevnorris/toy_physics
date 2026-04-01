# Review: Stage 137 — Coevolving core mouth map

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage137_coevolving_core_mouth_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage137_coevolving_core_mouth_sympy_audit.py`

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
**Notes Derivation Review:** Co-evolving map: source normalization, g and S moments, core R = (g-r_F1)^2/(1+r_F1^2), Green's functions T_s and T_q, Boltzmann fixed-point Sigma = exp(-Phi)/Z. R(g_*) = 1/4 verified algebraically. Exact shifted R = 1/4 - dg/sqrt(1+r^2) + dg^2/(1+r^2) correct. Slope identity Pi = Sigma_0(1-RS) from Green's function derivatives. Sigma_0 = (20/9)T^2 from Stage 123.
**Script Review:** R(g_*) = 1/4 via expect_zero. Shifted R formula verified. Linearized slope identity checked (multi-variable expansion). All genuine. All pass (exit code 0). Minor: dead code on line 46.
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The co-evolving core-mouth map is algebraically correct. The exact branch relation `R(g) = (g-r)^2/(1+r^2)` gives `R=1/4` on the compensated point, the shifted defect expands exactly as claimed, and the linearized `Pi` transport law is the proper first-order coupling of the source, core ratio, and traction response.

**Script Review:**

The audit script verifies the exact branch identity, the shifted-ratio formula, and the linearized slope law with symbolic expansion. The saved output matches the note and the algebraic identities reduce cleanly to zero.

**Issues Found:**

None.

---
