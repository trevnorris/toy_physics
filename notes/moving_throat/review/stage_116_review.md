# Review: Stage 116 — Coupled mouth fixedpoint

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage116_coupled_mouth_fixedpoint.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage116_coupled_mouth_fixedpoint_sympy_audit.py`

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
**Notes Derivation Review:** D/N BVP particular + homogeneous solution. Mouth derivative yields kernel S(Pi,kappa) = Pi[kappa tanh(kappa)+Pi(exp(-Pi)sech(kappa)-1)]/[(1-exp(-Pi))(kappa^2-Pi^2)]. Static limit S(Pi,0)=1 confirmed by Taylor. Diagonalization of coupled 2×2 into two scalar problems correct. Fixed-point Pi = sum M_alpha S(Pi,kappa_alpha) with M_alpha = L H_alpha G_alpha/Theta_sigma.
**Script Review:** Solution u(x) constructed from A,C formulas. ODE residual, BCs, mouth derivative all verified. S compared against independently typed target. Static limit via sp.limit. All genuine. All pass (exit code 0). Minor: removable singularity at kappa=Pi not discussed.
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The coupled mouth-layer reduction is correct. The explicit D/N solution satisfies the BVP, the mouth derivative kernel matches the stated `S(Pi,kappa)`, the static-shell limit is `1`, and the diagonalized two-channel fixed-point law is the right self-consistency statement for the mouth bias.

**Script Review:**

The audit script verifies the ODE residual, both boundary conditions, the mouth-derivative kernel, and the static limit. The saved output matches the derivation, and the fixed-point law is consistent with the diagonalized two-channel formulation.

**Issues Found:**

None.

---
