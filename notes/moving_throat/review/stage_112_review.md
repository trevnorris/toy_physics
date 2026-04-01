# Review: Stage 112 — Mouth boundary layer

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer_sympy_audit.py`

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
**Notes Derivation Review:** Free energy F_mouth with entropic + linear potential terms. Chemical potential mu = Theta ln(sigma/sigma_*) + V_m(z) correct. Onsager current J = -M sigma d_z(mu) correctly expanded. Zero-flux ODE sigma' + (V1/Theta) sigma = 0 gives exponential. Normalized sigma_Pi = Pi exp(-Pi z/L)/(L(1-exp(-Pi))). Pi→0 gives uniform 1/L, Pi→∞ concentrates at mouth. Pi_m = V1 L/Theta_sigma. All verified.
**Script Review:** Profile, normalization integral, zero-flux current, ODE residual — all genuine SymPy computations. All pass (exit code 0). Minor: no limit-behavior checks (Pi→0, Pi→∞).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The GNLS plus localized-Maxwell boundary-layer free energy yields the stated electrochemical potential and stationary exponential profile exactly.
2. The normalized source family `sigma_Pi(z)` is correct, and the interpretation of `Pi_m = V1 L / Theta_sigma` as the single remaining mouth-bias parameter is consistent with the derivation.
3. The stage cleanly upgrades the earlier ad hoc truncated exponential family into a concrete equilibrium law, without changing the branch-selection logic.

**Script Review:**

The script verifies the normalization, zero-flux current, and stationary ODE residual directly. I independently checked the normalization integral; it evaluates to 1 exactly.

**Issues Found:** None.

---
