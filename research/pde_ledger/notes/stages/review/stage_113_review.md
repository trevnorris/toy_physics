# Review: Stage 113 — Mouth bias map

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage113_mouth_bias_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage113_mouth_bias_map_sympy_audit.py`

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
**Notes Derivation Review:** g_Pi = 2Pi(2Pi exp(Pi)+pi)/[(4Pi^2+pi^2)(exp(Pi)-1)] verified by hand via standard exp*cos integral. Monotonicity via covariance: dg/dPi = -(1/L)Cov_Pi(cos,z) > 0 since cos decreasing and z increasing → negative covariance. Range: g→2/pi (Pi→0, uniform) to g→1 (Pi→∞, delta). Compensation point Pi_* ~ 1.5088 from g_Pi = g_minus via IVT + monotonicity. All correct.
**Script Review:** g_Pi computed from scratch by sp.integrate. Covariance identity verified (residual=0). Both limits computed and asserted. Pi_* by nsolve at 80-digit precision. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The exact mouth-bias factor is correct and matches the earlier penetration-family result after the change of variables `x = 1/Pi`.
2. The monotonicity argument is sound: the covariance identity has zero residual, and the limits `2/pi` and `1` bracket the lower Family-1 compensation target.
3. The compensation root `Pi_*` and the equivalent penetration depth `x_*` are numerically consistent with Stage 110.

**Script Review:**

The script integrates the exponential family directly, verifies the covariance identity, checks both limits, and solves for the compensation point at high precision. I independently checked the exact integral and it reproduces the stated `g_Pi` formula.

**Issues Found:** None.

---
