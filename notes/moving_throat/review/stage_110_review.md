# Review: Stage 110 — Penetration families

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage110_penetration_families.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage110_penetration_families_sympy_audit.py`

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
**Notes Derivation Review:** Slab g_slab = 2 sin(pi x/2)/(pi x) verified by integration. Exponential g_exp = 2(2+pi x exp(-1/x))/((4+pi^2 x^2)(1-exp(-1/x))) verified. Both families: g(0)→1 (point source), g(∞)→2/pi (uniform). Roots x_*^{slab} ~ 0.798, x_*^{exp} ~ 0.663 in expected range (moderate penetration).
**Script Review:** g_minus from exact radical. Both g functions symbolically defined. Roots by nsolve at 80-digit precision with residuals ~1e-82. No hardcoded roots. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The slab and truncated-exponential penetration families are integrated correctly, and both yield the same lower compensated branch with moderate penetration depths.
2. The limiting behavior is physically consistent: the slab family approaches the point-source and uniform-source bounds in the expected order, and the exact Family-1 target lies between them.
3. The stage’s conclusion is sound: the lower compensated branch is reachable with simple positive localized source laws, without sign changes or oscillatory forcing.

**Script Review:**

The script computes the mouth-bias functions from scratch, solves for the penetration depths with high precision, and checks the residuals directly. I independently rechecked both integrals symbolically; they match the stated formulas exactly.

**Issues Found:** None.

---
