# Review: Stage 109 — Positive source families

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage109_positive_source_families.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage109_positive_source_families_sympy_audit.py`

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
**Notes Derivation Review:** Self-matched D/N derivative sigma_match = k cos(kz), g_match = pi/4 ~ 0.785. Traction ratio g_match/g_- ~ 1.0361 (3.61% overshoot). Convex family sigma_xi = (1-xi) k cos(kz) + xi/L normalized and nonneg. g_xi = (1-xi) pi/4 + xi 2/pi linear. xi_* ~ 0.184 solves g_xi = g_-. All verified.
**Script Review:** Normalization by integration. g_match = pi/4 by integration. xi_* by sp.solve. g_xi(xi_*) - g_minus = 0 checked. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The self-matched derivative profile is normalized correctly and gives the exact mouth-bias value `g_match = pi/4`, which is already close to the physical lower compensated branch.
2. The traction comparison is consistent: reaching the exact Family-1 compensated point requires only about a `3.61%` enhancement relative to the self-matched profile.
3. The convex positive family is a clean constructive realization. Its bias is affine in `xi`, and the unique solution `xi_*` that hits `g_-^(F1)` lies inside `(0,1)`, so the compensated branch sits inside a simple positive family rather than requiring sign-changing sources.

**Script Review:**

The script genuinely integrates the matching profile, constructs the convex family, solves exactly for `xi_*`, and verifies that the resulting bias lands exactly on `g_-`. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
