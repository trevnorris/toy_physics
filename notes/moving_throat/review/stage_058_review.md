# Review: Stage 058 — Family1 threshold window

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage058_family1_threshold_window.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage058_family1_threshold_window_sympy_audit.py`

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

**Notes Derivation Review:** All carried values correct from Stages 56-57. Delta_0 ~ 1.733e-4, Delta_inf ~ 2.014e-2 consistent with large-alpha regime. Reduction V0 = alpha_r mu_* with alpha_r=10 giving Upsilon_w = 100 Theta_w correct. Theta thresholds follow trivially.

**Script Review:** Delta formulas from Stage 41 with exact rational kappa. Numerical evaluations match notes. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The explicit-branch numerics are consistent with the earlier stage locks: `Lambda_ell = 37`, `eta = 37`, `kappa = 12321/5`, and the reported `Delta_0` and `Delta_inf` values match the large-`alpha` regime.
2. The threshold window is correctly propagated into both `Upsilon` and `Xi`, and the reduction `Upsilon_w = 100 Theta_w` follows directly from the stated `V0 = alpha_r mu_*` freeze with `alpha_r = 10`.
3. The final `Theta_w` window is a straightforward rewrite of the same threshold band, not a new approximation.

**Script Review:**

The audit script computes the exact branch values, the threshold coefficients, and the `Theta_w` reduction cleanly. The saved output matches the notes, and the independent numeric check reproduces the reported coefficients.

**Issues Found:** None.

---
