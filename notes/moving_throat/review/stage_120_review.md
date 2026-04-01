# Review: Stage 120 — Core to mouth gain map

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage120_core_to_mouth_gain_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage120_core_to_mouth_gain_map_sympy_audit.py`

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
**Notes Derivation Review:** Free energy variation gives exponential Boltzmann source. Pi = (L/Theta_sigma)[rho_c U_s'(0) - sigma_c U_q'(0)] consistent with Stage 116. Core coefficients rho_c = g_s^2/K_s and sigma_c = (K_s g_q - lambda g_s)^2/(K_s(K_s K_q + lambda^2)) from Stage 97. M_s = L rho_c/Theta, M_q = -L sigma_c/Theta correct. Fixed-point Pi = M_s + M_q S_q(Pi) matches Stage 117. Sign structure correct (M_s > 0, M_q <= 0).
**Script Review:** rho_c, sigma_c, M_s, M_q symbolically defined. Denominator rewrite equivalence checked. Genuine. All pass (exit code 0). Coverage gaps: free-energy variation and sign properties not tested.
**Issues Found:** None (coverage gaps noted but not errors).

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The electrochemical free-energy setup leads to the stated exponential mouth source law, and the slope parameter `Pi = (L / Theta_sigma)[rho_c U_s'(0) - sigma_c U_q'(0)]` is consistent with the sign structure of the shell and mixed channels.
2. Importing the Stage-97 Schur coefficients gives the quoted explicit mouth gains directly: `M_s = L g_s^2 / (K_s Theta_sigma)` and `M_q = -L (K_s g_q - lambda g_s)^2 / [K_s (K_s K_q + lambda^2) Theta_sigma]`.
3. The resulting fixed-point law `Pi = M_s + M_q S_q(Pi)` is the right specialization of the earlier outlet-consistent mouth closure once the parent core quantities are made explicit.

**Script Review:**

The script faithfully reconstructs `rho_c`, `sigma_c`, `M_s`, and `M_q` from symbolic parent quantities and checks the alternate denominator form of `sigma_c`. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
