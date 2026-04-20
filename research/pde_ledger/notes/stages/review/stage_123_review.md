# Review: Stage 123 — Selfmatched mouth susceptibility

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility_sympy_audit.py`

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
**Notes Derivation Review:** Self-matched closure Theta_sigma = H_w J_s. Sigma_0 = L g_s^2/(K_s Theta_sigma) simplifies to 20 L ell^2 rho_w^2 T_m^2/(9 hbar^2 c_sw^2) — pi a^2 cancels, 4/3 and 5/3 combine to 20/9. Normalized That = rho_w ell sqrt(L) T_m/(hbar c_sw) gives Sigma_0 = (20/9)That^2. That_nat ~ 0.867, That_comp ~ 0.901, fractional enhancement ~4.04%. All verified.
**Script Review:** Full algebraic chain from primitives. Key assertion Sigma_0_hat = (20/9)That^2 is genuine multi-step simplification. Numerical values from Stage 122 inputs. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The self-matched mouth susceptibility closure is algebraically consistent: `Theta_sigma = H_w J_s`, and the `pi a^2` factors cancel to give `Sigma_0 = 20 L ell^2 rho_w^2 T_m^2/(9 hbar^2 c_{s,w}^2)`.
2. The normalized traction parameter `That` reduces the shell gain to `Sigma_0 = (20/9) That^2`, exactly as stated.
3. The computed natural and compensated traction amplitudes differ by about `4.04%`, which matches the note’s qualitative conclusion.

**Script Review:**

The script carries the full symbolic simplification from `Sigma_0` to `(20/9)That^2` and evaluates the two traction amplitudes from the Stage 122 gains. I independently checked the normalization reduction and it is exact.

**Issues Found:** None.

---
