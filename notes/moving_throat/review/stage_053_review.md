# Review: Stage 053 — Gnls wall shell

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage053_gnls_wall_shell.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage053_gnls_wall_shell_sympy_audit.py`

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

**Notes Derivation Review:**

1. **Equation-level correctness.** GNLS quadratic expansion standard. Shell norms N_phiphi and G_phiphi correctly defined. T_X and K_X from mode separation correct. kappa = 4(m c_sw L/hbar)^2 + (I_g/I_f)(L/ell)^2 — both terms verified by dividing K_X L^2/T_X. W_wall collapse: a and I_f cancel exactly, giving W_wall = 4 rho_w^2 V0^2 L^2/(hbar^2 c_sw^2 ell^2). Xi = W_wall identification on matched thin-wall branch correct.

2. **Logical flow.** Clean: GNLS shell reduction → kappa → W_wall collapse → Xi = W_wall identification.

3. **Physical interpretation.** Sound: independence from a and I_f is a genuine simplification.

**Script Review:** T_X, K_X, kappa, W_wall, Xi all built from first principles and checked against closed forms. No tautologies. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The wall-shell reduction is consistent with the GNLS quadratic density energy: projecting `delta rho = q(s) chi_phi(y)` gives the stated `T_X` and `K_X`, and the thin-shell measure produces the `N_(phi phi)` and `G_(phi phi)` formulas used later.
2. The reduced `kappa` expression is correct. One term comes from wall compressibility and the other from transverse curvature, with the expected `L^2 / ell^2` enhancement in the thin-shell part.
3. The collapse of `W_wall` and the identification `Xi = W_wall` on the matched thin-wall branch are algebraically sound. The cancellations of `a` and `I_f` are real, not cosmetic.

**Script Review:**

The script faithfully rebuilds the shell norms, support coefficients, `kappa`, `J_1`, `W_wall`, and `Xi` from symbolic inputs, then checks the stated closed forms. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
