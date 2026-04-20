# Review: Stage 008 — Minimal isotropic normalization

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (2× PASS, 2026-04-02)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage008_minimal_isotropic_normalization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage008_minimal_isotropic_normalization_sympy_audit.py`

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

**1. Equation-level correctness.** All formulas verified:

- Sec. 1 definitions of Delta, Q, P match Stage 7 exactly (single-mode truncation). B_0 = C^2/varpi^2, Z_0 = Q/Delta, N_0 = P^2/Delta^2, D_0 = K - B_0 - Z_0 all correct.
- P_0 = N_0/D_0 = P^2/[Delta^2(K - C^2/varpi^2 - Q/Delta)]. Algebra verified by multiplying through by Delta.
- Target equation mhat_rad^2 * P_0 = 54 G c_s^5/(5 a^5 c^5) follows from Stage 7 with mhat_ang = 1.
- Stability/positivity (Sec. 3): Delta > 0 (stable 2×2 block), D_0 > 0 (stiffness > softening), N_0 = P^2/Delta^2 >= 0 (perfect square). P_0 > 0 when P != 0. All correct.
- Monotonicity (Sec. 4): partial_K P_0 = -N_0/D_0^2 and partial_X P_0 = +N_0/D_0^2 verified by direct differentiation. Sum = 0. Correct.

**2-6.** Logical flow clean from Stage 7. All assumptions explicit. Five undetermined quantities correctly identified. Notation consistent. Physical interpretation (balance between port transfer and wall self-loading) sound.

**Script Review:**

**B.1-B.7.** Faithful implementation of all formulas. No bugs, no hardcoded values, no tautological checks. Symbol assumptions correct (K/varpi positive, couplings real, X nonneg). All expect_zero checks pass. Coverage includes coefficient definitions, compact P0 equivalence, stability identities, monotonicity derivatives.

**Issues Found:** None.

---

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. **Equation-level correctness.** The single-mode reduction is algebraically consistent. `Delta = Omega_U^2 Omega_W^2 - R^2`, `Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`, `P = Omega_U^2 G_W + R G_U`, and the derived coefficients `B_0`, `Z_0`, `N_0`, `D_0` all match the Stage 7 definitions after truncation. The normalization rewrite `P_0 = P^2 / [Delta (K Delta - Delta C^2 / varpi^2 - Q)]` is correct.

2. **Logical flow.** The stage cleanly narrows the problem from the isotropic angular theorem to the remaining radial/axial scalar. The jump from `mhat_ang = 1` to `mhat_0 = mhat_rad` is justified by the prior stage and the target equation is the right reduced form.

3. **Assumptions.** The stable-branch assumptions `Delta > 0` and `D_0 > 0` are stated explicitly and are the right admissibility conditions for the denominator. The support-softening variable `X = C^2 / varpi^2` is introduced consistently.

4. **Notation consistency.** The notation is consistent with Stage 7 and with the Stage 8 script. The reuse of `P_0` for the normalization prefactor is unambiguous in context.

5. **Physical interpretation.** The interpretation is sensible: stronger bare wall stiffness lowers the prefactor, while stronger support softening raises it. The stage correctly isolates the remaining unknowns to `X`, `Delta`, `Q`, `P`, and `mhat_rad`.

**Script Review:**

The audit script faithfully implements the notes. It recomputes the zero-frequency coefficients, proves the compact form of `P_0`, checks the stability identity `Delta*D_0 = K*Delta - Delta*C^2/varpi^2 - Q`, and differentiates with respect to `K` and `X`. The saved output matches the claims, and my independent SymPy spot-check of `P_0 - P_0_compact` also reduced to zero.

**Issues Found:** None.

**Questions:** None.
