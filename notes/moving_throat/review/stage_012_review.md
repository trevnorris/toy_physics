# Review: Stage 012 — Dynamic loading

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage012_dynamic_loading.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage012_dynamic_loading_sympy_audit.py`

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

1. **Equation-level correctness.** All equations consistent. Overlap constants kappa_0 = 2sqrt(2)/pi, kappa_1 = -4/(3pi), sigma = 88/(9pi^2) match Stage 10. Schur complement decomposes correctly into Xi(omega) I_2 + alpha(omega) v v^T (confirmed symbolically by script). Conservative static limits Xi_0, Delta_0, alpha_0 follow by direct substitution. Isotropic shift cancellation K1 - Xi_0 - (K0 - Xi_0) = DeltaK_ax is trivially correct and verified. tan(2 theta) formula consistent with Stage 11 after substitution. Refined softening threshold correctly uses shifted stiffnesses.

2. **Logical flow.** Clean progression: coupled operator setup → Schur complement → conservative static specialization → outgoing dressing → transfer factor → selected-mode projection. Each step builds on prior stages.

3. **Assumptions.** All stated: minimal internal block, sign convention for U-W cross coupling, passive/outgoing dressing, compact outgoing l=2 branch law.

4. **Completeness.** Conservative static, general dynamic, first-order outgoing, and selected-mode projection all covered. Refined softening threshold included. No missing branches.

5. **Notation consistency.** Consistent with Stages 10-11. New symbols (Xi, alpha(omega), Delta_UW, beta, beta_5, P_0) cleanly introduced.

6. **Physical interpretation.** Sound. Isotropic shift (from U doublet) vs directional rank-1 load (from BdG + mixed) has clear physical meaning. beta_0 >= 0 (ratio of squares) ensures nonnegative transfer factor.

**Script Review:**

**B.1 Faithful.** Internal block Mint and coupling C match notes Section 1. K1 = K0 + DeltaK_ax (recent fix) ensures K1 > K0.

**B.2 No bugs.** Matrix inversions, differentiation, series expansion, limits all correct. Exit code 0.

**B.3 Hardcoded values.** Only kappa_0 = 2sqrt(2)/pi and kappa_1 = -4/(3pi) — exact analytical values derived and verified in Stage 10.

**B.4 Tautological checks.** None. Sigma computed by C^T Mint^{-1} C vs. closed-form Xi*I + alpha*vv^T — independent computations. Beta check compares SymPy derivative vs hand-constructed rational expression. Hellmann-Feynman limits are non-trivial.

**B.5 Symbol assumptions.** All correct: omega real, stiffness/mass/frequency positive, coupling constants real (signs not fixed), Pi unconstrained (can be complex from outgoing port).

**B.6 Output agreement.** All expect_zero pass. All notes claims confirmed.

**B.7 Coverage.** Complete: Schur complement, conservative static data, angle law, softening threshold, outgoing dressing, transfer factor, selected-mode projection with limits.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The main structural claim checks out. Using the reduced internal block `(u, W, phi)`, the wall self-energy does decompose exactly as `Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`; I independently recomputed the Schur complement and got zero difference from that closed form.
- The conservative static limits `Xi_0`, `Delta_0`, and `alpha_0` are consistent with the Stage 11 loaded-angle law, and the isotropic shift cancels out of the wall-level splitting exactly as claimed.
- The outgoing transfer factor is also consistent with the earlier `P/Delta` structure. An independent SymPy check confirms that the `i omega^5` coefficient of the outgoing-substituted `alpha(omega)` is exactly the stated
  `beta_5 = Gamma_2^(port) (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2`.
- The selected-mode projection formula `delta D_-^(odd)(omega) = - i beta_5 kappa(theta_-)^2 omega^5 + O(omega^7)` is the correct rank-1 projection of the odd part onto the conservative lower wall mode.

**Script Review:**

- The script correctly verifies the Schur decomposition, the conservative profile-selection law, the refined softening threshold, and the selected-mode Hellmann-Feynman overlap formula. The saved output is consistent with the notes.
- I did find one coverage gap: after substituting the outgoing port law, the script prints `alpha_out(omega)` and `beta_5`, but it does not explicitly assert that the extracted imaginary `omega^5` coefficient equals the stated `beta_5`. That identity is true, but it is only displayed, not checked.

**Issues Found:**

- **(MINOR)** Missing explicit assertion for the final outgoing-series coefficient. The audit should isolate the `i omega^5` term of `alpha_out(omega)` and compare it programmatically to `beta_5` instead of only printing both objects.

**Questions:**

None.

---
