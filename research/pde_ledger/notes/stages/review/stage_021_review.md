# Review: Stage 021 — Reduced Maxwell/mixed one-port normal form

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md`
- **Script:** `scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`

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

**1. Equation-level correctness.** All displayed equations verified:

- **Gauge invariance (Sec. 2):** Gauge transformations on A_0, A_a, A_w applied to E_w and C_a yield cancellation via equality of mixed partials. Correct.

- **Conservative self-energy (Sec. 4):** 2×2 internal block Schur complement gives Sigma = (gA^2 W + 2 gA gW R + gW^2 A) / Delta where Delta = AW - R^2. Standard result, correct. Low-frequency coefficients z0, z2, z4 from geometric series expansion verified.

- **Outgoing dressing (Sec. 5):** N(omega) = (A gW + R gA)^2 / Delta^2 verified by quotient-rule differentiation of Sigma w.r.t. Pi. The perfect-square numerator was confirmed by hand: (-gA^2(AW-R^2) + A(gA^2 W + 2 gA gW R + gW^2 A)) / Delta^2 = (A gW + R gA)^2 / Delta^2. Correct.

- **Sign convention (Sec. 5.1):** delta D = -Pi*N with Pi = +i Gamma omega^5 gives negative imaginary part = positive dissipation. Standard passive convention. Correct.

- **l=2 fingerprint (Sec. 6):** Y_hat_2 from spherical Hankel h_2^(1)(z). Coefficient a^5/(27 c_s^5) for leading imaginary part correct.

- **Scalar compatibility (Sec. 8):** With gA=0, gW=eta*omega, N ~ eta^2 omega^2 A^2/Delta^2. Dangerous i*omega promoted to i*omega^3. Correct.

**2. Logical flow.** Follows naturally from Stage 3 (which left Maxwell/mixed sector as explicit next step).

**3. Assumptions.** All explicit: quadratic/linearized regime, perturbative outgoing dressing, scalar derivative coupling as conditional.

**4. Completeness.** Remaining gaps (multi-mode generalization, full P2 channelization) honestly listed in Sec. 10.

**5. Notation consistency.** Consistent with prior stages. Q_l for wall amplitude, Sigma for self-energy, D_l for wall operator.

**6. Physical interpretation.** Sound. Mixed A_w/F_{mu w} sector as natural outgoing carrier well-motivated. N(0) >= 0 ensures passive transfer.

**Script Review:**

**B.1 Faithful implementation.** All five sections covered: gauge invariance, conservative self-energy, outgoing transfer factor, l=2 fingerprint, scalar compatibility.

**B.2 Correctness of code.** No bugs. `euler_equations` used correctly, series expansions with proper orders.

**B.3 Hardcoded values.** None problematic. Integer coefficients (9, 81, 27) in Y2_hat arise from l=2 Bessel/Hankel and are verified against symbolic computation.

**B.4 Tautological checks.** None found. Gauge invariance tests cancellation from independent gauge transforms. E-L checks test Lagrangian against hand-derived EOM. Solution residuals plug back into original 2×2 system. N(omega) computed by series expansion in Pi. Y2_hat derived from fundamental Bessel/Neumann definitions.

**B.5 Symbol assumptions.** All correct: M,K positive; R,gA,gW real; omega real; Pi real (used as formal expansion parameter); k,a,c_s positive. D0 positive=True enforces stability condition (implicit but justified).

**B.6 Output agreement.** All expect_zero checks pass (exit code 0). Formulas match notes exactly.

**B.7 Coverage.** Complete coverage of all key results. Minor gap: combined wall coefficient (Gamma5_port × N(0)) not algebraically checked but is a trivial product of verified quantities.

**Issues Found:**

None.

---

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The gauge-invariance identities, conservative Maxwell/mixed Schur complement, outgoing transfer factor, compact outgoing `l=2` fingerprint, and scalar derivative-coupling check all line up with the note and the saved audit output.
2. The derivation is logically coherent and the stage stays honest about what is still deferred: the multi-mode generalization and the full grouped-`P2` Maxwell/mixed channelization.
3. I did not find a mathematical inconsistency or a missing assumption that would change the reduced result.

**Script Review:**

1. The script implements the same reduced claims the note makes and checks them independently: gauge invariance, exact elimination of the internal block, low-frequency coefficients, outgoing dressing, the `l=2` Hankel fingerprint, and the scalar derivative-coupling rescue.
2. The symbolic checks are nontrivial and the output matches the stated formulas.

**Issues Found:**

None.

**Questions:**

None.

---
