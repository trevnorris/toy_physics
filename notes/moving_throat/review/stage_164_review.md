# Review: Stage 164 — Coherent tracking defect

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage164_coherent_tracking_defect.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage164_coherent_tracking_defect_sympy_audit.py`

## Review Checklist

- [x] Equation-level correctness (signs, factors, indices, limits)
- [x] Logical flow from prior stage(s)
- [x] Assumptions stated and justified
- [x] Notation consistent with prior stages
- [x] Physical interpretation sensible
- [x] SymPy script faithfully implements notes
- [x] Script runs without error
- [x] Script output matches notes claims
- [x] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template: -->

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

A. Logical flow from Stage 163:

Stage 163 collapsed the many-port weighted sum into a single effective transfer-shape slope and gave the explicit one-port formula in terms of generic continuum variables (rho, epsilon_W). Stage 164 substitutes the actual coherent local D/N tracking-branch map into that formula, specializing rho -> chi_0 and epsilon_W -> the split blocking ratio epsilon. This is the correct and natural next step. The purpose statement accurately summarizes what was done in 163 and what 164 aims to accomplish.

B. Section 1 -- Exact coherent tracking-branch substitution:

The identification rho_0 = sigma_0 = chi_0 for the coherent branch is consistent with prior coherent-kernel stages. The split blocking ratio epsilon = epsilon_W [1 - (2/11) deltaU/(1+deltaU)] and the selected-branch demand ratio R_target = Lambda (1-epsilon_eta)(1-epsilon)^2 / [Z_W (1+chi_0)^2] are stated as established results.

Substituting R_target into T^2 = (27 pi^2 G c_s^5)/(20 a^5 c^5) * (1-epsilon_eta)/R_target:

  T^2 = [(27 pi^2 G c_s^5)/(20 a^5 c^5)] * (1-epsilon_eta) / [Lambda (1-epsilon_eta)(1-epsilon)^2 / (Z_W (1+chi_0)^2)]
      = [(27 pi^2 G c_s^5)/(20 a^5 c^5)] * Z_W (1+chi_0)^2 / [Lambda (1-epsilon)^2]

Using Lambda = (27 pi^2 G c_s^5)/(20 a^5 c^5) * Omega_W^2, the prefactors cancel and we get:

  T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2]

This is correct. The (1-epsilon_eta) factors cancel as stated.

C. Section 2 -- Support-blindness theorem:

T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2] contains no zeta-dependence. The variables Z_W, chi_0, Omega_W, epsilon (which depends on epsilon_W and deltaU) are all zeta-independent. So d/dzeta T^2 = 0 and d/dzeta ln T^2 = 0 are trivially correct. Similarly R_target = Lambda (1-epsilon_eta)(1-epsilon)^2 / [Z_W (1+chi_0)^2] has no zeta dependence, so d/dzeta R_target = 0. Since Xi_1 = delta ln T_A^2 / (epsilon lambda_A) is the logarithmic slope of T^2 with respect to weak-axisymmetric perturbations of those same zeta-independent variables, d/dzeta Xi_1 = 0 follows. Correct.

The physical interpretation is clear and well-stated: support enhancement changes the baseline M_tr but not the transfer shape, so it helps with steady normalization but cannot repair the weak-axisymmetric grouped defect.

D. Section 3 -- Exact weak-axisymmetric defect law:

From T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2], take the logarithmic derivative:

  delta ln T^2 = delta ln Z_W - delta ln Omega_W^2 + 2 delta chi_0/(1+chi_0) + 2 delta epsilon/(1-epsilon)

Using the perturbation definitions:
  = epsilon lambda_A [zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon)]

So Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon). Correct.

For the split-blocking drift, epsilon = epsilon_W * f(deltaU) where f = 1 - (2/11) deltaU/(1+deltaU). Differentiating:

  epsilon_1 = d(epsilon)/d(epsilon_W) * varepsilon_W + d(epsilon)/d(deltaU) * delta_{U,1}
            = f(deltaU) * varepsilon_W + epsilon_W * f'(deltaU) * delta_{U,1}

  f'(deltaU) = -(2/11) * 1/(1+deltaU)^2

So epsilon_1 = [1 - (2/11) deltaU/(1+deltaU)] varepsilon_W - (2/11) epsilon_W/(1+deltaU)^2 * delta_{U,1}. Correct.

Notation note: The notes use varepsilon_W for the weak-axisymmetric drift of epsilon_W and zeta_Z for the drift of Z_W. Stage 163 used zeta_W for the Z_W drift, not zeta_Z. Stage 164 switches to zeta_Z. This is a minor notation change. The meaning is clear from context and the definitions in Section 3, but worth noting.

E. Section 4 -- Selected-branch reformulation:

The identity Xi_1 = -eta_1/(1-epsilon_eta) - R_1 is carried over from Stage 163. Using R_target = Lambda (1-epsilon_eta)(1-epsilon)^2 / [Z_W (1+chi_0)^2]:

  R_1 = delta ln R_target / (epsilon lambda_A) = omega_W - eta_1/(1-epsilon_eta) - zeta_Z - 2 chi_1/(1+chi_0) - 2 epsilon_1/(1-epsilon)

(since delta ln Lambda = delta ln Omega_W^2 = epsilon lambda_A omega_W at this order, and the (1-epsilon_eta) factor contributes -eta_1/(1-epsilon_eta).)

Check: Xi_1 + eta_1/(1-epsilon_eta) + R_1 = [zeta_Z - omega_W + 2chi_1/(1+chi_0) + 2eps_1/(1-eps)] + eta_1/(1-eps_eta) + [omega_W - eta_1/(1-eps_eta) - zeta_Z - 2chi_1/(1+chi_0) - 2eps_1/(1-eps)] = 0. Correct.

F. Section 5 -- Tracking-factor drift:

R_tr = (1 + chi_0/(1+deltaU))/(1+chi_0) = (chi_0 + deltaU + 1)/[(1+chi_0)(1+deltaU)].

d/dchi_0 ln R_tr = 1/(chi_0 + deltaU + 1) - 1/(1+chi_0) = -deltaU/[(1+chi_0)(1+chi_0+deltaU)]

d/ddeltaU ln R_tr = 1/(chi_0 + deltaU + 1) - 1/(1+deltaU) = -chi_0/[(1+deltaU)(1+chi_0+deltaU)]

Theta_1 = [-deltaU/(1+chi_0)(1+chi_0+deltaU)] chi_1 + [-chi_0/(1+deltaU)(1+chi_0+deltaU)] deltaU_1
        = -[deltaU(1+deltaU) chi_1 + chi_0(1+chi_0) deltaU_1] / [(1+chi_0)(1+deltaU)(1+chi_0+deltaU)]

Wait -- checking the claimed formula more carefully:

  Theta_1 (claimed) = -[chi_0(1+chi_0) deltaU_1 + deltaU(1+deltaU) chi_1] / [(1+chi_0)(1+deltaU)(1+chi_0+deltaU)]

From my computation:
  Theta_1 = [-deltaU chi_1]/[(1+chi_0)(1+chi_0+deltaU)] + [-chi_0 deltaU_1]/[(1+deltaU)(1+chi_0+deltaU)]

Putting over the common denominator (1+chi_0)(1+deltaU)(1+chi_0+deltaU):
  = [-deltaU(1+deltaU) chi_1 - chi_0(1+chi_0) deltaU_1] / [(1+chi_0)(1+deltaU)(1+chi_0+deltaU)]

This matches the claimed formula exactly. Correct.

The observation that Theta_1 = 0 does not imply Xi_1 = 0 is well justified: Xi_1 also involves (zeta_Z, omega_W, varepsilon_W) which are absent from Theta_1.

G. Section 6 -- Summary theorem statement: Correctly summarizes all results. The continuation point is well-defined.

H. Assumptions and completeness: All inputs are clearly traced to prior stages. No hidden assumptions detected. The support-blindness result is exact (not perturbative in zeta), which is correctly emphasized. The stage is complete in the sense that it fully characterizes the defect on the coherent branch and establishes what remains to be computed.

**Script Review:**

1. Symbol setup (lines 36-38): Core symbols are correctly declared with appropriate assumptions. The zeta symbol is included for the support-blindness test.

2. Lambda construction (line 40): Lambda = (27/20) pi^2 G c_s^5 Omega_W^2 / (a^5 c^5). Matches the notes identity Lambda = (27 pi^2 G c_s^5)/(20 a^5 c^5) Omega_W^2. Correct.

3. Split blocking ratio (line 41): eps = epsW * (1 - Rational(2,11) * deltaU/(1+deltaU)). Matches notes. Correct.

4. R_target (line 42): Lam * (1-eps_eta) * (1-eps)^2 / (ZW * (1+chi0)^2). Matches notes. Correct.

5. Transfer-shape identity test (lines 44-47): T2_direct = Z_W(1+chi_0)^2 / [Omega_W^2 (1-eps)^2]; T2_selected = (27 pi^2 G c_s^5)/(20 a^5 c^5) * (1-eps_eta)/R_target. The script tests T2_direct - T2_selected = 0. This is a genuine symbolic identity check, not tautological. PASS in output.

6. Support-blindness tests (lines 48-49): d/dzeta ln T^2 = 0 and d/dzeta ln R_target = 0. These differentiate actual symbolic expressions with respect to zeta. Since none of the constituents depend on zeta, the result is correctly zero. Not tautological -- if zeta had been accidentally included in any expression, this would catch it. PASS in output.

7. Split-blocking drift (lines 54-62): eps1 is computed via symbolic differentiation of eps with respect to epsW and deltaU, then compared to the explicit formula from the notes. This is a genuine derivation-vs-formula check. PASS in output.

8. Xi_1 construction (lines 64-69): Built from the drift components as stated in the notes. The formula is correct.

9. R_1 construction (lines 71-77): Built from the notes formula. The selected-branch identity Xi_1 + eta_1/(1-eps_eta) + R_1 = 0 is tested at line 79. This is a genuine identity check across two independently constructed expressions. PASS in output.

10. Tracking-factor drift (lines 86-92): Theta_1 computed via symbolic differentiation of ln R_tr, compared to the closed-form formula from the notes. Genuine check. PASS in output.

11. Support-rigid specialization (lines 97-103): Sets chi_1 = deltaU_1 = 0 and checks that Xi_1 does NOT vanish (line 102-103). This confirms the notes' claim that tracking-factor rigidity does not kill the grouped defect. Also checks Theta_1 = 0 under this specialization. Both confirmed in output.

12. No hardcoded numerical results: All comparisons are symbolic. Rational constants (2/11, 27/20) match the notes. No floating-point approximations used.

13. No tautological assertions: Each expect_zero call compares two independently constructed expressions. The split-blocking drift test computes eps1 via differentiation and compares to the explicit formula. The transfer-shape identity compares the direct and selected-branch forms. The tracking-factor drift is likewise tested independently.

14. Coverage assessment: The script covers all major results of the stage:
    - Transfer-shape identity (Section 1): covered
    - Support-blindness (Section 2): covered
    - Split-blocking drift (Section 3): covered
    - Defect law Xi_1 (Section 3): constructed, used in identity check
    - Selected-branch reformulation (Section 4): covered via identity check
    - Tracking-factor drift (Section 5): covered
    - Non-sufficiency of tracking-factor rigidity (Section 5): covered

15. Output agreement: All expect_zero checks pass. The printed Xi_1 and R_1 expressions are consistent with the notes formulas (in SymPy's internal representation). Exit code 0.

**Issues Found:**

None of blocking severity.

Minor observations (no action required):

1. Notation shift: Stage 163 uses zeta_W for the Z_W drift, while Stage 164 switches to zeta_Z. Both are defined locally and the meaning is unambiguous, but if a later stage references both, the change in label should be noted. This is a documentation-level concern, not a correctness issue.

2. The script does not independently verify the intermediate formula R_tr = (chi_0 + deltaU + 1)/[(1+chi_0)(1+deltaU)] by checking it equals the original form (1 + chi_0/(1+deltaU))/(1+chi_0), but since the script uses the original form and SymPy handles the algebra, this is implicitly covered.

**Questions:**

None.

---

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The coherent tracking-defect derivation is consistent. The selected-branch transfer-shape identity, the support-blindness result, the split-blocking drift law, the selected-branch defect identity, and the tracking-factor drift formula all check out against the note.

**Script Review:** The audit script independently verifies the transfer-shape identity, the `zeta`-blind derivatives, the split-blocking drift formula, the selected-branch identity, and the tracking-factor drift relation. The saved output matches the note.

**Issues Found:** None.
