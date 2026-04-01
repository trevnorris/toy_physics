# Review: Stage 166 — Triangular normal form

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage166_triangular_normal_form.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage166_triangular_normal_form_sympy_audit.py`

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

<!-- Agents: append your review below this line using the template: -->

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. **Logical flow from Stage 165:** Stage 165 established the five microscopic slippages (Sigma_Z, Sigma_chi, Sigma_eps, Sigma_del, Sigma_eta), the tracking combination Sigma_tr, the Theta_1 formula, and the Xi_1 formula rewritten in terms of (Sigma_tr, Sigma_Z, Sigma_eps, Sigma_del). Stage 166 takes these as given and reorganizes into branch-adapted coordinates. The transition is clean and correctly cited.

2. **Definition of Sigma_nt (Section 1):** The nontracking slippage is defined by subtracting the A_tr*Sigma_tr piece from Xi_1. I verified the substitution step manually: starting from Stage 165's Xi_1 in terms of (Sigma_Z, Sigma_chi, Sigma_eps, Sigma_del) and replacing Sigma_chi = [Sigma_tr - (1+chi_0)*Sigma_del]/(1+deltaU), the Sigma_del coefficient collects to exactly -[2*chi_0/(1+deltaU) + 4*epsW*deltaU/(11*(1-eps)*(1+deltaU)^2)], matching the boxed formula. The definition is self-consistent.

3. **Triangular normal form (Section 2):** The three boxed equations for (Theta_1, Xi_1, R_1+Xi_1) correctly express the observables as functions of the three branch-adapted coordinates. The triangular structure is genuine: Theta_1 depends only on Sigma_tr, Xi_1 depends on (Sigma_tr, Sigma_nt), and R_1+Xi_1 depends only on Sigma_eta. The interpretation of the hierarchy (tracking, nontracking transfer-shape, dressing) is physically sensible.

4. **Inverse reconstruction (Section 3):** All three inverse formulas are correct.
   - Sigma_tr inversion from Theta_1: straightforward division by C_tr; sign and prefactor correct.
   - A_tr/C_tr = 2*(1+chi_0+deltaU)/deltaU: verified by hand. The common factors (1+chi_0)*(1+deltaU) cancel, and chi_0 cancels, leaving 2*(1+chi_0+deltaU)/deltaU.
   - Sigma_nt = Xi_1 + (A_tr/C_tr)*Theta_1: follows from Xi_1 = A_tr*Sigma_tr + Sigma_nt and Sigma_tr = -(1/C_tr)*Theta_1, giving Sigma_nt = Xi_1 + A_tr/C_tr * Theta_1 (sign correct since Theta_1 = -C_tr*Sigma_tr).
   - Sigma_eta = -[(1-eps_eta)/eps_eta]*(R_1+Xi_1): direct inversion, correct.

5. **Rigidity theorems (Section 4):** The biconditional statements are valid on the physical branch (chi_0 > 0, deltaU > 0, 0 < eps_eta < 1) because the prefactors C_tr and eps_eta/(1-eps_eta) are strictly nonzero there. The triple-rigidity theorem (Section 4.4) follows from the triangular structure plus invertibility.

6. **Notation:** Consistent with Stage 165 throughout. All symbols are either carried forward or explicitly defined.

7. **No missing steps or unjustified assumptions:** The only assumptions are the physical branch conditions (chi_0 > 0, deltaU > 0, 0 < eps_eta < 1), which are stated and necessary for invertibility.

**Script Review:**

1. **Faithful implementation:** The script correctly encodes the Stage 165 formulas for Theta_1, Xi_1, and R_1, the definition of Sigma_nt, the coefficients A_tr and C_tr, all three inverse reconstruction formulas, and the triple-rigidity forward check.

2. **Epsilon definition (line 39):** eps = epsW*(1 - (2/11)*deltaU/(1+deltaU)) matches the notes exactly.

3. **No hardcoded numeric values:** All expressions are built symbolically from the declared parameters. The only numeric constant is Rational(2,11) from the epsilon formula, which is structurally correct.

4. **Near-tautological check (line 81):** The check `Xi1 - (A_tr*SigmaTr + SigmaNT_def) == 0` is close to tautological because SigmaNT_def is constructed from exactly the same terms as Xi1 minus A_tr*SigmaTr. It confirms definition-level consistency but does not independently verify the formula against a separate computation. This is a minor coverage gap, not an error.

5. **Genuine cross-checks:** The inverse reconstruction checks (lines 93-104) are nontrivial -- they substitute the inverse formula back through Theta_1 and verify recovery of SigmaTr, and they combine Xi_1 with the ratio*Theta_1 to verify Sigma_nt recovery. The A_tr/C_tr ratio check (line 98) is also independently computed.

6. **Triple-rigidity (lines 109-114):** Only the forward direction (slippages=0 implies observables=0) is tested. The reverse direction is logically entailed by the already-verified inverse formulas, so coverage is complete in aggregate but not via a direct symbolic test.

7. **Output agreement:** The script output shows exit code 0, all expect_zero checks pass, and the printed formulas match the notes. No discrepancies.

8. **Missing coverage:** The script does not independently verify the substitution from (Sigma_chi, Sigma_delta) into Sigma_tr that was done in Stage 165, but this is within scope of Stage 165's audit, not Stage 166's.

**Issues Found:**

None blocking or substantive.

Minor observations (not requiring action):
- The near-tautological nature of the Sigma_nt definition check (line 81) means it tests self-consistency rather than independent correctness. A stronger check would build Xi_1 from the original five slippages (including Sigma_chi) and then verify the reduction. However, since the substitution was verified in Stage 165 and the script's goal is Stage 166 content, this is acceptable.
- Stage 165 has a minor LaTeX typo on line 405 (`\,rac` instead of `\,\frac`) but this is not a Stage 166 issue.

**Questions:**

None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The triangular normal form and inverse reconstruction are consistent with Stage 165. I rechecked the A_tr/C_tr ratio and the Sigma_eta inversion; both match exactly. The branch-adapted hierarchy of tracking, nontracking, and dressing variables is internally tight.

**Script Review:** The audit script compares the inverse formulas and rigidity statements against the saved output. It runs cleanly and the nontrivial inverse/back-substitution checks pass. The Sigma_nt definition check is close to the formula it is meant to confirm, but the stage still has genuine independent checks, so this does not affect the verdict.

**Issues Found:** None substantive.

---
