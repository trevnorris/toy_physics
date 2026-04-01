# Review: Stage 159 — Outgoing load factorization

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage159_outgoing_load_factorization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage159_outgoing_load_factorization_sympy_audit.py`

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
**Verdict:** PASS (one minor notation item)

**Notes Derivation Review:**

1. **Section 1 (Starting point).** The carry-forward from Stage 158 is consistent. The defect field definition Sigma_r^{(N)} = delta ln(Lambda_r^2 / K) matches the Stage-158 formula Sigma_r^{(N)} = 2 delta ln Lambda_r - delta_K = delta ln(Lambda_r^2 / K), and the load identity Xi_load = sum rho_r^{(N)} Sigma_r^{(N)} on conservative-shape-preserving branches is correctly inherited.

2. **Section 2 (Exact factorization).** Verified by direct algebra. Writing P_r = Omega_{U,r}^2 G_{W,r} + R_r G_{U,r} = Omega_{U,r}^2 G_{W,r} (1 + I_r) and Delta_r = Omega_{U,r}^2 Omega_{W,r}^2 - R_r^2 = Omega_{U,r}^2 Omega_{W,r}^2 (1 - H_r), one gets Lambda_r = P_r/Delta_r = [G_{W,r}(1+I_r)] / [Omega_{W,r}^2 (1-H_r)]. Squaring and dividing by K yields M_r^2 (1+I_r)^2 / (1-H_r)^2. Identity confirmed.

3. **Section 3 (Logarithmic decomposition).** Taking ln of the boxed identity from Section 2 and applying delta immediately gives the three-term logarithmic decomposition. The sign on the (1-H_r) term is correct: ln(Lambda_r^2/K) = 2 ln M_r + 2 ln(1+I_r) - 2 ln(1-H_r), so its variation has the minus sign in front of delta ln(1-H_r). Correct.

4. **Section 4 (First-order transport).** The three logarithmic drifts (dlnM, dlnI, dlnH) are correctly computed from the definitions (verified term by term). The chain-rule step from delta ln(1+I_r) = [I_r/(1+I_r)] delta ln I_r and delta ln(1-H_r) = -[H_r/(1-H_r)] delta ln H_r yields the boxed formula with the correct signs and prefactors.

5. **Section 5 (Expanded transport).** Verified by explicit expansion of the Section 4 formula into primitive variables. All six coefficients (dK, dGW, dGU, dR, dOU, dOW) check out:
   - dK: -1 (from the -1/2 dK in dlnM, times 2)
   - dGW: 2/(1+I) (from dlnM and dlnI cancellation)
   - dGU: 2I/(1+I) (from dlnI only)
   - dR: 2(I/(1+I) + 2H/(1-H)) (from dlnI and dlnH)
   - dOU: -2(I/(1+I) + H/(1-H)) (from dlnI and dlnH)
   - dOW: -2/(1-H) (from dlnM and dlnH cancellation)

6. **Section 6 (Conservative-shape theorem).** Correctly specializes to the case delta ln I_r = 0, delta ln H_r = 0, leaving Sigma = 2 delta ln M_r. The resulting Xi_load formula is consistent.

7. **Section 7 (Square-root mixed-leg law).** The conclusion G_{W,r}/Omega_{W,r}^2 proportional to sqrt(K) is the immediate consequence of delta ln M_r = 0. The sufficient conditions are properly stated as a conjunction of interference rigidity, hybridization rigidity, and the M_r scaling law.

8. **Section 8 (Dominant-port corollary).** Straightforward. Uses sum rho_r = 1, and the approximation rho_{r*} approx 1 collapses the sum to a single term. No issues.

9. **Section 9 (What the stage changes).** Accurate summary of the results and the next theorem gate.

**Notation issue:** Stage 158 uses lowercase g_{W,r}, g_{U,r} in the definition of P_r, while Stage 159 switches to uppercase G_{W,r}, G_{U,r} without any bridging remark or redefinition. The internal consistency of Stage 159 is unaffected (G is used uniformly throughout), but a one-line note acknowledging the relabeling would prevent confusion when cross-referencing.

**Script Review:**

1. **Symbol definitions (line 34).** All six primitive variables (K, OU2, OW2, R, GU, GW) are declared as positive real, which is correct for these physical quantities and avoids spurious sign issues.

2. **Lambda construction (line 36).** Lambda = (OU2*GW + R*GU)/(OU2*OW2 - R^2) faithfully encodes the notes' P_r/Delta_r. The script uses flat (un-subscripted) names for a single generic port, which is appropriate.

3. **Microscopic variables (lines 37-39).** M_r, I_r, H_r are constructed directly from definitions matching Section 2 of the notes. Verified.

4. **Test 1: Exact factorization (lines 41-44).** Tests Lambda^2/K - M^2 (1+I)^2/(1-H)^2 = 0. This is a genuine non-tautological algebraic identity verified by SymPy simplification. It does not hardcode the answer. PASS.

5. **Test 2: First-order perturbation (lines 54-72).** The perturbation method is sound: each variable is multiplied by exp(eps * d_var), the log ratio is series-expanded to O(eps), and the O(eps) coefficient is extracted. This independently computes the first-order defect without assuming the factored form. PASS.

6. **Test 3: Factored formula (lines 74-82).** The microscopic logarithmic drifts dlnM, dlnI, dlnH match the notes' Section 4 exactly. The factored defect formula is then tested against the independently-computed Sigma_exact. This is a non-tautological check. PASS.

7. **Test 4: Expanded formula (lines 84-92).** The expanded primitive-variable coefficients are tested against the same independently-derived Sigma_exact. This directly verifies the notes' Section 5. Non-tautological. PASS.

8. **Test 5: Rigidity corollary (lines 102-103).** Substitutes dlnI = 0 and dlnH = 0 into the factored formula and checks the result equals 2*dlnM. This is a straightforward substitution check, but still useful as it confirms the algebraic reduction is clean. PASS.

9. **Coverage assessment.** The script covers:
   - [x] Exact factorization identity (Section 2)
   - [x] First-order factored defect formula (Section 4)
   - [x] Expanded primitive-variable formula (Section 5)
   - [x] Rigidity corollary (Section 6)
   - [ ] The weighted sum Xi_load (Section 3) -- not tested, but this is just plugging Sigma into a sum, no nontrivial algebra
   - [ ] Dominant-port corollary (Section 8) -- not tested, but trivially follows from the single-term limit

   No critical gap. The untested items are trivial corollaries of tested results.

10. **No hardcoded values or tautological assertions found.** Every check compares an independently-derived quantity against a formula from the notes. The `expect_zero` function uses `sp.simplify(sp.expand(...))` which is adequate.

11. **Output agreement.** All five `expect_zero` calls return 0 in the output file. No errors.

**Issues Found:**

- **MINOR (notation):** Stage 159 switches from lowercase g_{W,r}, g_{U,r} (Stage 158 convention) to uppercase G_{W,r}, G_{U,r} without a bridging note. This is cosmetic -- the algebra is internally consistent -- but adding a one-line remark like "We write G_{W,r} for the coupling previously denoted g_{W,r}" would aid cross-stage readability.

**Questions:**

None. The derivation is clean, the factorization is algebraically verified, and the script is faithful and non-tautological.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The outgoing-load factorization is correct. The exact square-law factorization of `Lambda^2 / K` matches the microscopic definitions, the logarithmic decomposition into `M`, `I`, and `H` terms has the right signs, and the first-order transport law expands to the primitive-variable form stated in the notes. The same algebra yields the conservative-shape and dominant-port corollaries without extra assumptions.

**Script Review:**

The script checks the exact factorization, the first-order perturbation law, the factored and expanded transport formulas, and the rigidity corollary. I also rechecked the core factorization symbolically; it simplifies to zero, and the saved output agrees with the note.

**Issues Found:**

None.

---
