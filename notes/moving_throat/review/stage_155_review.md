# Review: Stage 155 — Physical slope collapse

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage155_physical_slope_collapse.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage155_physical_slope_collapse_sympy_audit.py`

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

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

Stage 155 rewrites the Stage-154 microscopic obstruction amplitudes (K_1, G_1) in terms of the physical grouped-response slopes u_2^{(1)} and P_1, then propagates the consequences through the outlet map and the 2.5PN normalization defect.

1. **Section 1 (Carry-forward physical grouped variables).** The definitions u_2^(A) = -D_{A,2}/D_{A,0}, u_4^(A) = (D_{A,2}^2 - D_{A,0} D_{A,4})/D_{A,0}^2, P_0^(A) = N_{A,0}/D_{A,0} are standard and correctly linearized with the weak-axisymmetric pattern (lambda_20, lambda_21, lambda_22) = (1, 1/2, -1). The grouped-defect decomposition into (a, b) pairs with coefficients epsilon/4 and 3*epsilon/4 is consistent with the three-lane pattern: a = (delta X_20 - delta X_22)/2 and b = (delta X_20 + delta X_22)/2 - delta X_21 type combinations, which with the given lambda pattern give the stated ratios. Correct.

2. **Section 2 (Exact collapse of the Stage-154 obstruction pair).** The first-order variation formulas delta u_2^(A) = -(dD_{A,2} + u_2 dD_{A,0})/D_0 and delta P_0^(A) = (dN_{A,0} - P_0 dD_{A,0})/D_0 follow immediately from differentiating the definitions. The rewrite K_A = -D_0 delta u_2^(A) + (1/9 - u_2) dD_{A,0} is verified algebraically. On the canonical branch u_2 = 1/9, the (1/9 - u_2) term vanishes exactly, giving K_A = -D_0 delta u_2^(A). The G_A = D_0 delta P_0^(A) identity holds for all u_2 (no cancellation needed). Both are confirmed by the SymPy script. Correct.

3. **Section 3 (Weak-axisymmetric amplitude collapse).** Dividing by epsilon * lambda_A gives K_1 = -D_0 u_2^{(1)} and G_1 = D_0 P_1. The grouped-defect forms follow from substituting a_{u,2} = epsilon/4 * u_2^{(1)} etc. Correct.

4. **Section 4 (Physical form of hidden-even consistency relation).** The variation delta u_4^(A) = -(5 dD_0 + 18 dD_2 + 81 dD_4)/(81 D_0) was verified by hand (see review computation below). Under the Stage-153 microscopic even-consistency relation dD_4 = (2/3) dD_2 + (1/27) dD_0, the numerator becomes 8 dD_0 + 72 dD_2 = 8(dD_0 + 9 dD_2), and the result is (8/9)(-(dD_0/9 + dD_2)/D_0) = (8/9) delta u_2. So delta u_4 = (8/9) delta u_2. Confirmed by SymPy. Correct.

5. **Section 5 (Direct outlet coefficients in physical grouped variables).** Substituting K_A = -D_0 delta u_2^(A) into delta kappa_W^(A) = 3(1-sigma_*)/(sigma_* D_0) K_A gives delta kappa_W^(A) = -3(1-sigma_*)/sigma_* * delta u_2^(A). Substituting G_A = D_0 delta P_0^(A) into delta gamma_W^(A) = -(1-sigma_*)/(9 sigma_* N_0) G_A gives delta gamma_W^(A) = -(1-sigma_*)/(9 sigma_*) * delta P_0^(A)/P_0. Both are straightforward substitutions. The weak-axisymmetric amplitude forms kappa_1 and gamma_1 follow by dividing through by epsilon * lambda_A. Correct.

6. **Section 6 (Remaining 2.5PN defect).** Setting delta kappa_W^(A) = 0 forces delta u_2^(A) = 0, then Delta_Q^(A) = -9 sigma_*/(1-sigma_*) * delta gamma_W^(A) = delta P_0^(A)/P_0. On the weak-axisymmetric branch with lambda values (1, 1/2, -1), the three defects Delta_Q^(20) = epsilon P_1/P_0, Delta_Q^(21) = (epsilon/2) P_1/P_0, Delta_Q^(22) = -epsilon P_1/P_0 follow directly. Correct.

7. **Section 7 (Summary and continuation point).** Accurately summarizes the stage results and identifies the narrowed theorem gate. Consistent with all preceding sections.

8. **Logical flow from Stage 154.** Stage 154 established (K_A, G_A) in terms of (dD_{A,0}, dD_{A,2}, dN_{A,0}) and reduced the weak-axisymmetric problem to the scalar pair (K_1, G_1). Stage 155 correctly takes over from there and rewrites these in terms of the physical slopes. The carry-forward definitions of K_A and G_A match Stage 154 exactly (lines 7-8 of Stage 154, lines 122-124 of Stage 155). The outlet map coefficients also match Stage 154's Section 6.

9. **Notation.** Consistent throughout. The introduction of grouped-defect (a, b) notation is new but clearly defined and internally consistent.

**Script Review:**

1. **Faithful implementation.** The script tests all major claims: (a) first-order variation formulas for u_2 and P_0, (b) the obstruction-pair rewrite with the residual (1/9 - u_2) dD_0 term, (c) canonical-even collapse at u_2 = 1/9, (d) the hidden-even consistency relation giving delta u_4 = (8/9) delta u_2, (e) the direct outlet formulas in physical variables, and (f) the Delta_Q identity.

2. **No hardcoded values used as outputs.** The script derives all intermediate quantities symbolically from the definitions. The only hardcoded values are the canonical branch parameters u_2 = 1/9 and u_4 = 4/81, and the Stage-153 even-consistency relation dD_4 = (2/3) dD_2 + (1/27) dD_0, all of which are justified carry-forward from prior stages.

3. **No tautological assertions.** Each expect_zero check compares independently computed quantities: e.g., K_A (built from dD_2, dD_0) is compared against D_0 * delta_u2 (computed via series expansion of the ratio -DA2/DA0). These are genuinely independent computational paths.

4. **Symbol assumptions.** D0, N0, u2 correctly marked nonzero (needed for denominators). eps and lam are real but not marked nonzero -- this is fine since they appear only as linear prefactors and are divided out in the delta_u2/delta_P0 extraction via series expansion. sigma is nonzero (needed for sigma in denominators). No issue.

5. **Series expansion correctness.** The script uses sp.series(..., eps, 0, 2) which captures O(eps^0) and O(eps^1) terms, then subtracts the background value and divides by eps*lam. This correctly extracts the first-order perturbation coefficient. The use of removeO() before algebra is standard.

6. **Coverage.** The script covers Sections 2, 3, 4, 5, and 6 of the notes. Section 1 (definitions and grouped-defect decomposition) is implicitly tested via the downstream checks. Section 3's grouped-defect forms (the a_{u,2}, b_{u,2} etc. identities from Section 1) are not directly tested, but these are simple substitutions of the lambda pattern. The weak-axisymmetric amplitude collapse (K_1 = -D_0 u_2^{(1)}) is tested indirectly through the per-lane identity. This is adequate coverage.

7. **Output agreement.** All six expect_zero checks produce " = 0" in the output. The intermediate printed values for delta u_2 and delta u_4 match the formulas in the notes. No errors or warnings.

**Issues Found:**

None.

**Minor Observations (not blocking):**

- The grouped-defect notation (a_{u,2}, b_{u,2}, etc.) introduced in Section 1 is not directly tested in the script. These are straightforward substitutions of the lambda pattern and do not warrant a separate check, but a future audit could include them for completeness.
- The script does not explicitly test the Section 3 identity K_1 = -(4 D_0 / epsilon) a_{u,2} = -(4 D_0 / (3 epsilon)) b_{u,2}, but this follows trivially from the definitions.

**Questions:**

None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The physical-slope rewrite is correct end to end. The first-order variations of `u_2^(A)` and `P_0^(A)` match the quotient-rule expansions, the Stage-154 obstruction pair collapses exactly on the canonical branch, and the weak-axisymmetric lane reduction gives `K_1 = -D_0 u_2^(1)` and `G_1 = D_0 P_1` as stated. The hidden-even relation also reduces to `delta u_4 = (8/9) delta u_2` with no leftover terms.

**Script Review:**

The audit script independently checks the variation formulas, the canonical collapse, the hidden-even relation, the direct outlet coefficients, and the final `Delta_Q` identities. The saved output matches the notes line by line, and I do not see any tautological check.

**Issues Found:**

None.

---
