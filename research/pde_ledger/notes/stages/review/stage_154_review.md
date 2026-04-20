# Review: Stage 154 — Microscopic grouped obstructions

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage154_microscopic_grouped_obstructions.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage154_microscopic_grouped_obstructions_sympy_audit.py`

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

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

**A.1 Equation-level correctness.**

All derivations verified by hand:

- **K_A decomposition (Section 1).** Starting from Stage 153 definition K_A = delta D_{A,2} + (1/9) delta D_{A,0}, substituting the full-bundle expansions delta D_{A,0} = delta K_A - delta B_{A,0} - delta Z_{A,0} and delta D_{A,2} = -delta M_A - delta B_{A,2} - delta Z_{A,2} yields K_A = (1/9 delta K_A - delta M_A) - (delta B_{A,2} + 1/9 delta B_{A,0}) - (delta Z_{A,2} + 1/9 delta Z_{A,0}). Confirmed: all signs and 1/9 factors correct. The three-bundle decomposition W_A, B_A, Z_A follows by definition.

- **G_A decomposition (Section 1).** Starting from G_A = delta N_{A,0} - P_0 delta D_{A,0} and substituting delta D_{A,0} = delta K_A - delta B_{A,0} - delta Z_{A,0} gives G_A = delta N_{A,0} - P_0 delta K_A + P_0 delta B_{A,0} + P_0 delta Z_{A,0} = -P_0 delta K_A + P_0 delta B_{A,0} + N_A, where N_A := delta N_{A,0} + P_0 delta Z_{A,0}. Correct.

- **BdG linearization (Section 2).** For B_{A,0} = c^2/w^2: d/dc gives 2c/w^2, d/dw gives -2c^2/w^3. For B_{A,2} = c^2/w^4: d/dc gives 2c/w^4, d/dw gives -4c^2/w^5. All signs and powers verified. The combined B_A = delta B_{A,2} + (1/9) delta B_{A,0} correctly collects dc and dw coefficients: 2c(1/w^4 + 1/(9w^2))dc - 2c^2(2/w^5 + 1/(9w^3))dw.

- **Maxwell/mixed linearization (Section 3).** The quotient-rule variations delta Z_{A,0} = (Delta dQ - Q dDelta)/Delta^2 and the four-term delta Z_{A,2} formula follow from standard linearization of Z_0 = Q/Delta and Z_2 = (QS - G*Delta)/Delta^2. Verified all four coefficient terms in Z_{A,2}: S/Delta^2 for dQ, Q/Delta^2 for dS, -1/Delta for dG, and (G/Delta^2 - 2QS/Delta^3) for dDelta. The combined Z_A^(r) = delta Z_{A,2} + (1/9) delta Z_{A,0} correctly yields the dDelta coefficient G_r/Delta_r^2 - Q_r/(9 Delta_r^2) - 2Q_r S_r/Delta_r^3.

- **Transfer bundle (Section 4).** From N_{A,0}^(r) = P_r^2/Delta_r^2, the variation gives 2P/Delta^2 dP - 2P^2/Delta^3 dDelta. Correct by quotient/power rule. The combination N_A^(r) = delta N_{A,0}^(r) + P_0 delta Z_{A,0}^(r) correctly cancels to give three-term form with P_0/Delta dQ + 2P/Delta^2 dP - (P_0 Q/Delta^2 + 2P^2/Delta^3) dDelta. Verified: delta S and delta G correctly drop out of the odd transfer bundle.

- **Primitive variations (Section 5).** All five variations (delta Delta, delta S, delta G, delta P, delta Q) verified against the definitions. delta Q has six terms from the quadratic form Q = g_U^2 W + 2 g_U g_W R + g_W^2 U; the partial derivatives are all correct.

- **Weak axisymmetric reduction (Section 6).** Substituting delta X_A = epsilon lambda_A X^(1) into both K_A and G_A microscopically decomposed forms yields the scalar amplitudes mathfrak K_1 and mathfrak G_1. The formulas for kappa_1 and gamma_1 correctly use the Stage 153 prefactors 3(1-sigma_*)/(sigma_* D_0) and -(1-sigma_*)/(9 sigma_* N_0). All lambda values (1, 1/2, -1) correctly inherited.

**A.2 Logical flow.** Clean progression: Section 1 substitutes the full-bundle decomposition into the Stage 153 obstructions, Sections 2-4 linearize each microscopic sector, Section 5 unpacks to primitive variables, Section 6 specializes to the weak axisymmetric branch, Section 7 summarizes the new theorem status. No steps skipped.

**A.3 Assumptions.** All assumptions explicitly stated: isotropic compensated branch for linearization (same as Stage 153), standard chain-rule linearization of all nonlinear expressions, portwise summation structure for Maxwell/mixed sector.

**A.4 Notation.** Consistent with Stages 6 and 153. New symbols (W_A, B_A, Z_A, N_A for the four microscopic bundles; mathfrak K_1, mathfrak G_1 for the weak scalar amplitudes) are clearly defined and do not collide with prior notation.

**A.5 Physical interpretation.** Sound. The decomposition into four microscopic defect bundles (wall baseline, BdG support, conservative Maxwell/mixed, outgoing transfer) provides clear physical attribution. The observation that delta S and delta G drop out of the odd transfer bundle is a nontrivial structural result with clear meaning: port-frequency sums and coupling-norm sums cannot generate transfer anisotropy at linear order.

**A.6 Completeness.** Both main obstructions (K_A and G_A) fully decomposed. All three microscopic sectors (BdG, Maxwell/mixed conservative, Maxwell/mixed transfer) treated. Primitive variations for the Maxwell/mixed port explicitly unpacked. Weak axisymmetric specialization included. Continuation point clearly stated.

**Script Review:**

**B.1 Faithful implementation.** The script covers all major formulas from the notes: the K_A/G_A decomposition (Section 1 of script), BdG linearization (Section 2), Maxwell/mixed portwise linearization including all primitive variations (Section 3), and the weak axisymmetric reduction (Section 4). Each boxed formula in the notes has a corresponding symbolic check.

**B.2 Code correctness.** All symbolic expressions use SymPy Rational for exact fractions (sp.Rational(1,9), sp.Rational(1,2)). The expect_zero helper correctly expands and simplifies before comparing to zero, and raises AssertionError on failure. No floating-point approximations.

**B.3 No hardcoded values.** The script does not hardcode any final answers. All results are derived from the symbolic definitions using SymPy's diff() for the linearization checks. The BdG variations (lines 56-57) are computed via sp.diff(B0, c)*dc + sp.diff(B0, w)*dw, not entered by hand. Similarly, the Maxwell/mixed variations (lines 87-94, 115-120) are computed from the symbolic expressions Delta_expr, S_expr, Q_expr, G_expr, P_expr, Z0, Z2, N0expr using SymPy differentiation.

**B.4 No tautological checks.** The script genuinely tests each claim. For example, K_exact and K_split are built from different algebraic arrangements, and the check K_exact - K_split = 0 verifies that the regrouping is correct. The BdG and Maxwell/mixed sections compute variations via SymPy's diff and compare against the manually written formula forms. The weak-axisymmetric section (lines 158-164) builds both the microscopically decomposed form and the direct substitution form and checks they agree.

**B.5 Symbol assumptions.** Appropriate: c, w declared nonzero (BdG frequencies and couplings); U, W, R declared nonzero (port frequencies and mixing); P0 declared nonzero; Delta, S, Q, G, P declared nonzero (port invariants). Variation symbols (dc, dw, dU, dW, dR, dgu, dgw, dDelta, dS, dQ, dG, dP) correctly left without positivity assumptions.

**B.6 Output agreement.** All 21 expect_zero checks return 0 in the output file. No assertion failures. Carry-forward formulas printed at the end match the notes.

**B.7 Coverage.** The script covers all four major derivation sections. The N_A^(r) combination formula (Section 4 of notes, combining delta N_{A,0}^(r) + P_0 delta Z_{A,0}^(r) into the three-term form) is tested symbolically on lines 143-148 as the "N obstruction bundle" check. The primitive-variation formulas including the full six-term delta Q expression from Section 5 are tested on lines 100-107. Coverage is complete.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The decomposition of the two grouped outlet obstructions into wall, BdG, conservative Maxwell/mixed, and transfer pieces is algebraically correct, including the `1/9`-weighted combinations carried from Stage 153.
2. The BdG linearization and the Maxwell/mixed primitive variations are all consistent with direct differentiation of the underlying microscopic expressions.
3. The weak-axisymmetric specialization is also correct: the grouped outlet obstruction really collapses to the scalar pair `(K_1, G_1)` once the lane dependence is reduced to a single amplitude.

**Script Review:**

The script is strong here. It symbolically differentiates the BdG and Maxwell/mixed primitives, verifies every bundled obstruction formula, and checks the weak-axisymmetric reduction at several representative `lambda` values. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
-->
