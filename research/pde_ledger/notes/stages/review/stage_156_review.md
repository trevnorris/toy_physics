# Review: Stage 156 — Axisymmetric loading mismatch

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage156_axisymmetric_loading_mismatch.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage156_axisymmetric_loading_mismatch_sympy_audit.py`
- **Output:** `scripts/moving_throat/output/moving_throat_pde_stage156_axisymmetric_loading_mismatch_sympy_audit.txt`
- **Prior stage:** `notes/moving_throat/moving_throat_pde_stage155_physical_slope_collapse.md`

## Review Checklist

- [x] Equation-level correctness (signs, factors, indices, limits)
- [x] Logical flow from prior stage(s)
- [x] Assumptions stated and justified
- [x] Notation consistent with prior stages
- [x] Physical interpretation sensible
- [x] SymPy script faithfully implements notes
- [x] Script runs without error
- [x] Script output matches notes claims
- [ ] No missing edge cases or branches (see coverage gap in script review)

## Agent Reviews

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

---

**A. Notes Derivation Review**

**A.1 Equation-level correctness**

Section 1 (physical slope transport). The definitions u2^(A) = -D_{A,2}/D_{A,0}, u4^(A) = (D_{A,2}^2 - D_{A,0} D_{A,4})/D_{A,0}^2, and P_0^(A) = N_{A,0}/D_{A,0} are carried forward correctly from Stage 155 Section 1. Expanding -D2A/D0A = -(D2 + eps lam D21)/(D0 + eps lam D01) to first order in eps and extracting the coefficient of eps lam gives u2^(1) = -(D21 + u2 D01)/D0 with u2 = -D2/D0. This is correct. The P1/P0 formula follows from the quotient rule on N0A/D0A and dividing through by P0 = N0/D0, giving N01/N0 - D01/D0. Verified algebraically.

Section 2 (canonical-even transport). Substituting (u2, u4) = (1/9, 4/81), i.e. D2 = -D0/9 and D4 = -D0/27, into the general u4^(1) formula is correctly carried out. The resulting expression -(5 D01 + 18 D21 + 81 D41)/(81 D0) is algebraically correct; verified by hand-expanding (D2A^2 - D0A D4A)/D0A^2 to O(eps). The hidden-even relation u4^(1) = (8/9) u2^(1) translating to D41 = (2/3) D21 + (1/27) D01 is verified by direct substitution: the residual u4^(1) - (8/9) u2^(1) simplifies to (D01/27 + 2 D21/3 - D41)/D0, which vanishes iff D41 = (2/3) D21 + (1/27) D01.

Section 3 (even-preserving collapse). Setting u2^(1) = 0 in the formula -(D21 + D01/9)/D0 gives D21 = -D01/9. Substituting into the hidden-even operator law: D41 = (2/3)(-D01/9) + (1/27) D01 = -2 D01/27 + D01/27 = -D01/27. Correct. The uniform factorization D_{A,n} = D_n(1 + eps lam delta_D) follows because D21/D2 = (-D01/9)/(-D0/9) = D01/D0 = delta_D, and likewise D41/D4 = (-D01/27)/(-D0/27) = D01/D0 = delta_D. Clean and correct.

Section 4 (static loading mismatch). Xi_load := delta_N - delta_D = N01/N0 - D01/D0. On the even-preserving branch, P1/P0 = N01/N0 - D01/D0 = Xi_load. The lane formulas Delta_Q^(A) = eps lambda_A Xi_load follow from P_0^(A)/P_0 = 1 + eps lambda_A Xi_load, matching the lambda values (1, 1/2, -1). All correct.

Section 5 (microscopic decomposition). The decomposition D_{A,0} = K_A - B_{A,0} - Z_{A,0}, D_{A,2} = -(M_A + B_{A,2} + Z_{A,2}), N_{A,0} = N_A is carried forward from prior stages. The first-order extraction D01 = K1 - B0^(1) - Z0^(1), D21 = -(M1 + B2^(1) + Z2^(1)), N01 = N1 is the natural weak-axisymmetric coefficient. Substituting into the u2^(1) formula with u2 = 1/9: D21 + (1/9) D01 = -(M1 + B2^(1) + Z2^(1)) + (1/9)(K1 - B0^(1) - Z0^(1)) = (1/9)K1 - M1 - (B2^(1) + (1/9)B0^(1)) - (Z2^(1) + (1/9)Z0^(1)). Matches the boxed formula.

Section 6 (portwise logarithmic form). The expansion N_1^(r) = (2 P_r P_{1r})/Delta_r^2 - (2 P_r^2 Delta_{1r})/Delta_r^3 follows from differentiating N_0^(r) = P_r^2/Delta_r^2 with the product/chain rule. The resulting log-slope formula N1/N0 = sum_r w_r (2 P_{1r}/P_r - 2 Delta_{1r}/Delta_r) with w_r = N_0^(r)/sum_s N_0^(s) is the standard weighted-average form. The Z_0^(1) formula follows similarly from Z_0^(r) = Q_r/Delta_r. All correct.

**A.2 Logical flow from prior stage**

Stage 155 concluded with the theorem gate: "compute u2^(1) and P1/P0 from the actual grouped moving-throat response." Stage 156 picks up exactly at this point, computing those slopes from the actual weak-axisymmetric grouped response coefficients D_{A,n} and N_{A,0}. The link is clean. The hidden-even relation u4^(1) = (8/9) u2^(1) from Stage 155 Section 4 is correctly inherited and translated into operator-slope language in Section 2. The canonical branch values (u2, u4) = (1/9, 4/81) and the lambda values (1, 1/2, -1) from Stage 7 are consistently used throughout.

**A.3 Assumptions stated and justified**

The two main assumptions are: (i) canonical compensated even branch (u2 = 1/9, u4 = 4/81), and (ii) weak-axisymmetric Y_{20} deformation pattern with lambda values from Stage 7. Both are explicitly stated. The even-preserving condition u2^(1) = 0 is imposed as a conditional branch, not universally -- this is correctly handled by separating the general formulas (Sections 1-2) from the even-preserving specialization (Section 3).

**A.4 Notation consistency**

Notation is fully consistent with Stage 155. The new objects delta_D, delta_N, Xi_load are cleanly defined. The use of subscript (A) for lanes and superscript (1) for first-order slopes follows the established convention.

**A.5 Physical interpretation**

The central physical interpretation -- that the entire remaining linear grouped defect reduces to a single scalar mismatch between two static logarithmic slopes -- is well stated and follows from the algebra. The portwise decomposition in Section 6 gives concrete microscopic meaning to the two slopes, which is a useful reduction.

---

**B. Script Review**

**B.1 Faithful implementation**

The script implements Sections 1-4 of the notes. It constructs the weak-axisymmetric expansions D0A, D2A, D4A, N0A from symbolic variables (no hardcoded numerics except the canonical branch values 1/9 and 1/27, which are correct). It then:
- Extracts u2^(1), u4^(1), P1 via series expansion and differentiation (lines 42-48).
- Checks the u2^(1) formula against the notes' boxed claim (lines 54-57).
- Substitutes canonical branch values and checks the u4^(1) and P1/P0 formulas (lines 60-75).
- Verifies the hidden-even operator identity and solves for D41 (lines 78-85).
- Imposes u2^(1) = 0 and verifies D21 = -D01/9, D41 = -D01/27 (lines 88-95).
- Computes Xi_load and the lane defects (lines 98-111).

All checks use generic symbolic variables. No tautological assertions found.

**B.2 Code correctness**

The series expansion approach (sp.series to O(eps^2), then differentiate and evaluate at eps=0) is the correct method for extracting first-order coefficients. Dividing by lam after differentiation correctly isolates the slope coefficient. The expect_zero helper properly simplifies and expands before testing. No issues.

**B.3 Hardcoded values**

The only numerical values in the script are sp.Rational(1,9), sp.Rational(1,27), and sp.Rational(8,9), all of which correspond to the canonical branch parameters and the hidden-even relation from prior stages. These are not hardcoded results -- they are input parameters from the physical model. No hardcoded answers or backdoor confirmations.

**B.4 Symbol assumptions**

D0, D01, D2, D21, D4, D41, N0, N01 are declared with `real=True, nonzero=True`. The nonzero assumption is necessary for the series expansion of 1/D0A to be valid (avoids division-by-zero ambiguity in SymPy). u2, u4, P0 are declared `real=True` only, which is appropriate since they are derived quantities. No problematic assumptions.

**B.5 Output agreement**

The output file matches all notes claims:
- u2^(1) general = (-D0*D21 + D01*D2)/D0^2, which equals -(D21 + u2 D01)/D0 with u2 = -D2/D0. Correct.
- u4^(1) canonical = (-5*D01 - 18*D21 - 81*D41)/(81*D0). Matches the notes.
- P1/P0 = N01/N0 - D01/D0. Matches.
- All six expect_zero checks pass (confirmed zero).
- D41 from hidden-even = D01/27 + 2*D21/3. Matches.
- D21 from u2^(1)=0 = -D01/9. Matches.
- D41 on even-preserving = -D01/27. Matches.
- Lane defects: Delta_Q^(20)/eps = Xi_load, Delta_Q^(21)/eps = Xi_load/2, Delta_Q^(22)/eps = -Xi_load. Correct.

**B.6 Coverage**

The script does not audit Sections 5 and 6 of the notes (microscopic decomposition and portwise logarithmic form). These sections introduce additional symbolic objects (K1, M1, B, Z terms; portwise P_r, Delta_r, Q_r decompositions) that would benefit from symbolic verification. This is a coverage gap, though the uncovered algebra is straightforward substitution and standard quotient-rule differentiation, and was verified by hand in this review.

---

**Issues Found:**

No blocking issues. No minor algebraic errors.

One coverage note: the SymPy audit does not cover Sections 5 (microscopic decomposition of D01, D21, N01) or 6 (portwise logarithmic form of N1/N0 and Z0^(1)). These were verified by hand during this review and are correct, but a future script extension could add symbolic checks for completeness.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The axisymmetric loading transport is consistent with the prior stage. The general `u_2^(1)`, `u_4^(1)`, and `P_1/P_0` formulas follow directly from the first-order expansions, the canonical branch specialization is correct, and the even-preserving reduction gives the expected `Xi_load` transport law. The microscopic decompositions and portwise logarithmic forms are consistent with the saved output and the quotient-rule algebra.

**Script Review:**

The script checks the general lane expansions, canonical branch values, hidden-even relation, even-preserving collapse, and the loading-mismatch lanes. I independently checked the core `u_2^(1)` and `P_1/P_0` identities against SymPy; they simplify as expected, and the output matches the notes.

**Issues Found:**

None.

---
