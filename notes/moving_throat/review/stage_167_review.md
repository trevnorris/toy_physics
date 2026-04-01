# Review: Stage 167 — Branch invariant coordinates

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage167_branch_invariant_coordinates.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage167_branch_invariant_coordinates_sympy_audit.py`

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

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The stage is a clean repackaging step that takes the three branch-adapted scalars (Sigma_tr, Sigma_nt, Sigma_eta) from Stage 166 and promotes them to exact logarithmic drifts of three branch composites. Verified the following:

1. **R_tr definition and algebra (Section 1.1).** The two forms of R_tr are algebraically equivalent:
   (1 + chi0/(1+deltaU)) / (1+chi0) = (1+chi0+deltaU) / ((1+chi0)(1+deltaU)).
   Multiply numerator and denominator of the left form by (1+deltaU) to get the right form. Correct.

2. **C_* definition and tracking invariant (Section 2).** Stage 166 gave Sigma_tr = -[(1+chi0)(1+deltaU)(1+chi0+deltaU)/(chi0*deltaU)] * Theta1. The notes define C_* as the magnitude of this prefactor, then T_* = R_tr^{-C_*}, so delta ln T_* = -C_* * delta ln R_tr = -C_* * Theta1 = Sigma_tr. The sign chain is correct: the negative exponent times Theta1 yields the negated product, which matches the Stage-166 inverse formula.

3. **B_* definition and nontracking composite (Section 3).** Stage 166 gave Sigma_nt = Xi1 + [2(1+chi0+deltaU)/deltaU]*Theta1. The notes define B_* = 2(1+chi0+deltaU)/deltaU and N_* = T^2 * R_tr^{B_*}, so delta ln N_* = Xi1 + B_*Theta1 = Sigma_nt. This matches the Stage-166 reconstruction formula at line 238 exactly.

4. **Dressing coordinate (Section 4).** delta ln eps_eta = Sigma_eta is inherited directly from Stage 166. The selected-branch complement identity R_target*T^2/Lambda0 = 1-eps_eta leads to delta ln(1-eps_eta) = -(eps_eta/(1-eps_eta))*Sigma_eta, matching the Stage-166 dressing relation at line 173-183.

5. **Equivalence chain (Section 5).** The three zero-defect conditions in the Sigma variables are equivalent to the three observable-drift conditions by the Stage-166 triangular normal form, now re-expressed as invariance of (R_tr, N_*, eps_eta). The chain of biconditionals is valid because C_* > 0 on the constructive branch (chi0 > 0, deltaU > 0).

6. **Minor typo (Section 6, line 393).** The notes contain "\rac{2..." which is missing a backslash -- should read "\frac{2...". This is cosmetic only and does not affect any derivation.

7. **Assumptions.** The positivity conditions chi0 > 0, deltaU > 0, 0 < eps_eta < 1 are inherited from Stage 166 and stated where needed. No new assumptions are introduced.

8. **Notation.** Consistent with Stage 166 throughout. The new symbols (T_*, N_*, D, E, B_*, C_*) are all cleanly defined at first use with boxed equations.

**Script Review:**

1. **Faithful implementation.** The script constructs B_* and C_* symbolically from (chi0, deltaU) matching the notes definitions. It defines SigmaTr = -C_* * Theta1 and SigmaNT = Xi1 + B_* * Theta1, matching the Stage-166 formulas that Stage 167 inherits.

2. **No hardcoded values.** B_* and C_* are computed from the symbolic variables, not hardcoded numerically. All checks use the symbolic expressions.

3. **No tautological checks.** The script constructs T_* = R_tr^{-C_*} from the exponential parametrization R_tr = Rtr0*exp(small*Theta1), takes the series expansion of the logarithm, and compares the coefficient against SigmaTr. This is a genuine computation -- the series expansion and simplification are non-trivial. Similarly for N_* and the dressing coordinate.

4. **Symbol assumptions correct.** chi0, deltaU, eps_eta are declared positive=True (appropriate for the constructive branch). Theta1, Xi1, SigmaEta are declared real=True. The expansion parameter "small" is declared real=True. These are all correct.

5. **Coverage.** The script tests:
   - The exact product identity R_target*T^2 = Lambda0*(1-eps_eta) [Section 1.3]
   - delta ln T_* = Sigma_tr [Section 2]
   - delta ln N_* = Sigma_nt [Section 3]
   - delta ln eps_eta = Sigma_eta [Section 4]
   - The selected-branch complement identity delta ln E = -(eps_eta/(1-eps_eta))*Sigma_eta [Section 4]
   - The zero-map forward direction (setting Theta1=0, Xi1=0, SigmaEta=0 kills the respective drifts) [Section 5]

6. **Potential gap: inverse direction of equivalence.** The script tests the forward zero map (setting observables to zero kills the branch-invariant drifts), but does not test the reverse direction (setting the Sigma coordinates to zero forces the observables to zero). However, this reverse direction was already established in Stage 166 and is algebraically obvious from the definitions, so this is not a blocking gap.

7. **eps_eta linearization model.** The script models eps_eta_var = eps_eta*(1 + small*SigmaEta), which is a first-order multiplicative perturbation. This correctly represents delta ln eps_eta = SigmaEta at first order. No issue.

8. **Output agreement.** All expect_zero calls return 0. Exit code is 0. The printed intermediate expressions (B_*, C_*, delta ln T_*, delta ln N_*, delta ln eps_eta, delta ln E) are all consistent with the notes.

**Issues Found:**

1. **COSMETIC (notes line 393):** "\rac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}" is missing the leading backslash on \frac. Should be "\frac{2...}". Does not affect correctness.

No blocking or substantive issues found.

**Questions:**

None. The stage is a clean repackaging of Stage-166 results into exact branch-invariant coordinates. No new physical content is introduced; the gain is purely organizational.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The branch-invariant coordinates are correctly derived from Stage 166. I independently verified the R_tr algebra and the exact first-order log-drift identities for T_*, N_*, and eps_eta. The selected-branch complement relation is consistent with the same branch assumptions.

**Script Review:** The saved SymPy output matches the notes, and the series-based checks are genuine rather than mechanical. The only discrepancy I saw is the cosmetic `\rac` typo in the note text, which does not change the math.

**Issues Found:** None substantive.

---
