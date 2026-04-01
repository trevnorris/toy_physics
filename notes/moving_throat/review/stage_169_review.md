# Review: Stage 169 — Similarity orbit closure

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage169_similarity_orbit_closure.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage169_similarity_orbit_closure_sympy_audit.py`

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

A. Logical flow from Stage 168:
- The notes correctly inherit the three direct microscopic monomials (C_tr,*, C_nt,*, epsilon_eta) and their definitions (E_*, F_*) from Stage 168. Cross-checked against Stage 168 notes -- all definitions, exponents, and coherent-kernel expansions match exactly.
- The transition from "three linear compatibility equations" to "tangent space of a finite similarity orbit" is a natural and well-motivated algebraic completion step.

B. Matrix M_* (Section 2):
- Verified all 24 entries of the 3x8 matrix row by row against the logarithmic derivatives of the three monomials using the coherent-kernel definitions (chi_0 = gamma*c_etaU/K_U, delta_U = pi^2*T_U/(L^2*K_U), etc.).
- Row 1 (tracking): exponent vector of chi_0^{1+delta_*} * delta_U^{1+chi_*}. All 8 entries correct.
- Row 2 (nontracking): exponent vector of (Z_W/Omega_W^2)*epsilon_W^{E_*}*delta_U^{-F_*}. All 8 entries correct. In particular: lam1 coefficient is 2(1+E_*) from the 2 in Z_W/Omega_W^2 plus 2E_* from epsilon_W; kU coefficient is (F_*-E_*) from -E_* in epsilon_W and +F_* from -F_*(-1) in delta_U; kW coefficient is -(2+E_*) from -2 in Z_W/Omega_W^2 and -E_* from epsilon_W.
- Row 3 (dressing): exponent vector of c_etaU^2/(K_U*K_eta^eff) = [0,2,0,-1,-1,0,0,0]. Trivially correct.

C. Minor determinant (Section 2):
- The 3x3 minor selecting columns (tau1, kEta, mu1) = columns (7, 4, 6) is:
  | 1+chi  0   0 |
  |  -F   -1   1 |
  |   0   -1   0 |
  Cofactor expansion along the first row: (1+chi)*[(-1)(0)-(1)(-1)] = 1+chi. Confirmed.
- The conclusion dim ker M_* = 8-3 = 5 follows from rank 3 (guaranteed by the non-vanishing determinant on the constructive branch chi_{0,*} > 0).

D. Finite similarity orbit (Section 3):
- The five free parameters (Lambda, C, Gamma, U, W) scale (lambda_W, c_etaU, gamma, K_U, K_W^eff).
- The three derived scalings for (K_eta^eff, T_U, mu_W) are determined by requiring exact preservation of C_tr,*, C_nt,*, and epsilon_eta. This is a well-posed system (three equations, three unknowns) precisely because the 3x3 minor is nonsingular.
- K_eta^eff exponent: 2C - U. Comes from epsilon_eta = c_etaU^2/(K_U*K_eta^eff) being invariant: c_etaU scales by e^C, K_U by e^U, so K_eta^eff must scale by e^{2C-U}. Correct.
- T_U exponent: U - (1+delta_*)/(1+chi_*)*(Gamma+C-U). Comes from tracking monomial preservation. Verified algebraically.
- mu_W exponent M(Lambda,C,Gamma,U,W): the formula in the notes matches the requirement for nontracking monomial preservation. Verified that the formula is self-consistent.

E. Monomial preservation proofs (Section 4):
- Section 4.1 (tracking): the cancellation is explicit and correct. The (1+chi_*) factor in the second term absorbs the (1+delta_*)/(1+chi_*) from T_U's exponent, giving exact cancellation.
- Section 4.2 (dressing): trivially 2C - U - (2C - U) = 0.
- Section 4.3 (nontracking): follows by construction of M. The notes defer the explicit expansion, which is reasonable given the script verifies it symbolically.

F. Tangent-space equivalence (Section 5):
- The linearization (Lambda,C,Gamma,U,W) -> (lam1,c1,gam1,kU,kW) directly yields the three compatibility formulas from Stage 168 Section 6. Cross-checked against Stage 168: all three formulas match exactly.
- The conclusion ker M_* = T_id G_* follows immediately.

G. Final closure theorem (Section 6):
- The chain of equivalences (Theta_1=Xi_1=R_1=0) <=> (Sigma_tr=Sigma_nt=Sigma_eta=0) <=> (M_* delta x = 0) <=> (delta x in T_id G_*) is correctly assembled from Stage 166, Stage 168, and the current stage's results.

H. Interpretation (Section 7):
- The codimension-3 / dimension-5 orbit characterization is physically sensible: five microscopic co-scalings are free, three are constrained.
- The explicit failure-mode decomposition (rows of M_* corresponding to tracking, nontracking, and dressing failure) is a clean diagnostic.
- The continuation point ("does the true branch lie tangent to G_*?") is correctly identified as the remaining open question.

I. Notation: consistent with all prior stages. No new unexplained symbols.

J. No missing cases: the notes correctly note L is inert in this reduction. The positivity assumption chi_{0,*} > 0 (constructive coherent branch) is stated and is the standing assumption throughout the derivation.

**Script Review:**

A. Faithful implementation:
- The script independently constructs the monomial drift expressions (dlog_Ctr, dlog_Cnt, dlog_eta) from the microscopic kernel definitions (lines 44-46). These are consistent with but algebraically independent from the matrix M_* (lines 52-56). Both are present in the script, providing a cross-check.
- The 3x8 matrix M_* is entered manually and matches the notes' boxed formula entry-by-entry.
- The finite orbit exponents (Eta_exp, Tau_exp, Mu_exp) on lines 101-103 faithfully encode the notes' formulas. The Mu_exp construction uses intermediate variables (Eta_exp, Tau_exp) rather than the notes' closed-form expression, but algebraic expansion confirms equivalence (verified by hand above).

B. No hardcoded values:
- All expressions are purely symbolic. The parameters chi, delta, E, F are declared as free SymPy symbols with appropriate assumptions (chi, delta positive; E, F real). No numerical values are substituted anywhere.

C. No tautological checks:
- The minor determinant check (lines 60-64) is a genuine computation: SymPy computes the 3x3 determinant and the script asserts it equals 1+chi.
- The linear solve (lines 67-71) uses sp.solve on the three drift equations, producing symbolic solutions that are then compared against the notes' closed-form formulas. This is a genuine verification, not a restatement.
- The finite orbit preservation checks (lines 110-116) are partially tautological for Cnt_orbit, since Mu_exp was defined to make Cnt_orbit vanish. However, (a) this still verifies the transcription is correct, and (b) the Ctr_orbit and Eta_orbit checks are genuinely verifying that Tau_exp and Eta_exp were correctly derived.
- The linearization comparison (lines 119-122) is the strongest check: it verifies that the finite orbit's tangent space, computed by substitution, matches the independently-solved linear compatibility formulas. This would catch any transcription error in either the finite orbit or the linear system.

D. Symbol assumptions:
- chi and delta are declared positive=True, which is correct for the constructive coherent branch. This ensures SymPy can simplify 1+chi != 0.
- E and F are declared real=True without sign assumptions. This is appropriate since their signs depend on the branch parameters.

E. Output agreement:
- All expect_zero checks pass (7 total).
- The minor determinant equals chi+1 as claimed.
- The linear solve produces formulas matching the notes exactly.
- The finite orbit preserves all three monomials exactly.
- The linearization reproduces the compatibility ledger exactly.
- Exit code 0, status PASS.

F. Coverage:
- The script covers: matrix construction, rank verification via minor determinant, linear compatibility solve with closed-form comparison, finite orbit construction, exact monomial preservation (all three), and linearization-to-ledger equivalence.
- One thing the script does not check: it does not verify that the kernel of M_* is exactly 5-dimensional by computing the full rank of M_*. However, the nonsingular 3x3 minor already implies rank >= 3, and since M_* is 3x8, rank = 3 follows, giving dim ker = 5. This is logically complete even without an explicit rank computation.

**Issues Found:**

None.

**Questions:**

None. This is a clean and well-structured final stage that correctly identifies the geometric content of the Stage-168 compatibility ledger as the tangent space of an exact finite similarity orbit.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The rank-3 compatibility matrix, the nonvanishing 3x3 minor, and the finite similarity orbit closure all line up with Stage 168. I independently checked the minor determinant and the orbit-preservation identities; they are consistent with the claimed 5-dimensional kernel / orbit picture.

**Script Review:** The SymPy audit verifies the matrix entries, determinant, solve, orbit preservation, and tangent-space equivalence, and the saved output reports zero residuals. I saw no hidden assumptions or unsupported branches.

**Issues Found:** None.

---
