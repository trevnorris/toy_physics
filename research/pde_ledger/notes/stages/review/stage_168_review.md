# Review: Stage 168 — Microscopic monomials

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Hardened and verified (dual-CAS PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage168_microscopic_monomials.md`
- **Script:** `scripts/moving_throat_pde_stage168_microscopic_monomials_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage168_microscopic_monomials_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
- [ ] Script runs without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template: -->

### Agent: Claude Opus 4.6 (1M context) — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

The notes derive three direct microscopic monomials whose logarithmic drifts reproduce the Stage-166 branch-adapted coordinates (Sigma_tr, Sigma_nt, Sigma_eta). The logical structure is clean and the progression from Stage 167 is well-motivated.

1. **Section 1 (microscopic log coordinates):** The five microscopic drift variables (Sigma_chi, Sigma_delta, Sigma_Z, Sigma_eps, Sigma_eta) and their expansions in elementary drifts (gamma1, c1, kU, tau1, lam1, mu1, keta, kW) are consistent with Stage 165 definitions. Verified against the coherent-kernel definitions of chi0, delta_U, epsilon_W, Z_W/Omega_W^2, epsilon_eta.

2. **Section 2 (tracking monomial):** The definition C_{tr,*} = chi0^{1+deltaU*} deltaU^{1+chi0*} is correctly constructed so that d ln C_{tr,*} = (1+deltaU*) Sigma_chi + (1+chi0*) Sigma_delta = Sigma_tr, matching the Stage-166 definition Sigma_tr = (1+chi0) Sigma_delta + (1+deltaU) Sigma_chi. The symmetry in the sum is exact (addition is commutative), so both orderings are equivalent. Correct.

3. **Section 3 (nontracking monomial):** The frozen coefficients E* and F* are the reference-branch evaluations of the Sigma_nt coefficients from Stage 166. The monomial C_{nt,*} = (Z_W/Omega_W^2) eps_W^{E*} delta_U^{-F*} then has d ln C_{nt,*} = Sigma_Z + E* Sigma_eps - F* Sigma_delta = Sigma_nt by construction. The explicit microscopic form at lines 264--278 is consistent with plugging in the coherent-kernel definitions.

4. **Section 4 (dressing monomial):** Trivially correct -- epsilon_eta was already a direct microscopic ratio, so its log drift is Sigma_eta by definition.

5. **Section 5 (observable triangular law):** Re-expresses the Stage-166 triangular relations (Theta1, Xi1, R1+Xi1) in terms of the new monomial coordinates. The coefficients C_{tr,*} and A_{tr,*} match Stage 166. The sign convention Theta1 = -C_{tr,*} Sigma_tr is preserved. Correct.

6. **Section 6 (compatibility ledger):** The three zero-defect equations are correctly obtained by setting Sigma_tr = 0, Sigma_eta = 0, Sigma_nt = 0 and substituting the microscopic drift expansions. The solves for (tau1, kappa_eta, mu1) are algebraically verified:
   - tau1 solve from Sigma_tr=0: correct.
   - kappa_eta from Sigma_eta=0: gives keta = 2c1 - kU. Correct.
   - mu1 from Sigma_nt=0: the initial form and the fully-substituted form (after eliminating tau1, keta) both check out.

7. **Notation consistency:** All symbols match prior stages. The switch from script-T / script-N composites (Stage 167) to fraktur-C monomials (Stage 168) is clearly signaled and motivated.

8. **Physical interpretation:** The claim that the continuation point has been further sharpened -- from "drifts of branch composites" to "invariance of microscopic kernel monomials" -- is justified. The monomial decomposition is a genuine factoring that isolates which microscopic couplings control each defect channel.

**Script Review:**

1. **Faithful implementation of microscopic drifts (lines 37--41):** The five Sigma variables match the notes exactly. Correct.

2. **Branch-adapted coordinates (lines 44--47):** Sigma_tr, E_star, F_star, Sigma_nt all match the notes. Correct.

3. **Tracking monomial test (lines 55--56):** Computes dln_Ctr = (1+deltaUs)*Sigma_chi + (1+chi0s)*Sigma_delta and checks dln_Ctr - Sigma_tr = 0. This is a genuine algebraic identity (both sides are built from the same elementary drifts but through independent formulas). Passes.

4. **Nontracking monomial test (lines 60--61):** Computes dln_Cnt = Sigma_Z + E_star*Sigma_eps - F_star*Sigma_delta and checks dln_Cnt - Sigma_nt = 0. This is structurally the same expression as Sigma_nt itself, so the test is tautological -- both dln_Cnt and Sigma_nt are defined by the same formula on lines 47 and 60. See Issues below.

5. **Dressing monomial test (line 64):** `expect_zero("d ln epsilon_eta - Sigma_eta", Sigma_eta - Sigma_eta)`. This is an outright tautology -- it subtracts Sigma_eta from itself. It does not test anything. See Issues below.

6. **Observable triangular law (lines 67--74):** Correctly defines C_tr_star and A_tr_star matching the notes. Prints Theta1, Xi1, R1+Xi1 but does not assert any identity about them. The Rcombo expression introduces a fresh symbol `epseta_s` that is not tied to any other variables in the script, so it is effectively cosmetic. This section is display-only, not a verification.

7. **Zero-defect solve (lines 77--96):** Uses sp.solve to get tau1 and keta, then constructs mu1 algebraically. The back-substitution tests (lines 87--96) are the strongest part of the script: they substitute the solved values back into Sigma_tr, Sigma_eta, and Sigma_nt and verify all three vanish. These are genuine non-trivial checks. All pass.

8. **No hardcoded values:** The script uses purely symbolic computation throughout. No numerical constants are smuggled in.

9. **Output agreement:** The script output shows all expect_zero checks pass and exit code is 0. The printed formulas for tau1, keta, mu1 are consistent with the notes.

**Issues Found:**

1. [MINOR] **Tautological dressing test (line 64):** The expression `Sigma_eta - Sigma_eta` is identically zero by construction. A non-tautological test would compute d ln epsilon_eta from the microscopic drift definitions independently (e.g., define dln_epseta = 2*c1 - kU - keta and check that against Sigma_eta). Since Sigma_eta is already defined as exactly that expression on line 41, even this is somewhat circular, but the current form is a literal self-subtraction.

2. [MINOR] **Tautological nontracking test (lines 60--61):** The dln_Cnt formula on line 60 is character-for-character identical to the Sigma_nt definition on line 47. The test `dln_Cnt - Sigma_nt = 0` therefore proves nothing about the monomial construction -- it only confirms that the same symbolic expression equals itself. A more meaningful test would construct dln_Cnt from the monomial definition's exponents acting on the elementary Sigma variables independently (e.g., by computing the log-derivative of C_{nt,*} = (Z_W/Omega_W^2) eps_W^{E*} delta_U^{-F*} through explicit partial derivatives).

3. [INFO] **Observable triangular law section is display-only:** Lines 67--74 print the observables but assert nothing. Not a defect per se (the relevant identities are already verified in Stage 166's script), but worth noting for coverage.

**Mitigation:** Despite issues 1 and 2, the zero-defect solve and back-substitution tests (lines 77--96) provide genuine non-trivial coverage for the stage's main claims. The tracking monomial test (line 56) is also meaningful because dln_Ctr and Sigma_tr, while algebraically equivalent, are constructed through different coefficient orderings. The overall verification confidence remains adequate because the back-substitution tests exercise the full system of equations.

**Questions:**

None -- the derivation is clean and the continuation point is clearly stated.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:** The microscopic monomial factorization is mathematically sound. I independently checked the Y20 average and the axisymmetric invariant scaling, and both match the notes. The zero-defect solve and back-substitution are the right continuation of Stage 167.

**Script Review:** The script runs and the back-substitution checks are real, but two of the advertised micro-checks are tautological: `d ln C_nt,* - Sigma_nt` repeats the definition, and `d ln epsilon_eta - Sigma_eta` is a self-subtraction. That makes the audit less independent than the notes suggest, though not wrong.

**Issues Found:** Minor coverage gap only; no mathematical defect.

---

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

The Stage 168 note still states the right local result: the grouped
weak-axisymmetric continuation point is compressed to the first-order invariance
of the three direct microscopic monomials
`(C_tr,*, C_nt,*, epsilon_eta)`, together with the zero-defect compatibility
ledger for the microscopic drifts. It remains a local compiler stage, not an
overclaim of a new global theorem.

**Script Review:**

The historical tautology issue is now repaired in both CAS layers.

The updated SymPy and Mathematica audits now:

1. reconstruct the primitive microscopic ratios for
   `(gamma, c_etaU, T_U, K_U, K_eta^(eff), K_W^(eff), lambda_W, mu_W)`,
2. rebuild `chi_0`, `delta_U`, `epsilon_W`, `Z_W/Omega_W^2`, and
   `epsilon_eta` from those primitive ratios,
3. rebuild the tracking, nontracking, and dressing monomial compilers from the
   primitive ratios rather than only from the carried `Sigma_*` packet,
4. and then rerun the observable complement law and the zero-defect
   substitution checks.

That removes the earlier concern that the nontracking and dressing checks were
just replaying the same symbolic definition. The symbolic theorem path is now
clean in both CAS layers, with the existing `168--170` numerical stress still
present as secondary support.

**Issues Found:**

None. Historical MINOR findings from `2026-04-03` are resolved by the current
scripts.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

No new theorem-path issues surfaced in this pass. The follow-up work here was
operational rather than mathematical: the stage content was already sound, but
the numerical stress harnesses had stale default config paths after the
`research/pde_ledger/` move.

**Script Review:**

The Python and Mathematica numerical harness defaults now resolve the live JSON
fixture directly:

- Python:
  `research/pde_ledger/scripts/numerical/stage168_170_orbit_stress.py`
- Mathematica:
  `research/pde_ledger/mathematica/numerical/stage168_170_orbit_stress.wl`

Both now load
`research/pde_ledger/scripts/numerical/stage168_170_orbit_samples.json`
without requiring a CLI override, and both complete successfully from their
default invocation path.

**Issues Found:**

None in the current Stage `168--170` numerical harness path.

---
