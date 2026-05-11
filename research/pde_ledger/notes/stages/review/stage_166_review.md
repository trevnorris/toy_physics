# Review: Stage 166 — Selected mode normalization

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage030_selected_mode_normalization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`

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

1. **Equation-level correctness.** All formulas verified. Eigenvalue formula for 2×2 loaded wall block correct (signs checked via diagonal difference and off-diagonal elements). Hellmann-Feynman overlap s_- = -d lambda_-/d alpha_0 verified. Closed form for s_- matches direct differentiation. Weak-loading limit s_-(0) = kappa_0^2 and strong-loading limit s_- → sigma both verified. Normalized response expansion u2, u4, Gamma5 follow from standard geometric series. Prefactor algebra Gamma_{5,-} = Gamma_2^{port} * P_{0,-} correct. Target equation mhat^2 * P0 = 54 G c_s^5/(5 a^5 c^5) verified. Determinant identity via matrix determinant lemma correct. Critical threshold alpha_crit = AB/(B*x0 + A*x1) correct.

2. **Logical flow.** Clean chain: eigenvalue → overlap → normalization → prefactor → target → determinant. Connects Stage 12 operator-level results to Stage 022 normalized-response language.

3. **Assumptions.** All explicit: DeltaK_ax > 0, KappaProd > 0, beta_0 > 0, lambda_- > 0 on stable branch.

4. **Completeness.** O(omega^6) remainder not addressed, consistent with project convention.

5. **Notation consistency.** Matches Stages 12 and 5. New symbols P_{0,-}, D_{-0}, N_Q^{target} well-defined.

6. **Physical interpretation.** Sound: normalization reduces to spectral ratio — sensitivity of selected eigenvalue to loading.

**Script Review:**

**B.1 Faithful.** Same 2×2 block, eigenvalues, overlap, prefactor, determinant identity, target equivalence.
**B.2 No bugs.** Series expansion, im() extraction, expect_zero all correct.
**B.3 No problematic hardcoding.** Constants 27, 54 are derived (2×27=54).
**B.4 No tautologies.** HF vs closed-form overlap is two independent computations. Determinant identity expands product and compares. Normalization equivalence verifies two condition forms are proportional.
**B.5 Symbol assumptions correct.** A, DK positive; alpha nonneg; x0, x1 positive; beta0, G5 positive.
**B.6 All pass.** Exit code 0.
**B.7 Coverage.** Complete. Strong-loading limit not in script but verified analytically.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The selected-mode normalization chain is internally consistent. The generic response expansion gives the expected `u2`, `u4`, and odd `Gamma5` coefficients, and the selected lower eigenvalue / overlap formulas match the loaded `2x2` wall block exactly.
2. The selected static prefactor `P_{0,-} = beta_0 s_- / lambda_-` is correctly translated into the normalized-response language, and the `mhat_-^2 P_{0,-}` target is algebraically equivalent to the quadrupole normalization condition.
3. The determinant identity and the spectral rewrite are both correct. I independently checked the core normalization identities with SymPy, and the results reduced to zero.

**Script Review:**

1. The script faithfully implements the note structure: series expansion, exact eigenvalues, Hellmann-Feynman overlap, prefactor identity, determinant identity, and target equivalence.
2. The saved output matches the note claims, and the checks are nontrivial rather than tautological.

**Issues Found:**

None.

**Questions:**

None.

---
