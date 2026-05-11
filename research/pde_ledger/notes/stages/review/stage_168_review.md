# Review: Stage 168 — Source map from mode integrals

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage032_source_map_from_mode_integrals.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`

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

1. **Equation-level correctness.** All verified: overlap integrals match prior stages, sigma/kappa_0^2 = 11/9 correct. Reduced coupling structure follows from mode expansions. Schur complement Sigma_wall = Xi I_2 + alpha v v^T correct (Xi from diagonal U-U, alpha from BdG + mixed). Source coupling and projection J_- = g_Q Q_STF (v.e_-) correct. Source map mhat_-^2 = s_-/kappa_0^2. Bound 1 <= mhat_-^2 < 11/9 from monotonic s_- growth. Elimination mhat_-^2 * P_{0,-} = beta_0 s_-^2/(kappa_0^2 lambda_-) by direct substitution.

2. **Logical flow.** Clean: basis → reduced couplings → Schur elimination → source map → bound → elimination.

3. **Assumptions.** All explicit: local isotropic couplings, D/N source branch, same finite-throat basis.

4. **Completeness.** Both bound endpoints handled. Upper bound correctly noted unreachable on stable branch.

5. **Notation consistent** with Stages 12-14.

6. **Physical interpretation.** Sound: source amplification modest (~22%), normalization burden on lambda_- and beta_0.

**Script Review:**

**B.1-B.7.** Faithful: 5 test blocks covering overlaps, mode reductions, 4×4 Schur complement (genuine matrix inversion), source-map endpoints, elimination identity. No bugs. No hardcoded values (constants derived by integration). No tautologies (4×4 inversion is substantive). All pass (exit code 0). Coverage complete.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The finite-throat axial integrals are correct: the `N/N` basis is orthonormal, the `D/N` half-wave overlaps give `kappa_0 = 2 sqrt(2)/pi` and `kappa_1 = -4/(3 pi)`, and the squared ratio `sigma/kappa_0^2 = 11/9` is right.
2. The local-kernel reduction, Schur-complement decomposition, and source projection all line up with the note. I independently checked the overlap integrals and the exact matrix reduction with SymPy.
3. The elimination of the abstract source-map factor is also correct: `mhat_-^2 P_{0,-} = beta_0 s_-^2 / (kappa_0^2 lambda_-)` follows directly from the Stage 13 prefactor formula and the natural source-branch map.

**Script Review:**

1. The audit script is faithful and the saved output matches the stage claims in each block.
2. The matrix inversion in the Schur complement is substantive, and the source-map limit checks are nontrivial.

**Issues Found:**

- **[MINOR] The opening summary is slightly too strong about the source-map bound.** The stage later states the exact stable-branch window as `1 <= mhat_-^2 < 11/9`, but the purpose section phrases it as `1 <= mhat_-^2 <= 11/9`. The strict upper bound is the mathematically precise statement for the branch analyzed here; if the non-strict form is intentional, it should be labeled as a supremum rather than an attained endpoint.

**Questions:**

None.

---
