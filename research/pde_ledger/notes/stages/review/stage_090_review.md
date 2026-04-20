# Review: Stage 090 — General dtn deformation

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage090_general_dtn_deformation.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage090_general_dtn_deformation_sympy_audit.py`

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
**Notes Derivation Review:** General deformation Lambda_def = S Lambda_out(beta z) + Sigma_i z^i. L0-L5 coefficients correct. Y_hat via geometric series inversion verified through O(z^5). chi_Q = 3(S beta^5+9 Sigma5)/(3S-Sigma0) correct. Even matching for Sigma2, Sigma4 verified by back-substitution. chi_Q-1 formula correct.
**Script Review:** Series expansion, coefficient extraction, even matching via sp.solve, chi_Q check. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The general isotropic deformation ansatz is expanded correctly. The extracted coefficients `L0`, `L2`, `L4`, and `L5` match the outgoing model with scale `S`, argument rescaling `beta`, and additive self-energy data `Sigma_i`.
2. The normalized branch expansion is right through `O(z^5)`, and the resulting deformation-normalized odd factor `chi_Q = 3(S beta^5 + 9 Sigma5)/(3S - Sigma0)` follows cleanly.
3. Solving the canonical-even matching conditions for `Sigma2` and `Sigma4` gives the stated formulas, and the equivalent `chi_Q - 1` expression is algebraically consistent.

**Script Review:**

The script faithfully reconstructs the deformation model, expands the normalized response, solves the even-matching equations, and verifies the closed-form `chi_Q` expression. I reran the audit locally and independently checked the `chi_Q - 1` rearrangement.

**Issues Found:** None.

---
