# Review: Stage 023 — Generalized selected branch

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage023_generalized_selected_branch.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage023_generalized_selected_branch_sympy_audit.py`

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

1. **Equation-level correctness.** Rank-1 loaded 2×2 solve: alpha_req = A_0 xi(delta+xi)/(z_0^2(delta+(1+q^2)xi)) verified from determinant condition. Eigenvector ratio r = q xi/(delta+xi) verified from null-space. Overlap formulas for loading and source channels verified. F_{q,eta} and G_q correct. Split-U specialization with q=-(sqrt(2)/3)R_U, eta=(2/9)R_U correct. Recovery at R_U=1 matches Stage-18/19 exactly. Deformation coefficients H_F, H_G verified by differentiation.

2. **Logical flow.** Clean: generic 2×2 eigenvalue → overlaps → branch functions → split-U specialization → flat recovery → deformation.

3. **Physical interpretation.** R_U < 1 (constructive branch) lowers F_U and raises G_U — correct from sign analysis.

**Script Review:** Constructs full eigenvalue problem from scratch. alpha_req from determinant, eigenvector from first row, overlaps from unnormalized eigenvector. Stage-18/19 recovery by substitution. Deformation by sp.diff. No tautologies. All pass (exit code 0). Thorough coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The rank-1 selected-branch solve is correct. I independently checked the determinant condition and recovered the exact loading formula `alpha_req = A_0 xi (delta + xi) / [z_0^2 (delta + (1 + q^2) xi)]`, together with the lower-eigenvector ratio `e_1/e_0 = q xi/(delta + xi)`.
2. The generalized two-vector normalization law is algebraically sound. The loading/source overlap factors and the combined function `F_(q,eta)` match the symbolic reduction, and the one-dimensional `G_q` function is the correct baseline-loading factor.
3. The split-`U` specialization is consistent with Stage 22. The substitution `q = -(sqrt(2)/3) R_U`, `eta = (2/9) R_U` reproduces the Stage-18/19 flat-`U` branch at `R_U = 1`, and the first-order deformation coefficients `H_F`, `H_G` are correct.
4. Physically, the stage does what it claims: it turns the previously collinear branch into a one-parameter deformation without breaking the selected-branch theorem geometry.

**Script Review:**

1. The script faithfully builds the rank-1 eigenproblem from scratch, computes the overlaps, and checks the flat-branch recovery and first-order deformation by symbolic expansion.
2. The saved output matches the note claims, and I did not find a coding bug or a tautological check.

**Issues Found:**

None.

**Questions:**

None.

---
