# Review: Stage 024 — Rank2 support completion

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage024_rank2_support_completion.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage024_rank2_support_sympy_audit.py`

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

1. **Equation-level correctness.** Loaded matrix determinant correctly derived: xi(delta+xi) - m[delta+(1+q^2)xi] - n[delta+(1+r^2)xi] + mn(q-r)^2 = 0. Key observation that determinant is linear in n is correct. n_req formula from solving linear-in-n equation verified. Monotonicity dn_req/dm = -[delta+(1+qr)xi]^2/[...]^2 is negative definite (perfect square with minus sign). Tracking collapse r=q correctly recovers Stage-23 formula G_q - m. Source-tied specialization with q=t*R_U, r=t, t^2=2/9 all correct.

2. **Logical flow.** Clean: determinant construction → support-loading theorem → monotonicity → tracking collapse → source-tied specialization.

3. **Assumptions.** All explicit. Rank-2 loading structure clearly introduced.

4. **Notation consistent** with Stages 22-23.

5. **Physical interpretation.** Sound: tracking vs source-tied fork clearly identified.

**Script Review:**

**B.1-B.7.** Faithful implementation of all 5 notes claims. Determinant from SymPy matrix, n_req from sp.solve, dn/dm from sp.diff — all compared against closed forms. No bugs, no hardcoded values (only lam0=2/9), no tautologies. All pass (exit code 0). Complete coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The rank-2 support determinant is exact. I independently re-derived the `2 x 2` loaded determinant with two rank-1 directions and got the same linear-in-`n` decomposition and exact support theorem `n_req = [xi(delta+xi) - m(delta + (1+q^2)xi)] / [delta + (1+r^2)xi - m(q-r)^2]`.
2. The monotonicity theorem is correct. Differentiating the closed form gives the negative square `d n_req / d m = -[delta + (1 + q r) xi]^2 / [...]^2`, so increasing the mixed baseline always lowers the support load required for the same softening depth.
3. Both specializations are faithful. The tracking case `r=q` collapses back to `G_q - m`, and the source-tied case `r=t`, `q=t R_U` produces the stated support-feasibility window and threshold formulas.
4. This is the right structural extension of Stage 23: it isolates exactly when the support kernel follows the deformed mixed direction versus when it remains tied to the original source direction.

**Script Review:**

1. The script computes the determinant from the matrix, solves for `n_req`, differentiates for monotonicity, and verifies both the tracking and source-tied reductions against the closed forms.
2. The output matches the note, and I did not find a coding or algebraic defect.

**Issues Found:**

None.

**Questions:**

None.

---
