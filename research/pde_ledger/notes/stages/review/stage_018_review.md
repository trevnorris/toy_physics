# Review: Stage 018 — Dimensionless normalization locus

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage018_dimensionless_normalization_locus.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage018_dimensionless_normalization_locus_sympy_audit.py`

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

1. **Equation-level correctness.** All verified: eta = kappa_1^2/kappa_0^2 = 2/9 correct. Shape function F(xi,delta) = (9d+11xi)^4/[81(1-xi)(9d^2+18dxi+11xi^2)^2] derived from Stage 17 N_-(x) by substituting D/N constants and factoring N_-(0). F(0,delta) = 1 verified. F→∞ as xi→1 correct. Softening constant C(delta) = (9d+11)^4/[81(9d^2+18d+11)^2] correct. Monotonicity: dF/dxi numerator polynomial has all positive coefficients for d>0, xi>0. Required loading alpha_req = 9 pi^2 A xi(xi+d)/[8(9d+11xi)] verified. Near-onset F = 1 + (1+8/(9d))xi + O(xi^2) confirmed. Near-softening 1-xi_req ~ C(delta)/R_target correct.

2. **Logical flow.** Clean: D/N constants → dimensionless variables → shape function → monotonicity → uniqueness theorem → asymptotics.

3. **Assumptions.** D/N constants from earlier stages. delta > 0 stated.

4. **Completeness.** Three R_target cases handled. Both near-onset and near-softening limits given.

5. **Notation consistent.** New symbols (xi, delta, F, R_target, C) cleanly introduced.

6. **Physical interpretation.** Sound: shape function encodes universal branch geometry.

**Script Review:**

**B.1-B.7.** Faithful: F derived from Stage 17 formulas with D/N constants (genuine derivation check, not restatement). Monotonicity derivative checked via sp.diff. Limits verified. Series expansions via sp.series. 8 expect_zero checks, all pass. No tautologies. D/N constants as exact rationals. Minor: xi declared positive=True (excludes xi=0) but F(0) check works via explicit substitution. Near-softening asymptotic not tested (trivial consequence).

**Issues Found:**

1. **(MINOR, script)** Near-softening asymptotic (Sec 6.2) not explicitly tested. Trivial consequence of verified C(delta).

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The reduction is algebraically exact. I independently rechecked the substitution of the D/N constants into the Stage 17 normal form and got the same closed form for `F(xi,delta)`, the same monotone derivative, and the same `R_target` normalization ratio.
2. The endpoint and asymptotic statements are consistent with the closed form: `F(0,delta)=1`, `F -> +infinity` as `xi -> 1^-`, and `lim (1-xi)F = C(delta)` matches the stated softening constant.
3. The required total loading and support-loading formulas are consistent with the earlier split `alpha_0 = alpha_mix + g_B^2/varpi^2`. The stage is doing the right reduction work and cleanly isolates the uniqueness condition `R_target > 1`.

**Script Review:**

1. The audit script faithfully derives `F`, `dF/dxi`, `alpha_req`, `alpha_crit`, the support-loading formula, and the near-onset series from symbolic definitions. The saved output matches the note exactly on those checks.
2. I also ran a small independent SymPy verification of the closed-form `F - F_target` identity and the manifestly positive derivative form; both reduced to zero.

**Issues Found:**

- **[MINOR] The branch parameter `xi` is declared `positive=True` even though the note’s stable branch includes the boundary `xi = 0`.** The script still checks `F(0,delta)` explicitly, so the result is not wrong, but the symbol assumptions should be relaxed to `nonnegative=True` or the boundary convention should be stated explicitly.

**Questions:**

None.

---
