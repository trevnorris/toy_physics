# Review: Stage 184 — Selected branch reachability

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage031_selected_branch_reachability.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`

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

1. **Equation-level correctness.** All formulas verified. Overlap derivative ds_-/d alpha_0 = 2 DK^2 KP/R^3 confirmed by explicit differentiation (difference-of-products identity yields 4 KP DK^2, factor of 1/2 gives result). Sign/positivity confirmed (all factors positive). Prefactor derivative dP/dalpha = beta_0(ds/dalpha * lambda_- + s_-^2)/lambda_-^2 > 0 on stable branch. Initial values at alpha_0=0: lambda_-(0)=A, s_-(0)=x0=kappa_0^2, P_{0,-}(0)=beta_0*kappa_0^2/(K_0-Xi_0). Softening threshold alpha_crit correct. Crossing theorem follows from IVT: continuous, strictly increasing, finite start, divergent endpoint.

2. **Logical flow.** Clean: overlap derivative → prefactor monotonicity → starting value → divergence → crossing theorem. Main result follows from constituent lemmas.

3. **Assumptions.** All explicit: DeltaK_ax > 0, KappaProd > 0, beta_0 > 0, lambda_- > 0. Continuity on [0, alpha_crit) follows from lambda_- > 0 there.

4. **Completeness.** Edge case P_target <= P_{0,-}(0) correctly handled by conditional statement.

5. **Notation consistency.** Matches Stages 12-13.

6. **Physical interpretation.** Clear: normalization is a one-dimensional ordered crossing problem on the stable branch.

**Script Review:**

**B.1 Faithful.** Eigenvalue, overlap, prefactor, derivative, initial values, threshold, divergence factorization all implemented.
**B.2 No bugs.** Abstract function approach for quotient-rule check is clean.
**B.3 Hardcoded values.** Line 88 has expanded polynomial for radcrit — would be cleaner computed symbolically, but the check (verifying perfect-square form) is valid and passes.
**B.4 No tautologies.** Overlap derivative: SymPy diff vs closed form. Quotient/HF identity: abstract functions with substitution. Initial values: alpha=0 substitution. Divergence factorization: non-trivial algebraic identity.
**B.5 Symbol assumptions correct.** Same as Stage 13.
**B.6 All pass.** Exit code 0.
**B.7 Coverage.** Overlap derivative, prefactor monotonicity, initial values, threshold, determinant factorization, divergence factorization all covered.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The monotonicity chain is correct. I independently rechecked the overlap derivative, the quotient-rule identity for `dP_{0,-}/d alpha_0`, the zero-loading start value, and the softening-threshold factorization.
2. The crossing theorem follows cleanly from the exact symbolic structure: on the stable branch `P_{0,-}` is continuous, strictly increasing, finite at the origin, and divergent at the threshold.
3. The stage stays aligned with Stage 13 notation and with the saved audit output. No hidden sign or factor error showed up in the symbolic checks.

**Script Review:**

1. The script is faithful to the note. It proves the derivative formula directly, checks the quotient/Hellmann-Feynman identity, and verifies the threshold factorization and stable-side divergence form.
2. The output matches the claims and does not rely on a trivial identity.

**Issues Found:**

None.

**Questions:**

None.

---
