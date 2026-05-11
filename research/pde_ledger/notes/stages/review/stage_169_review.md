# Review: Stage 169 — Microscopic normalization equation

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage033_microscopic_normalization_equation.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`

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

1. **Equation-level correctness.** All verified: microscopic abbreviations (A, Delta_0, Chi, beta_0, alpha_0) correctly transcribed from Stage 15. N_- = beta_0 s_-^2/(kappa_0^2 lambda_-) correct. Three-part stability gate: Delta_0 > 0, A > 0, alpha_crit = 9 pi^2 A(A+DK)/[8(11A+9DK)] verified by substituting D/N constants. Monotonicity dN_-/dalpha_0 > 0 from quotient rule + Hellmann-Feynman (all factors positive). Zero-loading onset N_-(0) = beta_0 kappa_0^2/A. Weak-loading first-order coefficient 64(8A+9DK)/(9 pi^4 A^2 DK) verified term by term.

2. **Logical flow.** Clean: abbreviations → normalization equation → stability gate → onset criterion → weak-loading expansion → summary.

3. **Assumptions.** Three stability conditions explicit. Zero-loading start standard. Local isotropic kernel from Stage 15.

4. **Completeness.** Gate complete for 2×2 conservative wall problem. Four microscopic "lanes" correctly identified.

5. **Notation consistent** with Stages 13-15. Chi newly introduced, clearly defined.

6. **Physical interpretation.** Sound: onset criterion = necessary condition, four lanes control different aspects.

**Script Review:**

**B.1-B.7.** Faithful: 5 blocks covering normalization product, monotonicity formula (genuine diff check), alpha_crit simplification, zero-loading onset, weak-loading series, microscopic coupling rewrite with onset stiffness (sp.solve). No bugs. D/N constants as exact symbolic expressions. No tautologies (monotonicity check compares sp.diff output against closed form; onset stiffness uses sp.solve). All pass. Complete coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The microscopic abbreviations are consistent with Stage 15 and the saved audit output: `A`, `Delta_0`, `Chi`, `beta_0`, and `alpha_0` are all carried through correctly.
- I independently rechecked the key algebraic steps. The zero-loading limit gives `N_-(0) = 8 beta_0/(pi^2 A)`, the closed-form stability threshold reduces to `9 pi^2 A(A + DeltaK)/(8(11A + 9DeltaK))`, and the weak-loading coefficient matches the stated `64 beta_0 (8A + 9DeltaK)/(9 pi^4 A^2 DeltaK)`.
- The monotonicity identity is also sound: differentiating the selected-branch product gives the exact formula used in the note, and the positivity assumptions on the factors are sufficient for the intended branch-admissibility argument.
- Physically, the stage does what it claims: it turns the selected-branch normalization problem into one explicit microscopic equation plus a stability gate and an onset test, rather than leaving the source-map factor abstract.

**Issues Found:**

- **[MINOR] The fully substituted branch-admissibility inequality is not independently checked.** The note states the exact microscopic gate `alpha_0 = g_B^2/varpi^2 + Chi^2/(Omega_U^2 Delta_0) < alpha_crit`, but the audit script only verifies the closed form for `alpha_crit` and the onset rearrangement. That is enough for the algebraic formulae, but it leaves the final inequality as a documented assumption rather than a directly validated microscopic bound.

---
