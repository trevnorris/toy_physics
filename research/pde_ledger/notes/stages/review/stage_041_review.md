# Review: Stage 041 — Coupled support source operator

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage041_coupled_support_source_operator.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage041_coupled_support_source_sympy_audit.py`

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

1. **Equation-level correctness.** Nondimensionalization consistent with Stage 40. Support-drop kernel K_{kappa,eta}(x) as Green's function difference satisfies homogeneous ODE. Wronskian normalization W = alpha sinh(alpha) + eta cosh(alpha) correct. Kernel derivative verified: all three numerator terms positive on (0,1). Delta_0 = int K dx with formula eta(cosh alpha - 1)/(alpha^2 W) correct. Delta_inf = K(1) formula correct. IVT bracket argument for Pe_* existence sound. Weak-coupling law Pe_* = Xi Delta_0 + O(Xi^2) correct.

2. **Logical flow.** Clean: nondimensional operator → kernel → Delta formulas → IVT bracket → weak-coupling.

3. **Physical interpretation.** Three-way coupling (eta, kappa, Pe) clearly described. Pe now output rather than input.

**Script Review:** Kernel, source, integrals Ic/Is built from primitives. Antiderivative checks genuine (sp.diff). Delta_0 and Delta_inf formulas verified. Bracket endpoints and weak-coupling constant checked. All pass (exit code 0). Comprehensive.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The nondimensionalization and closure `Pe = Xi Delta(Pe; kappa, eta)` are consistent with the Stage-40 support variables and with the Robin/Neumann support problem written here. The kernel `K_(kappa,eta)` and its positive derivative are algebraically correct, so the constructive-branch monotonicity and the bracket `Xi Delta_0 <= Pe_* <= Xi Delta_inf` follow.
2. The endpoint formulas for `Delta_0` and `Delta_inf` match direct evaluation of the kernel integral and the bottom-end kernel value. I independently spot-checked both identities numerically.
3. The weak-coupling statement `Pe_* = Xi Delta_0 + O(Xi^2)` is the right first-order consequence of the small-`Pe` expansion.

**Script Review:**

The audit script faithfully implements the kernel, normalized source family, auxiliary antiderivatives, exact `Delta(Pe; alpha, eta)`, both endpoint identities, and the weak-coupling constant term. I reran the Stage 41 audit locally; it completed successfully and matched the saved audit output.

**Issues Found:** None.

---
