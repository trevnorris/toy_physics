# Review: Stage 043 — Entropic microclosure

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage043_entropic_microclosure.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage043_entropic_microclosure_sympy_audit.py`

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

1. **Equation-level correctness.** Free energy F[sigma,phi] well-constructed: entropic term Theta sigma(log(sigma/sigma_*)-1), coupling -Lambda_phi sigma phi, elastic terms. Chemical potential mu = Theta log(sigma/sigma_*) - Lambda_phi phi correct. Onsager current J = -M sigma partial_s mu expands to drift-diffusion with Einstein D = M Theta. Recovery of Stage 39 exponential family under affine drop: J=0 gives Pe = Lambda Delta_phi/Theta. Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X with three equivalent forms. Passivity/Lyapunov identity via integration by parts: (partial_s mu)J = -J^2/(M sigma) correct.

2. **Logical flow.** Clean: free energy → chemical potential → transport → exponential recovery → Xi identification → passivity.

3. **Physical interpretation.** Sound: microscopic closure resolves mobility/diffusion ambiguity — only Lambda^2/(Theta T_X) appears.

**Script Review:** Chemical potential from sp.diff(f,sigma). Onsager decomposition. J=0 exponential family. Pe-parameterized normalization. Xi_micro three forms. Dissipation density identity. All pass (exit code 0). Minor: SymPy Piecewise/nan artifact in normalization constant (cosmetic, correctly handled via alternative path).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The free energy, chemical potential, Onsager current, and Einstein relation are consistent with the stated positive-density closure. The affine-drop reduction reproduces the exponential zero-flux family exactly, and the identification `Pe = (Lambda_phi/Theta_sigma) Delta phi` is algebraically correct.
2. The microscopic coupling `Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X) = chi_sigma Lambda_phi^2 L^2/T_X` matches the support normalization used in Stage 41 and cleanly explains the earlier phenomenological form via `D_sigma = M_sigma Theta_sigma`.
3. The passivity statement is correct under no-flux boundaries: the integration-by-parts identity yields the standard nonpositive dissipation rate.

**Script Review:**

The script faithfully implements the derivation: chemical potential, Onsager expansion, affine-drop branch, normalized `Sigma_Pe`, coupling identities, dissipation identity, and the support Euler-Lagrange equation all check out. The `Piecewise(..., nan)` normalization artifact is a SymPy edge case, but the alternate normalized branch computation is correct and the audit exits cleanly.

**Issues Found:** None.

---
