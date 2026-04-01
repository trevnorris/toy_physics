# Review: Stage 054 — Tanh wall branch

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage054_tanh_wall_branch.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage054_tanh_wall_branch_sympy_audit.py`

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

1. **Equation-level correctness.** Tanh wall f(xi) = (1+tanh xi)/2: I_f = 1/3 (integral of sech^4/4 = 1/3), I_g = 4/15 (verified via tanh substitution), I_g/I_f = 4/5. Specialization of Stage-53 formulas correct. Mouth closure K_m = T_X/ell gives eta = L/ell = Lambda_ell. Three-parameter reduction to (chi_s, Lambda_ell, Upsilon_w) correct.

2. **Logical flow.** Clean: wall moments → Stage-53 specialization → mouth closure → parameter reduction.

**Script Review:** f, f', f'' symbolic. I_f by sp.integrate. I_g via tanh substitution (correct but hand-performed — minor, documented). kappa_reduced and W_wall_reduced checked. All pass (exit code 0).

**Issues Found:**

1. **(MINOR)** I_g computed via hand-performed tanh substitution rather than direct integration of f''(xi)^2. Correct and documented, but direct integral would be stronger independent check.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The canonical wall choice `f(xi) = (1 + tanh xi)/2` and the resulting support profile `chi_phi = f'(xi) = (1/2) sech^2 xi` are consistent with Stage 53. The moments `I_f = 1/3`, `I_g = 4/15`, and `I_g / I_f = 4/5` are the correct specialization.
2. The local mouth closure `K_m = T_X / ell` gives `eta = K_m L / T_X = L / ell` exactly, so the branch really does collapse to the three parent ratios advertised in the notes.
3. The reduced formulas `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` and `W_wall = Upsilon_w Lambda_ell^2` follow directly from the Stage-53 identities with the tanh-wall moments inserted.

**Script Review:**

The script matches the note structure and the saved output. Its `I_g` computation uses a documented tanh-substitution form rather than the raw `f''(xi)^2` integral, but that is still a genuine symbolic derivation, not a hardcoded answer. I reran the audit locally and it passed.

**Issues Found:** None.

---
