# Review: Stage 138 — Frozen traction fixedpoint

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage138_frozen_traction_fixedpoint.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage138_frozen_traction_fixedpoint_sympy_audit.py`

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
**Notes Derivation Review:** Frozen-Sigma_0 fixed point via genuine nonlinear iteration. g_fp ~ 0.6934 (broadened from 0.758). R_fp ~ 0.2827 > 1/4 (drifts off compensation). Pi_fp ~ 1.489 (small downshift from 1.509). Transport law cross-check: delta_R from direct R and from -dg/sqrt(1+r^2)+dg^2/(1+r^2) agree to ~4.5e-16. Physical interpretation sound: broadening pushes R above 1/4 but Pi shift small.
**Script Review:** Genuine fixed-point iteration from canonical seed (400 max iterations, 1e-13 convergence). T_s, T_q Green's functions on grid. Transport law assertion is non-trivial cross-check. All pass (exit code 0).
**Issues Found:** None (minor: delta_g and delta_Pi not explicitly printed; uniqueness not tested from second seed).

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The frozen-traction co-evolving solve is consistent with Stage 137. The fixed point shifts to a broader source profile, the mixed ratio rises above `1/4`, and the mouth bias remains finite.
2. The numerical values in the note match the saved audit output, and the transport-law cross-check is the right nontrivial consistency test.
3. The stage is scoped correctly: it reports the effect of holding the canonical traction fixed, without pretending the canonical compensation has survived exactly.

**Script Review:**

The script genuinely iterates the nonlinear fixed-point map, computes the moments, and verifies the exact transport-law relation numerically. The saved output agrees with the note.

**Issues Found:** None.

**Questions:** None.

---
