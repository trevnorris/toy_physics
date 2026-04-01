# Review: Stage 146 — Off family normal coordinate

**Batch:** 18 — Linear Defect Transport & Final
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage146_off_family_normal_coordinate.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage146_off_family_normal_coordinate_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** Off-family normal coordinate delta_perp from F = 1+r^2-4(g-r)^2. delta R_q = -delta_perp/sqrt(1+r^2). Microscopic two-channel formula. Outlet transport, tangent/normal Pi split. All correct.
**Script Review:** Linearizations, microscopic formula, outlet cross-check, Pi split. Numerical evaluation. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The off-family normal coordinate `delta_perp` is the right defect scalar, and the compensation-ratio transport `delta R_q = -delta_perp/sqrt(1+r_*^2)` is exact.
2. The microscopic expression for `delta_perp` in terms of the parent logarithmic variations is consistent with the note and cleanly isolates the two imbalance channels.
3. The transport into the outlet defects and the tangent/normal split of the mouth-bias deformation follow directly from the linearization.

**Script Review:**

The script verifies the linearized parent defect, the compensation-ratio transport, the microscopic expression, and the outlet transport formulas. I independently checked the key transport identity, and it simplifies exactly to zero.

**Issues Found:** None.

---
