# Review: Stage 066 — Family1 direct operator window

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage066_family1_direct_operator_window.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage066_family1_direct_operator_window_sympy_audit.py`

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

**Notes Derivation Review:** Implicit differentiation dPe_*/dXi correct. Monotonicity chain valid. Pe bounds from Xi*Delta match Stage 61 exactly. Zeta and Pi windows match Stages 63-64. Full cross-stage consistency confirmed.

**Script Review:** End-to-end: Xi → Pe → zeta → Pi for both wall data sets. Implicit diff via sp.solve(sp.diff(...)). All genuine, all match prior stages. Exit code 0.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The implicit-differentiation formula for `dPe_*/dXi` is correct, and the monotonicity conclusion on the stable branch follows as stated.
2. The explicit Family-1 constants reproduce the Stage-61/63/64 windows directly, so the stage is consistent with the earlier threshold chain.
3. The conclusion is appropriately scoped: the explicit branch is nearly saturated, but the residual open question is still the outgoing quadrupole product `Pi_tr`.

**Script Review:**

The script genuinely computes the fixed-point derivative, evaluates the Family-1 support/source constants, and reproduces the direct operator-selected `zeta` and `Pi/C_mix` windows. The output agrees with the notes.

**Issues Found:** None.

**Questions:** None.

---
