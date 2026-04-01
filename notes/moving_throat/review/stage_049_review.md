# Review: Stage 049 — Wall figure of merit

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage049_wall_figure_of_merit.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage049_wall_figure_of_merit_sympy_audit.py`

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

1. **Equation-level correctness.** W_wall = kappa * G_eq^{tw} = 4 pi a^2 L^2 J_1 V0^2/(T_X ell) verified — K_X cancels in product. Threshold translations W_fail = Pe_req/Delta_inf, W_suff = Pe_req/Delta_0 correct. All monotonicity claims trivially correct from explicit product/quotient structure. Constant-compressibility form W_H with J_1 = I_f/H_w correct.

2. **Logical flow.** Clean: define control variable → thresholds simplify → monotonicity → constant-compressibility.

3. **Physical interpretation.** Sound: entire support/source problem reduces to one number against a window.

**Script Review:** Threshold substitution round-trips, monotonicity via sp.diff, constant-compressibility reduction. All pass (exit code 0). Clean and correct.

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The wall figure-of-merit compression is correct: `W_wall = kappa G_eq^(tw)` reduces the explicit branch to one dimensionless number, and the fail/succeed windows follow directly from the Stage-48 thresholds. The monotonicity statements are also correct from the explicit formula.

**Script Review:**

The script matches the notes and the output file. Its threshold round-trips and derivative checks are genuine, and the constant-compressibility reduction is consistent with the earlier wall-layer formulas.

**Issues Found:**

None.

---
