# Review: Stage 114 — Parent mouth threshold

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage114_parent_mouth_threshold.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage114_parent_mouth_threshold_sympy_audit.py`

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
**Notes Derivation Review:** Pi_m = V1 L/Theta_sigma from Stage 112. Canonical condition Pi_m = Pi_* from Stage 113 monotone compensation. V1 = Pi_* Theta_sigma/L correct. Split V1 = T_m - q_* A_0'. Linearized deviation law via Taylor expansion of g(Pi) around Pi_*. g'(Pi_*) ~ 0.0714 from Stage 113. All correct.
**Script Review:** Pi_* via nsolve at high precision. V1_star symbolic. g'(Pi_*) via sp.diff. All genuine. All pass (exit code 0). Minor: no back-check of g(Pi_*) = g_minus residual.
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The parent mouth-threshold derivation is consistent with Stages 112-113: the canonical compensation point is `Pi_*`, the implied `V1 = Pi_* Theta_sigma/L` is correct, and the local split `V1 = T_m - q_* A_0'` is just the parent bias in mechanical-plus-Maxwell form. The linearized deviation law around `Pi_*` is also correct.

**Script Review:**

The audit script solves for `Pi_*`, computes `V1_*`, and evaluates the derivative at the compensation point. The independent residual check is effectively zero, and the saved output matches the notes.

**Issues Found:**

None.

---
