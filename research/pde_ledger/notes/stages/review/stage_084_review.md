# Review: Stage 084 — Natural source map reduction

**Batch:** 13 — Outgoing DtN
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage084_natural_source_map_reduction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage084_natural_source_map_reduction_sympy_audit.py`

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

**Notes Derivation Review:** mhat_0 = 1 on natural branch gives N_Q = 1/chi_Q. chi_Q = 1 implies N_Q = 1 (closure). Small-deviation: N_Q - 1 = -Delta_Q + O(Delta_Q^2). Clean one-step reduction.

**Script Review:** sp.solve for N_Q. Natural branch substitution. Series expansion of 1/(1+Delta_Q). All genuine. All pass (exit code 0). Minor: DeltaQ declared positive (cosmetic, doesn't affect checks).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. Combining the natural source-map branch `mhat_0 -> 1` with Stage 83 is straightforward and correct: the last reduced 2.5PN obstruction becomes purely outgoing, `N_Q = 1 / chi_Q`.
2. The canonical compact outgoing branch then closes the reduced theorem immediately when `chi_Q = 1`, so the remaining finish-line question is rightly phrased as whether the actual passive/outgoing branch realizes that canonical value.
3. The small-defect expansion `N_Q - 1 = -Delta_Q + O(Delta_Q^2)` is the correct linearization of `1 / (1 + Delta_Q) - 1`.

**Script Review:**

The script correctly solves the exact odd-normalization condition for `N_Q`, applies the natural-branch and canonical-branch substitutions, and expands the small outgoing defect. I reran the audit locally and it passed unchanged.

**Issues Found:** None.

---
