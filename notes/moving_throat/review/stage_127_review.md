# Review: Stage 127 — Unique regular canonical branch

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage127_unique_regular_canonical_branch.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage127_unique_regular_canonical_branch_sympy_audit.py`

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
**Notes Derivation Review:** Four-section argument: (1) g_+ ~ 2.798 > 1 → impossible by positive-source bound; (2) 2/pi < g_- ~ 0.758 < 1 → unique Pi_* exists by IVT + Stage 113 monotonicity; (3) finite moderate traction at Pi_* ~ 1.509 (That ~ 0.901) vs self-matched Pi_match ~ 1.908 (That ~ 1.011); (4) combined theorem: upper impossible, equal-normalized singular (Stage 126), lower uniquely regular. All numerics match Stages 113/125 to 14+ sig figs.
**Script Review:** All symbolic objects from Stages 113-125. Both nsolve roots (Pi_*, Pi_match). Ordering assertion genuine. 30-digit precision. All pass (exit code 0). Minor coverage gaps: no explicit g_+ > 1 or 2/pi < g_- < 1 assertions.
**Issues Found:** None (minor script coverage gaps noted).

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The branch-selection logic is correct: `g_+^{F1} > 1` is incompatible with the positive-source bound, while `g_-^{F1}` lies strictly inside the finite-mouth interval `(2/pi, 1)`.
2. The numerical values and ordering are consistent with the note and the saved audit output. The finite canonical point `Pi_*` and the comparison point `Pi_match` are both regular and correctly ordered.
3. The stage is scoped properly: it identifies the unique regular canonical branch without claiming the full PDE branch-selection theorem is already proved.

**Script Review:**

The script genuinely computes the compensated branches, the mouth-layer kernel, the canonical root, and the comparison point. The saved output matches the note, and the ordering assertion is a real check rather than a restatement.

**Issues Found:** None.

**Questions:** None.

---
