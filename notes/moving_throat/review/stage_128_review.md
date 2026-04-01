# Review: Stage 128 — Mouth branch selection status

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage128_mouth_branch_selection_status.md`
- **Script:** None (status/consolidation stage)

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
**Notes Derivation Review:** All five claims traced: g_Pi range from Stage 113, self-consistent gain law from Stage 125, upper branch impossible from Stage 127, equal-normalized singular limit from Stages 113/127, lower branch uniquely reachable (Pi_* ~ 1.509, T_hat ~ 0.901) from Stages 113/125. Progression from Stage 124 logically consistent. Remaining open items (finite corrections to canonical point) well-scoped.
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The consolidation is faithful to the prior stages: it correctly summarizes the positive mouth-layer family, the self-consistent gain law, the impossibility of the upper branch, the singular nature of the equal-normalized branch, and the regular lower canonical point.
2. The note appropriately narrows the remaining uncertainty to finite corrections around the canonical point rather than reopening the branch-selection question.
3. The progression from Stage 124 through 127 to 128 is logically tight and the stage does not overstate what has been proved.

**Script Review:**

No script is attached for this consolidation note, which is appropriate. The stage is a status summary, and its claims are inherited consistently from the preceding derivations.

**Issues Found:** None.

**Questions:** None.

---
