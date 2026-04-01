# Review: Stage 124 — Mouth gain status

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage124_mouth_gain_status.md`
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
**Notes Derivation Review:** All formulas verified: M_s and M_q match Stages 121/120. Normalized form M_q = -M_s(g_c-r)^2/(1+r^2) matches Stage 121. Core-balance M_q=-M_s/4 matches Stage 121. Family-1 numerics match Stage 122 (5 sig figs). Self-matched M_s=(20/9)T_hat_m^2 matches Stage 123. ~4% branch difference matches Stage 123. No overclaiming.
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The status summary is faithful to Stages 120–123: the explicit mouth gains are derived from the throat-core ansatz, normalized to the `R_q` ratio, and reduced to the self-matched traction parameter.
2. The Family-1 numerical values match the underlying branch analysis, and the `~4%` traction difference is consistent with the Stage 123 normalization formula.
3. The note correctly narrows the remaining ambiguity to the physical branch choice and its corresponding traction amplitude.

**Script Review:**

No companion audit script is listed for this consolidation stage. The note is consistent with the already-verified derivations and the saved outputs from Stages 122 and 123.

**Issues Found:** None.

---
