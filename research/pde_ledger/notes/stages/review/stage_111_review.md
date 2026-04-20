# Review: Stage 111 — Mouth source bias status

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage111_mouth_source_bias_status.md`
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
**Notes Derivation Review:** All values verified against Stages 108-110. g bound 0<=g<=1 from positivity theorem (Stage 108). Upper branch g_+ ~ 2.798 excluded, lower g_- ~ 0.758 uniquely selected. Self-matched D/N g=pi/4 ~ 0.785 with 3.61% traction (Stage 109). Slab x_*^{slab} ~ 0.798, exponential x_*^{exp} ~ 0.663 (Stage 110). Boxed conclusion faithful. No overclaiming.
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The status summary is faithful to the earlier stages: positivity confines the mouth bias to `[0,1]`, which excludes the upper branch and leaves the lower compensated branch as the unique physical candidate.
2. The numerical comparisons to the self-matched derivative source and the positive penetration families are consistent with the earlier derivations and do not overstate what is proved.
3. The note appropriately narrows the remaining question to the detailed source profile, not the branch sign.

**Script Review:**

No companion audit script is listed for this status stage. The note is a correct consolidation of the already-verified Stages 108-110 results.

**Issues Found:** None.

---
