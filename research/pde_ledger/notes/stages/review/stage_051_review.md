# Review: Stage 051 — Resonance thresholds

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage051_resonance_thresholds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage051_resonance_thresholds_sympy_audit.py`

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

1. **Equation-level correctness.** G_res = C^2 * G_match, W_res = C^2 * W_wall correct. Threshold translation: requiring C^2 W_wall >= Pe_req/Delta_0 gives W_wall >= Pe_req/(C^2 Delta_0). At resonance C^2 = 1/P_res, thresholds become P_res-scaled. Band widths: profile-sensitive zone relative width P_res - 1 ~ 0.56%. All correct.

2. **Logical flow.** Clean: resonance-corrected gain → threshold translation → penalty at resonance → band width.

3. **Physical interpretation.** Correctly notes that resonance family has C^2 < 1 so fail zone is slightly wider. Section 4 handles this properly.

**Script Review:** Symbolic threshold formulas, P_res substitution, band width, relative width checks. All pass (exit code 0). Algebraically lightweight but appropriate for a repackaging stage.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The threshold translation is correct: replacing the matched-branch factor by `C^2(r)` gives `W_res(r) = C^2(r) W_wall`, so the fail/succeed bounds scale by `1 / C^2(r)`.
2. At the resonance point, using `P_res = 1 / C_res^2` correctly yields the `P_res`-shifted thresholds and the narrow `P_res - 1` profile-sensitive band. The stated `0.56%` width is consistent with the benchmark values from Stage 50.
3. The stage does not overclaim: it treats the sech-Gaussian family as a refinement of the existing theorem window, not a replacement for it.

**Script Review:**

The script is faithful and non-tautological for a synthesis stage. It derives the resonance-corrected thresholds symbolically, substitutes `P_res` cleanly, and computes the success/failure band widths exactly. The saved output agrees with the notes.

**Issues Found:** None.

**Questions:** None.

---
