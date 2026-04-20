# Review: Stage 108 — Positive source theorem

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage108_positive_source_theorem.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage108_positive_source_theorem_sympy_audit.py`

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
**Verdict:** MINOR
**Notes Derivation Review:** Mouth-bias g[sigma] = integral sigma cos(pi z/(2L)) dz. Positivity bound 0 <= g <= 1 correct (cos in [0,1] on [0,L], sigma >= 0, normalized). g_+^{F1} ~ 2.798 > 1 excluded; g_-^{F1} ~ 0.758 in (0,1) admissible. Branch values internally consistent (37^2*3 = 4107 verified). Reduction logic sound.
**Script Review:** Branch values computed and position relative to [0,1] checked. Output matches notes. Moderate coverage — positivity bound stated in prose not verified symbolically; balance relation not cross-checked.
**Issues Found:**
1. **(MINOR)** Notes reference "Stage 105" for Family-1 branch values but stages 104-107 don't exist in repository. Traceability gap.
2. **(MINOR)** Positivity bound (core theorem) not symbolically verified in script.
3. **(MINOR)** Balance relation 1+r^2 = 4(g-r)^2 not cross-checked for branch values.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The theorem itself is correct. For any nonnegative normalized source profile on `[0,L]`, the mouth-bias factor is the cosine moment `g[sigma]`, so the bound `0 <= g[sigma] <= 1` follows immediately from `0 <= cos(pi z / (2L)) <= 1`.
2. With the carried Family-1 branch values, the physical conclusion is right: `g_-` lies inside `(0,1)` while `g_+ > 1`, so the upper compensated branch is excluded by positivity and the lower one is the only admissible positive-source branch.
3. The issue is not the algebra but the presentation trail: the note cites a nonexistent Stage `105` for the branch values, so the provenance is wrong even though the numbers themselves are consistent.

**Script Review:**

The script reruns cleanly and checks the numeric location of `g_-` and `g_+` relative to the positive-source window, but it leaves the core positivity theorem and the underlying branch relation at the prose level rather than verifying them symbolically.

**Issues Found:**

1. **(MINOR)** The note references Stage `105`, but stages `104-107` do not exist in the repository. That is a real traceability error.
2. **(MINOR)** The script does not independently verify the core positivity bound or the carried branch relation; it only checks the final numeric inequalities.

---
