# Review: Stage 147 — Microscopic log channels

**Batch:** 18 — Linear Defect Transport & Final
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage147_microscopic_log_channels.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage147_microscopic_log_channels_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** MINOR
**Notes Derivation Review:** Two channels identified with d ln(g/r) and -d ln(r_c). Expanded into 7 microscopic variables. Uniform-overlap and healing-lock closures applied. Product formulas (27/40 and 27/320 prefactors) verified. delta_perp with 7 numerical coefficients all match. Tangency law correct.
**Script Review:** Section 1-5 all tested. Product formulas checked against independently typed expected forms. Numerical coefficients at Family-1 point. All pass (exit code 0).
**Issues Found:**
1. **(MINOR)** Notes Section 4: intermediate K_s formula has spurious m_psi factor in numerator. Does not affect any downstream formula.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The microscopic log-channel reduction is correct. The two exact channels reduce cleanly to `dln(g/r)` and `-dln(r_c)`, the healing-locked formulas match the reported numerical coefficients, and the final `delta_perp` expression is the right weighted sum of the two branch drifts.

**Script Review:**

The audit script checks the exact channel identities, the uniform-overlap and healing-locked reductions, and the final numeric coefficients. The saved output matches the note, and the branch drift formulas reduce as expected.

**Issues Found:**

None.

---
