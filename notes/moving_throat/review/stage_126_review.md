# Review: Stage 126 — Equal normalized singular limit

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage126_equal_normalized_singular_limit.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage126_equal_normalized_singular_limit_sympy_audit.py`

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
**Notes Derivation Review:** Three-part argument: (1) strict inequality g_Pi < 1 for all finite Pi (from Stage 113 range), so equal-normalized g_c = 1 never reached at finite Pi; (2) g → 1 only as Pi → ∞ (singular limit); (3) traction That diverges as Pi → ∞ (unphysical). Sound and complete.
**Script Review:** Genuine symbolic checks. Output matches all numerical values. All pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The singular-limit argument is correct. The explicit positivity decomposition shows `0 < g_Pi < 1` for every finite positive `Pi`, so the equal-normalized branch cannot occur at finite bias. The endpoint limit `g_Pi -> 1` only appears as `Pi -> infinity`, and the traction divergence follows from the same asymptotics.

**Script Review:**

The audit script performs the inequality decomposition, checks both endpoint limits, and confirms the large-`Pi` behavior of the traction scaling. The saved output matches the stated numerics and asymptotics.

**Issues Found:**

None.

---
