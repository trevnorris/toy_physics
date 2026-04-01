# Review: Stage 085 — Higher odd irrelevance

**Batch:** 13 — Outgoing DtN
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage085_higher_odd_irrelevance.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage085_higher_odd_irrelevance_sympy_audit.py`

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

**Notes Derivation Review:** Generalized denominator with tau_Q omega^7. tau_Q first appears at O(omega^7) in geometric series, so O(omega^5) is tau_Q-free. Series through omega^5 matches Stage 083 identically. At 2.5PN order only chi_Q matters. Five-point reduction ledger faithful.

**Script Review:** Series at O(omega^5) and O(omega^7). tau_Q coefficient at omega^5 = 0 verified. tau_Q coefficient at omega^7 = 1/4 verified (includes cross-term factor of 2). Canonical odd coefficient confirmed. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The generalized retarded one-pole expansion is consistent through `O(omega^5)`: the `tau_Q omega^7` term cannot contribute before `omega^7`, so the 2.5PN-relevant odd structure is controlled only by `chi_Q`.
2. The factorization statement matches the earlier stages: the even conservative defect remains separate from the odd outgoing normalization, and the reduced obstruction is still one scalar `chi_Q`.
3. The finish-line interpretation is correct: at 2.5PN the only live retarded obstruction is the leading outgoing normalization factor, not any higher odd tail.

**Script Review:**

The audit script computes the `omega^5` and `omega^7` series directly and checks the odd coefficients against the canonical `chi_Q` term. I independently rechecked the series expansion; the `tau_Q` contribution is zero at `omega^5` and appears with coefficient `1/4` at `omega^7`, exactly as claimed.

**Issues Found:** None.

---
