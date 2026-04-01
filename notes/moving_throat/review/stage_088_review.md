# Review: Stage 088 — ChiQ fix from outgoing dtn

**Batch:** 13 — Outgoing DtN
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage088_chiQ_fix_from_outgoing_dtn.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage088_chiQ_fix_from_outgoing_dtn_sympy_audit.py`

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

**Notes Derivation Review:** Retarded module with chi_Q expanded. sigma_can = 4a^5/(27 c_s^5) verified. omega^5 imaginary coefficient = chi_Q a^5/(27 c_s^5). Matching to Stage 087 DtN fingerprint i a^5/(27 c_s^5) fixes chi_Q = 1. Deformed DtN parametrization: xi_Q deformation gives chi_Q = xi_Q. Clean three-step derivation.

**Script Review:** sigma_can computed and verified. Retarded module series expansion. chi_Q = 1 from sp.solve. Deformed DtN coefficients checked. All genuine. All pass (exit code 0).

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The `chi_Q` fixing is algebraically consistent. With `sigma_Q^{can} = 4 a^5/(27 c_s^5)`, the retarded module reproduces the outgoing `\omega^2`, `\omega^4`, and `i \omega^5` coefficients, and matching the Stage 87 fingerprint fixes `chi_Q = 1` exactly. The deformed DtN parametrization also behaves as claimed.

**Script Review:**

The script checks the canonical `\Omega_Q` substitution, the `sigma_Q^{can}` normalization, the retarded series expansion, the exact solve for `chi_Q`, and the deformed DtN coefficient map. The saved output matches the notes and the checks are genuine.

**Issues Found:**

None.

---
