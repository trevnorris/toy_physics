# Review: Stage 087 — Outgoing dtn fingerprint

**Batch:** 13 — Outgoing DtN
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage087_outgoing_dtn_fingerprint.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage087_outgoing_dtn_fingerprint_sympy_audit.py`

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

**Notes Derivation Review:** Lambda_2^{out}(z) = z h_2^{(1)}'(z)/h_2^{(1)}(z) is standard l=2 DtN. Static slot -3 = -(l+1) correct. Normalized admittance Y_hat = 1 + z^2/9 + 4z^4/81 + iz^5/27 + ... verified from Hankel expansion. Physical-frequency series through omega^5 correct. Gamma_{5,can}^{DtN} = a^5/(27 c_s^5).

**Script Review:** h_2^{(1)} from jn+I*yn via expand_func. Series through z^7. All even and odd coefficients checked against targets. Omega-domain coefficients verified. All genuine two-sided checks. All pass (exit code 0). Complete coverage.

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The outgoing DtN fingerprint is correct: the `l=2` spherical Hankel DtN operator has static slot `-3`, and the normalized admittance expands to the stated even and odd coefficients through `O(z^7)`. The restored physical-frequency series also matches the claimed `a^5/(27 c_s^5)` odd coefficient.

**Script Review:**

The audit script faithfully builds `h_2^{(1)}` from `j_2 + i y_2`, computes the DtN logarithmic derivative, expands both the DtN operator and the normalized admittance, and checks the frequency-domain coefficients against the notes. The saved output matches the quoted series.

**Issues Found:**

None.

---
