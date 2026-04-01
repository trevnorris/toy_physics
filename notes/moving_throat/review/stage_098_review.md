# Review: Stage 098 — Core balance compensation

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage098_core_balance_compensation.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage098_core_balance_compensation_sympy_audit.py`

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
**Notes Derivation Review:** Balance law g_s^2(K_s K_q+lambda^2) = 4(K_s g_q-lambda g_s)^2 verified. g_q solutions correct. sigma_* = g_s^2/(4K_s). Bare channel conditions kappa_0=(1+r_c)/3, gamma_0=(1+r_c)/9. Key result: Y-hat fingerprint preserved for any sigma_* (compensation is sigma_*-independent). Nontrivial.
**Script Review:** Balance law, g_q solve, sigma_*, exact collapse, fingerprint preservation — all genuine symbolic checks. All pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The balance-surface theorem is correct: the codimension-one coupling law produces the exact compensated branch, and the bare-channel conditions `kappa_0 = (1+r_c)/3`, `gamma_0 = (1+r_c)/9` land the core model on the canonical outgoing fingerprint.

**Script Review:**

The script checks the balance equation, solves for the two `g_q` branches, verifies the reduced `sigma_*`, and confirms the exact collapse to the compensated outlet plus the preserved normalized fingerprint. The saved output matches those steps.

**Issues Found:**

None.

---
