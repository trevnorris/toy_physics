# Review: Stage 094 — Mixed sidechannel pole

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage094_mixed_sidechannel_pole.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage094_mixed_sidechannel_pole_sympy_audit.py`

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
**Notes Derivation Review:** Mixed pole Lambda_mix = Lambda_out - sigma_W/(1-kappa_W z^2-i gamma_W z^5). z^2 matching forces kappa_W=-1/9 (sigma-independent). z^4 matching forces sigma_W=0 (even-branch no-go). chi_Q = 3(1-9 sigma gamma)/(3+sigma). Linearization correct.
**Script Review:** kappa and sigma solved sequentially from even matching. chi_Q computed and checked. All genuine (sp.solve). All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The mixed side-channel model expands exactly as stated, and the even matching conditions force `kappa_W = -1/9` and then `sigma_W = 0`, so the standalone isotropic mixed pole cannot preserve the canonical even fingerprint.
2. The outgoing normalization factor `chi_Q^mix = 3(1 - 9 sigma_W gamma_W)/(3 + sigma_W)` is algebraically correct, with the linearized defect `1 - sigma_W(1/3 + 9 gamma_W) + O(sigma_W^2)`.
3. The no-go conclusion is appropriately sharp: if this sector survives on the real branch, it has to be more structured than a naive passive pole.

**Script Review:**

The script solves the even matching equations in the right order, extracts the normalization factor, and checks the linearized form. I independently rechecked the even-branch conditions and they indeed collapse to the no-go result.

**Issues Found:** None.

---
