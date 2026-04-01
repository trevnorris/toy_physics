# Review: Stage 097 — Concrete core schur

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage097_concrete_core_schur.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage097_concrete_core_schur_sympy_audit.py`

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
**Notes Derivation Review:** 2×2 core Schur complement verified. Partial fraction decomposition: rho_c (constant) + sigma_c pole term. Reduced coefficients kappa_c = kappa_0/(1+r_c), gamma_c = gamma_0/(1+r_c) correct — hybridization ratio rescales uniformly.
**Script Review:** M^{-1} computed by SymPy, compared to claimed Schur form. z-expansion identity verified. All genuine, all pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The core Schur reduction is correct. The two-channel inverse yields the stated static term and pole term, and the reduced coefficients `kappa_c` and `gamma_c` are rescaled uniformly by `1/(1+r_c)` as claimed.

**Script Review:**

The script performs the matrix inversion and compares the exact Schur form to the reduced low-frequency form. The saved output matches the displayed coefficient identifications.

**Issues Found:**

None.

---
