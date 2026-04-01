# Review: Stage 050 — Sech gaussian resonance

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage050_sech_gaussian_resonance.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage050_sech_gaussian_sympy_audit.py`

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

1. **Equation-level correctness.** Norms N_ss = 2 w_f, N_pp = w_g sqrt(pi/2) correct. Coherence C^2(r) = I(r)^2/(r sqrt(2pi)) verified. Duality law I(r) = (r/sqrt(pi)) I(pi/r) from Parseval on FT of sech correct. Self-duality gives C^2(r) = C^2(pi/r) — verified by substitution. Stationary point at r_* = sqrt(pi) from chain rule: C^2'(r) = -C^2'(r) → C^2'(r) = 0. Numerical C_res^2 = 0.994418836... consistent with 0.56% deficit.

2. **Logical flow.** Excellent: norms → coherence → duality → symmetry → stationary point → numerical evaluation.

3. **Physical interpretation.** Appropriately cautious: benchmark is real and sharp but does not alone prove threshold survival.

**Script Review:** Mixed symbolic/numeric audit. Duality verified algebraically and numerically (5 sample points at 60-digit mpmath precision, differences < 1e-40). Overlap integral via mp.quad. Monotonicity scan. Stationary point identity. All pass (exit code 0). Strongest script in the batch.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The exact norms, overlap normalization, and coherence formula are consistent: `N_ss = 2 w_f`, `N_pp = w_g sqrt(pi/2)`, and `C^2(r) = I(r)^2 / [r sqrt(2 pi)]`.
2. The duality statement `I(r) = (r/sqrt(pi)) I(pi/r)` correctly implies `C^2(r) = C^2(pi/r)`, so the self-dual point `r = sqrt(pi)` is a valid stationary point. The numerical maximum `C_res^2` and penalty `P_res` match the notes.
3. The interpretation is appropriately constrained: this benchmark nearly saturates the matched branch, but it does not by itself prove threshold survival.

**Script Review:**

The script faithfully implements the notes at the right level of abstraction. It checks the algebraic duality implication, evaluates the non-elementary overlap numerically, and runs a monotonicity scan showing the sampled constructive branch rises to the self-dual point and then falls. Output matches the saved audit file.

**Issues Found:** None.

**Questions:** None.

---
