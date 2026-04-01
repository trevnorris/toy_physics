# Review: Stage 091 — Robustness classes

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage091_robustness_classes.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage091_robustness_classes_sympy_audit.py`

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
**Notes Derivation Review:** Class A (pure scale): chi_Q=1 trivially. Class B (scale+argument): even matching forces beta=±1, chi_Q=1. Class C (additive core, beta=1): Sigma2=-Sigma0/9, Sigma4=-Sigma0/27. Preservation submanifold Sigma5 = S(1-beta^5)/9-Sigma0/27 correct. At beta=1: Sigma5=-Sigma0/27.
**Script Review:** Each class independently verified. sp.solve for beta in Class B. Preservation locus checked. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The robustness classification is correct. Pure overall scaling is invisible after normalization, and pure scale-plus-argument deformation cannot move `chi_Q` on the natural positive branch once the canonical even moments are enforced.
2. The additive core channel is the genuine obstruction class: with `beta = 1`, the even-matching conditions reduce to `Sigma2 = -Sigma0/9` and `Sigma4 = -Sigma0/27`, while `chi_Q` can still shift through `Sigma0` and `Sigma5`.
3. The exact preservation submanifold is also correct. At `beta = 1`, it collapses to `Sigma5 = -Sigma0/27`, which is the precise condition for preserving `chi_Q = 1`.

**Script Review:**

The script checks all three classes independently, solves the even-matching condition in the scale-plus-argument branch, and verifies the additive-core preservation locus. I reran the audit locally and it matches the saved output.

**Issues Found:** None.

---
