# Review: Stage 092 — Linearized branch selection

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage092_linearized_branch_selection.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage092_linearized_branch_selection_sympy_audit.py`

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
**Notes Derivation Review:** Linearization S=1+eps s, beta=1+eps b, etc. chi_Q = 1+eps(5b+a0/3+9a5)+O(eps^2). Scale invariance (s drops out) confirmed. Preservation: 5b+a0/3+9a5=0 → a5=-5b/9-a0/27.
**Script Review:** chi_Q from exact formula, sp.series in eps. Coefficient checked against notes. Scale independence verified. Preservation solved and checked. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The linearization is correct: expanding `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)` about `S=1`, `beta=1`, `Sigma_0=0`, and `Sigma_5=0` gives `1 + eps(5b + a_0/3 + 9a_5) + O(eps^2)`.
2. The overall scaling parameter `s` cancels to first order, so the branch-selection data really are the three deformation scalars `(b, a_0, a_5)`.
3. The preservation condition `5b + a_0/3 + 9a_5 = 0` is the exact first-order constraint stated in the notes.

**Script Review:**

The script computes the exact series expansion in `eps`, extracts the first-order defect coefficient, and solves the linear preservation condition. I independently rechecked the coefficient algebra and it matches the note exactly.

**Issues Found:** None.

---
