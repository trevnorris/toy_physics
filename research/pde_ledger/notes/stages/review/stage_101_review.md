# Review: Stage 101 — Parent core extraction

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage101_parent_core_extraction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage101_parent_core_sympy_audit.py`

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
**Notes Derivation Review:** Tanh-wall moments I_f=1/3, I_g=4/15 verified. D/N half-wave data correct. K_s = 6 pi a^2 c_sw hbar/(5 rho_w) under healing lock verified. Bilinear hybridization lambda from GNLS kinetic energy. K_q, g_q, g_s from Maxwell sector. All overlap integrals traced.
**Script Review:** Five sections: shell moments by integration, D/N normalization, K_s (general + healing-locked), bilinear hybridization, mouth couplings. All genuine. All pass (exit code 0). Thorough.
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The parent-action extraction is internally consistent: the tanh-wall moments, D/N normalization, shell stiffness, bilinear hybridization, and mouth couplings all line up with the formulas stated in the note.
2. The healing-locked specialization reproduces the quoted closed form for `K_s`, so the branch reduction is mathematically coherent.
3. The note is careful about scope: it presents a concrete parent-level realization of the core parameters without claiming the full PDE branch is uniquely selected.

**Script Review:**

The script genuinely computes the shell moments, D/N integrals, `K_s`, the bilinear `sq` coefficient, and the mixed couplings. The saved output matches the note and the checks are nontrivial.

**Issues Found:** None.

**Questions:** None.

---
