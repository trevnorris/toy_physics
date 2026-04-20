# Review: Stage 099 — Dn mixed tube realization

**Batch:** 14 — General DtN & Core Extraction
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage099_dn_mixed_tube_realization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage099_dn_mixed_tube_realization_sympy_audit.py`

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
**Notes Derivation Review:** D/N half-wave Omega_W = pi c_s/(2L_W). kappa_0 = 4L_W^2/(pi^2 a^2) from eigenfrequency. L_W from compensation: L_W = (pi a/2)sqrt((1+r_c)/3). Bare branch D_W^bare = (1+r_c)(1-z^2/3-iz^5/9) gives kappa_c=1/3, gamma_c=1/9 after renormalization. Geometrically concrete.
**Script Review:** kappa_0 from tube length, L_W solved, kappa_c and gamma_c after hybridization, scaled-canonical renormalization. All pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The finite D/N tube realization is correct. The mixed-channel half-wave gives the stated `kappa_0`, the compensation-selected `L_W` follows directly, and the scaled-canonical bare outlet renormalizes to the canonical `kappa_c = 1/3`, `gamma_c = 1/9`.

**Script Review:**

The script computes the tube-length relation, checks the final reduced coefficients, and verifies the scaled-canonical renormalization directly. The saved output matches the notes.

**Issues Found:**

None.

---
