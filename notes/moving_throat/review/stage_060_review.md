# Review: Stage 060 — Family1 theta extraction

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage060_family1_theta_extraction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage060_family1_theta_extraction_sympy_audit.py`

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

**Notes Derivation Review:** Profile rho_r = [1-10 S^2]_+^{1/4}. Cut point xi_* ~ -0.3856 verified. chi_phi = sech^2/2, I_f = 1/3. Shell moments: <rho>_chi ~ 0.1926, <rho^2>_chi ~ 0.1627. Theta_chi ~ 4.069, Theta_J ~ 0.928. Jensen <rho^2> >= <rho>^2 correctly applied.

**Script Review:** SymPy exact + mpmath 50-digit quadrature. All values match. Jensen check verified numerically.

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The shell-weighted extraction is consistent with the explicit Family-1 wall profile. The cut point, canonical support normalization, and the support-weighted averages of `rho` and `rho^2` all match the stated numerical branch values, and the Jensen floor is correctly ordered below the quadratic average.

**Script Review:**

The script is faithful to the notes: it computes the exact support normalization, uses the explicit cut point for the active shell, and evaluates the branch averages numerically at high precision. The saved output matches those numbers.

**Issues Found:**

None.

---
