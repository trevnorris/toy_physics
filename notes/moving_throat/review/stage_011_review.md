# Review: Stage 011 — Loaded profile selection

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage011_loaded_profile_selection.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage011_loaded_profile_selection_sympy_audit.py`

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
**Verdict:** MINOR

**Notes Derivation Review:**

1. **Equation-level correctness.** All verified: Bare wall matrix K_bare = diag(K_0, K_1) with K_1 = K_0 + T_w pi^2/L^2. Loaded matrix K_eff = K_bare - alpha v v^T correct. Trace, determinant, eigenvalue formulas all standard 2×2 results. tan(2 theta) formula verified with correct sign analysis: numerator negative, denominator positive for alpha > 0, giving theta < 0. Weak loading theta ~ alpha kappa_0 kappa_1/DK correct. Strong loading limit matches tan(2 theta_max) = -6sqrt(2)/7 (verified by double-angle formula). alpha_crit = K_0 K_1/(K_1 kappa_0^2 + K_0 kappa_1^2) correct.

2. **Logical flow.** Clean progression from free angle (Stage 10) to dynamically selected angle.

3. **Assumptions.** All explicit: energy functional with alpha > 0, rank-1 structure, two N/N modes.

4. **Completeness.** Eigenvalue spectrum, profile angle, both limits, blind-angle disfavoring, softening threshold all treated.

5. **Notation consistent** with Stage 10.

6. **Physical interpretation.** Sound: rank-1 attractive loading rotates toward max-coupling direction.

**Script Review:**

**B.1-B.7.** Faithful. No bugs. No hardcoded answers. Substantive checks: trace/det from SymPy matrix vs manual forms, characteristic polynomial factorization, stationarity from differentiation, weak/strong limits via sp.series/sp.limit, threshold via sp.solve. All pass.

**Issues Found:**

1. **(MINOR)** Positivity assumption mismatch: Stage 010 declares K_eta, T_Omega as real=True; Stage 011 declares them positive=True. Harmless (both scripts pass) but inconsistent.

2. **(MINOR)** No explicit programmatic assertion that tan(2 theta) < 0 for alpha > 0 (the main physical conclusion). Follows trivially from verified sign data, but an assert would strengthen the audit.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The loaded two-mode model is internally consistent. Starting from `E_eff[q] = (1/2) q^T K_bare q - (alpha/2) (v.q)^2`, the matrix `K_eff = K_bare - alpha v v^T` and the quoted trace, determinant, eigenvalue, and `tan(2 theta)` formulas are the standard exact `2 x 2` results.
- I independently re-derived the angle equation from `dE/dtheta = 0`; it matches the note and script exactly:
  `tan(2 theta) = 2 alpha kappa_0 kappa_1 / (Delta K_ax + alpha (kappa_0^2 - kappa_1^2))`.
  Since `kappa_0 kappa_1 < 0` and the denominator is positive for `alpha > 0`, the selected lower branch indeed rotates in the negative `u_1` direction.
- The weak-loading and strong-loading limits are also correct. The `alpha -> infinity` limit agrees with `tan(2 theta_max)`, so the max-coupling branch is the asymptotic loaded eigenvector, not a hand-picked direction.
- The softening threshold from `det(K_eff)=0` is derived correctly and has the expected interpretation as the exact onset of quadrupole marginality in this reduced model.

**Script Review:**

- The script faithfully checks the matrix formulas, characteristic-factorization identity, stationarity equation, weak/strong-loading limits, and exact softening threshold. The cached output matches the note.
- I did not find a mathematical or coding error. The previously noted sign conclusion is supported by the verified sign ledger even though it is not wrapped in a dedicated boolean assertion.

**Issues Found:**

None.

**Questions:**

None.

---
