# Review: Stage 129 — Positive deformation expansion

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage129_positive_deformation_expansion.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage129_positive_deformation_expansion_sympy_audit.py`

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
**Notes Derivation Review:** Convex deformation Sigma_eps = (1-eps)Sigma_* + eps varsigma preserves positivity and normalization. Two moments (g_bar, S_bar) affine in eps from linearity of integration. Compensation delta_Pi = -eps(g_bar-g_*)/g'_* correct at first order. delta_S formula follows. All six canonical constants match prior stages.
**Script Review:** S_q cross-checked against direct numerical integration at 3 Pi values (strong check). Pi_* via nsolve. Derivatives by sp.diff. All pass (exit code 0). Minor: affine-law tests (lines 74-75) are tautological identities; no end-to-end test with concrete deformation profile.
**Issues Found:**
1. **(MINOR)** Notes call formula "exact" but it is first-order (O(eps^2) dropped).

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The convex deformation family is set up correctly and preserves positivity and normalization exactly, so the linearized shape problem is legitimately parametrized by the two moments `gbar` and `Sbar`.
2. The first-order retuning law `delta Pi = -eps (gbar - g_*) / g'_*` follows directly from holding the overlap fixed at `g_*`, and the induced `delta S_q` correction has the right compensating structure.
3. The numerical canonical data carried into the stage are consistent with the saved audit output, including `Pi_*`, `g_*`, `S_*`, and the derivatives `g'_*`, `S'_*`.

**Script Review:**

The script genuinely cross-checks `S_q(Pi)` against direct numerical integration, solves for the canonical `Pi_*`, and evaluates the needed derivatives and first-order correction formulas. I independently verified the sign and value of `g'_*` at the canonical point. The affine moment checks are tautological, but the substantive kernel computations are real.

**Issues Found:**

1. **(MINOR)** The note calls the correction law "exact" even though it is explicitly a first-order expansion in `eps`.

---
