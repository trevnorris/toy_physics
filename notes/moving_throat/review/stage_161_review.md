# Review: Stage 161 — Outgoing port coloading theorem

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage161_outgoing_port_coloading_theorem.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage161_outgoing_port_coloading_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** Transfer slope nu_r = 2(p_r-d_r) correct. Weighted collapse Xi_1 = nubar_N - kappa_1. Numerator/detuning slopes with convex decomposition alpha+beta=1. Identity chi_r-zeta_r=1 verified. Combined nu_r formula. Stage 160 equivalence via slippage variables. Co-loading theorem: Xi_1=0 iff nubar_N = kappa_1. Naive rigidity fails (Xi_1=-kappa_1 ≠ 0).
**Script Review:** Series expansion of perturbed ratios. Weighted collapse with 3 ports. Microscopic port slopes via independent computation. Stage 160 bridge verified. All genuine. All pass (exit code 0).
**Issues Found:** None.

---

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** The outgoing-port co-loading theorem is correct. The microscopic slope law `nu_r = 2(p_r - d_r)` follows from the port data, the weighted collapse `Xi_1 = nubar_N - kappa_1` is exact, and the equivalence to the Stage-160 slippage form `nu_r = kappa_1 + sigma_r` is consistent with the note.
**Script Review:** The audit script independently derives the exact slope of `N_{A,0}^{(r)}`, recomputes the weighted collapse, solves the actual `p_r` and `d_r` formulas from the port data, and checks the Stage-160 equivalence. The saved output is clean and matches the note.
**Issues Found:** None.
