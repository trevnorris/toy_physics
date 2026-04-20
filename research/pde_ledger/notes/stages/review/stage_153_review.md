# Review: Stage 153 — Linear grouped outlet map

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage153_linear_grouped_outlet_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage153_linear_grouped_outlet_map_sympy_audit.py`

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
**Verdict:** MINOR
**Notes Derivation Review:** Linearizations of u_2, u_4, P_0 around canonical branch verified. Outlet map inversion correct. Even-consistency delta u_4 = (8/9) delta u_2 verified. Compatibility constraint delta D_4 = (2/3) delta D_2 + (1/27) delta D_0 correct. Grouped trace/anomaly formulas follow by linearity. Full linear problem collapses to pair (K_A, G_A). All algebra verified.
**Script Review:** Series expansions genuine. Outlet inversion via sp.solve. Even-consistency checked.
**Issues Found:**
1. **(MINOR)** Script lines 110-111: tautological assertions — compare expression to itself for trace/anomaly even-consistency check. Genuine verification exists elsewhere in script for the scalar case.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The grouped-lane linear transport is correct: the Stage-6 expansions for `delta u_2`, `delta u_4`, and `delta P_0` lead directly to the two microscopic obstruction combinations `K_A` and `G_A`.
2. The direct outlet inversion is also correct. Using the Stage-142 hybrid-outlet formulas with `delta C = 0` yields the stated maps from `delta u_2` and `delta P_0 / P_0` into `delta kappa_W` and `delta gamma_W`.
3. The one-parameter even-consistency relation `delta u_4 = (8/9) delta u_2`, equivalently `delta D_4 = (2/3) delta D_2 + delta D_0 / 27`, is the right compatibility condition for a single hidden even-pole deformation.

**Script Review:**

The script genuinely verifies the grouped conservative/output transport formulas, the outlet inversion, and the even-pole compatibility relation. I reran the audit locally and it matches the saved output. The only weak spot is the final trace/anomaly “checks,” which compare expressions to themselves and therefore add no independent verification.

**Issues Found:**

1. **(MINOR)** The `aD4 relation` and `bD4 relation` checks at the end of the script are tautological self-equalities rather than real tests.

---
