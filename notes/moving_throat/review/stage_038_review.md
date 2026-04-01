# Review: Stage 038 — Explicit lowest lane reachability

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage038_explicit_lowest_lane_reachability.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage038_explicit_reachability_sympy_audit.py`

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

1. **Equation-level correctness.** Combined zeta_0^{exp+R} = Omega_exp(alpha)^2 * A_K(eta) correct product. Closure range: lower (alpha=0, eta=∞) gives 1×1=1; upper (alpha→∞, eta→0) gives (pi/2)^2 × 4/(4-x) = pi^2/(4-x). Reachability criterion zeta_req <= pi^2/(4-x) ⟺ x >= 4-pi^2/zeta_req. Three regimes (A: overlap alone, B: combined, C: fails) correctly classified. Stiffness-ratio form K_X/K_W <= pi^2/(4 zeta_req) verified.

2. **Logical flow.** Clean combination of Stage-36 overlap and Stage-37 softening.

3. **Physical interpretation.** Sound. Independence of alpha and eta in closure range correctly noted as a property of the family, not a claim about the physical system.

**Script Review:** Combined family from Stage-36/37 formulas. Twin and closure limits by sp.limit. Reachability floor by sp.solve. KX/KW equivalence cross-checked. All genuine, all pass (exit code 0).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The explicit reachability family is algebraically correct. The combined branch
  `zeta_0^(exp+R) = Omega_exp(alpha)^2 * A_K(eta)`
  has the stated closure range, with the twin point at `1` and the upper closure `pi^2/(4-x)`.
- I independently rechecked the two load-bearing reachability identities:
  `lim_(alpha -> oo, eta -> 0^+) zeta_0^(exp+R) = pi^2/(4-x)`
  and
  `x_floor = 4 - pi^2/zeta_req`.
  Both reduce to zero exactly. The stiffness-ratio form `K_X/K_W <= pi^2/(4 zeta_req)` is therefore consistent.
- The regime split is physically sensible and matches the notes: overlap alone can cover the pure-overlap ceiling, while the combined exp+Robin family extends the reachability window up to `pi^2/(4-x)`.
- The prior Stage 37 results are carried forward correctly; the current stage is a clean composition of the overlap and compliance branches, not a new derivation trick.

**Script Review:**

- The audit script checks the twin point, the closure maximum, the reachability floor, and the stiffness-ratio equivalence directly from symbolic definitions. The saved output matches the note.
- I did not find a code bug or a tautological check. The script is short but it does verify the actual reachability algebra.

**Issues Found:**

None.

**Questions:**

None.

---
