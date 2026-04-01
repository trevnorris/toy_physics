# Review: Stage 040 — Physical parameter map

**Batch:** 7 — Physical Parameter Placement [CP]
**Status:** Verified (2× PASS, 2026-04-03)

**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage040_physical_parameter_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage040_physical_parameter_map_sympy_audit.py`

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

1. **Equation-level correctness.** x-to-kappa substitution: x = pi^2/(kappa+pi^2/4) verified, K_W^{eff} = (T_X/L^2)(kappa+pi^2/4). Robin softening simplifies to A_K = (kappa+pi^2/4)/(kappa+y^2) — verified by expanding denominator. Combined zeta = Omega_Pe^2(kappa+pi^2/4)/(kappa+y^2) correct. Monotonicity: partial_y zeta < 0 (manifestly negative), partial_kappa zeta < 0 for y < pi/2 (numerator y^2-pi^2/4 < 0), partial_Pe zeta > 0 (inherited from Stage 039). Closure ceiling zeta_max = (pi^2/4)(kappa+pi^2/4)/kappa verified by double limit. kappa_max = pi^4/(4(4 zeta_req-pi^2)) verified by solving. All three threshold formulas verified by back-substitution.

2. **Logical flow.** Clean: variable substitution → physical softening → combined family → monotonicity → ceiling → thresholds.

3. **Notation consistent** with Stages 036-039.

**Script Review:**

**B.1-B.7.** Faithful: x-to-kappa substitution, A_K simplification, full zeta formula, both derivatives (sp.diff vs closed form), closure ceiling (sp.limit double limit), kappa ceiling (sp.solve), all three thresholds (sp.solve). No tautologies — each involves genuine CAS computation. All pass (exit code 0). Only unchecked item is Pe monotonicity (inherited from Stage 039 audit).

**Cross-Stage Consistency:**

- **Stage 036:** Omega_Pe formula identical with alpha=Pe substitution. Bound A_I <= pi^2/4 propagated correctly.
- **Stage 037:** A_K formula algebraically equivalent after x = pi^2/(kappa+pi^2/4). Ceiling 4/(4-x) → (kappa+pi^2/4)/kappa verified.
- **Stage 038:** Combined family and closure range [1, pi^2/(4-x)] correctly mapped to physical variables. pi^2/(4-x) = (pi^2/4)(kappa+pi^2/4)/kappa = zeta_max verified.
- **Stage 039:** alpha = Pe identification carried forward correctly.
- **Stage 035:** Reachability condition zeta_req <= zeta_max(kappa) is direct specialization. Three-regime structure preserved.
- **No silently papered-over issues.** Regime boundaries and edge cases correctly maintained.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The physical reparameterization is algebraically correct. I independently checked the key substitution
   `x = pi^2/(kappa + pi^2/4)`,
   and it gives the stated Robin softening factor
   `A_K(eta;kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2)`.
2. The physical lowest-lane family
   `zeta_0^(Pe+R) = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`
   is consistent with Stages 37-39 after `alpha -> Pe` and the `x -> kappa` identification.
3. The monotonicity claims are correct on the constructive branch: the map increases with `Pe` and decreases with both `eta` and `kappa`. The derivative with respect to `kappa` has the expected sign because `0 < y < pi/2`.
4. The closure ceiling and reachability criterion also check out. My independent algebra reproduces
   `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa`
   and
   `kappa_max(zeta_req) = pi^4/[4(4 zeta_req - pi^2)]`.

**Script Review:**

1. The script faithfully verifies the `x -> kappa` substitution, the simplified Robin factor, the physical `zeta_0^(Pe+R)` family, the derivative identities, the closure ceiling, and the three threshold formulas.
2. The saved output matches the notes, and I did not find a tautological check or coding bug.

**Cross-Stage Consistency (Checkpoint):**

1. **Stage 036:** The `Pe` dependence is inherited consistently from `Omega_Pe`.
2. **Stage 037:** The Robin softening factor is exactly the Stage-37 expression rewritten in physical variables.
3. **Stage 038:** The explicit lowest-lane reachability window is preserved, now as a physical ceiling in `kappa`.
4. **Stage 039:** The transport-bias identification `alpha = Pe` is carried through without changing the branch meaning.
5. I did not find any place where Stage 40 changes the earlier threshold logic rather than reparameterizing it.

**Issues Found:** None.

**Questions:** None.

---
