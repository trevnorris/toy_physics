# Review: Stage 104 — Support feasibility frontier

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (current PASS after full frontier closure, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage036_support_feasibility_frontier.md`
- **SymPy:** `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
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

1. **Equation-level correctness.** All equations correct. G(xi,delta) = 9 xi(xi+delta)/(9 delta + 11 xi) verified by direct algebra from Stage 18 alpha_req. Derivative dG/dxi confirmed by quotient rule; numerator quadratic has discriminant -72 < 0, confirming strict monotonicity. G(0,delta) = 0, G_max = 9(1+delta)/(9 delta+11) by direct evaluation. Loading split g_B,req^2/varpi^2 = (pi^2 A/8)(G - M_mix) algebraically exact. Near-onset series G = xi - 2 xi^2/(9 delta) + O(xi^3) confirmed by Taylor expansion.

2. **Logical flow.** Clean: G defined from Stage 18 → monotonicity/range → admissible region by combining F and G → near-onset asymptotics → theorem gate. No gaps.

3. **Assumptions.** All inherited from Stages 17-18 (positive delta, xi in [0,1), rank-1 selected-branch model). No new assumptions.

4. **Completeness.** Both endpoints covered, near-onset asymptotic provided. No missing edge cases.

5. **Notation consistency.** Consistent with Stages 17-18. New symbols G and M_mix cleanly introduced.

6. **Physical interpretation.** Coherent: G measures directional loading capacity; G >= M_mix ensures BdG sector doesn't need unphysical negative loading. Parametric frontier in (R_target, M_mix) space is natural.

**Script Review:**

**B.1 Faithful.** All three quantitative sections implemented: G definition and loading split, monotonicity/derivative/endpoints, near-onset series.

**B.2 No bugs.** Standard SymPy calls throughout.

**B.3 Hardcoded values.** Numbers (9, 11, 8, 2, 18, 81) all derive from D/N overlap constants established in Stage 18. Structural, not arbitrary.

**B.4 No tautologies.** Loading split checks two independent algebraic paths. Derivative compared against separately defined target. Series check compares SymPy expansion against manual target.

**B.5 Symbol assumptions.** All positive+real. Correct for stable branch parameters.

**B.6 Output agreement.** All expect_zero pass. Exit code 0. Matches notes.

**B.7 Coverage.** All quantitative claims covered. Qualitative/geometric statements (admissible region, theorem gate summary) reasonably not script-tested.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The support-feasibility split is correct. Recomputing `G(xi,delta)` from Stage 18 gives the same `9 xi(xi+delta)/(9 delta + 11 xi)` expression, and the factorization through `G - M_mix` matches the note exactly.
2. The monotonicity and endpoint claims are sound. I independently checked `dG/dxi`, `G(0,delta)`, and `G_max(delta)`, and they agree with the note and saved output.
3. The near-onset series is also correct, so the geometric interpretation of the combined admissible region in `(R_target, M_mix)` space is internally consistent.

**Script Review:**

1. The script faithfully verifies the support-loading factorization, the derivative form, the endpoint values, and the near-onset expansion. The saved output matches the printed formulas.
2. I ran a small independent SymPy check of the derivative identity and the near-onset series; both reduced to zero.

**Issues Found:**

- **[MINOR] The branch parameter `xi` is again declared `positive=True` even though the note treats `xi = 0` as the onset boundary.** This does not change the algebra, but it is the same boundary-assumption mismatch as Stage 18 and should be made explicit or loosened in the script.

**Questions:**

None.

---

### Agent: GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

1. The checkpoint claim is the exact support-feasibility frontier inside the
   selected D/N normal form: closed-form `G(xi,delta)`, the exact loading split
   through `G - M_mix`, strict monotonicity on `0 <= xi < 1`, the endpoint
   values, and the near-onset asymptotic. I did not find a mathematical
   mismatch between the current note and the live CAS outputs.
2. The note already treats `xi = 0` as the onset boundary and writes the
   parametric frontier with `xi in [0,1)`, so the theorem statement itself is
   internally aligned.

**Script Review:**

1. The earlier boundary-assumption caveat is stale. The current SymPy audit
   declares `xi` as `nonnegative=True`, and the Mathematica mirror assumes
   `0 <= xi < 1` explicitly.
2. Both CAS layers now check the same live surface: exact `G`, loading split,
   manifestly positive derivative form, `G(0,delta)`, the `xi -> 1^-` endpoint,
   and the near-onset series.
3. I did not find a remaining checkpoint-level trust defect. The stage is a
   selected-branch normal-form theorem, and the current scripts match that
   scope cleanly.

**Issues Found:**

- None. The older `xi > 0` boundary mismatch is closed in the live audit
  surface.

**Questions:**

- None.

---

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The stage-level frontier claim is now fully represented in the executable
surface. The scripts no longer stop at the `G`-only loading split; they now
cover the parametric frontier together with the final admissibility gate.

**Script Review:**

1. Both CAS layers now import the carried
   `F(\xi,\delta)` frontier from Stage `035`, define
   `R_target = N_Q A / (\beta_0 (8/\pi^2))`, and verify the exact
   support-loading equivalence already used for `G`.
2. Both scripts now assemble the full three-conjunct admissibility surface:
   `R_target >= 1`,
   `F(\xi_req,\delta) = R_target`,
   and
   `M_mix <= G(\xi_req,\delta)`.
3. The positive/negative witness cases make this materially stronger than the
   previous `G`-only checkpoint.

**Issues Found:**

- None.

---
