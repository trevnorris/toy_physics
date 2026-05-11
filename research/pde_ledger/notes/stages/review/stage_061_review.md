# Review: Stage 061 — Nonconstant axial family

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage027_nonconstant_axial_family.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py`

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

1. **Equation-level correctness.** All equations verified: N/N basis functions u_0, u_1 correct with proper BCs and normalization. D/N half-wave f_0 correct. Overlap kappa_0 = 2sqrt(2)/pi verified by integration. kappa_1 = -4/(3pi) verified using product-to-sum identity. rho = 2sqrt(22)/(3pi) from sigma = 88/(9pi^2). Blind angle tan(theta_blind) = 3sqrt(2)/2 correct. Max-coupling tan(theta_max) = -sqrt(2)/3 correct with sin^2(theta_max) = 2/11. Wall stiffness K_geo(theta) correctly includes T_w(pi/L)^2 sin^2(theta) from first cosine mode. Section 4 substitutions all algebraically direct from Stage 9.

2. **Logical flow.** Clean: basis → overlaps → wall stiffness → Stage-8/9 substitution → normalization → special cases.

3. **Assumptions.** All explicit: same profile for wall and brane, D/N half-wave for support/mixed, one-parameter coherent family.

4. **Completeness.** Three special cases (Stage-9 recovery, blind angle, max-coupling) all treated.

5. **Notation consistent** with Stages 8-9.

6. **Physical interpretation.** Sound: angular family trades overlap vs stiffness.

**Script Review:**

**B.1-B.7.** Faithful implementation. No bugs. No hardcoded final answers (blind/max substitutions derived from tangent values). No tautologies — overlaps computed by sp.integrate, wall stiffness by integrating chi*G_eta*chi with sp.diff. All expect_zero pass. Complete coverage of all sections. Minor style issue: nested function calls produce repeated output blocks (cosmetic only).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The first nonconstant family is set up correctly. `u_0`, `u_1`, and `f_0` satisfy the quoted N/N and D/N boundary conditions and normalizations, and my independent SymPy check reproduces `kappa_1 = -4/(3 pi)` and the closed overlap law `kappa(theta) = 2(-2 sin theta + 3 sqrt(2) cos theta)/(3 pi)`.
- The blind-angle and max-coupling statements are algebraically consistent: `kappa(theta_blind)=0` gives the stated no-go branch, while `tan(theta_max)=kappa_1/kappa_0=-sqrt(2)/3` and `sin^2(theta_max)=2/11` match the amplitude-maximizing direction.
- The geometry-side lifting is also correct. Evaluating `G_eta = -T_w d_s^2 + K_eta + 6 T_Omega` on `chi_theta` gives `K_geo(theta)=K_eta + 6 T_Omega + T_w pi^2 sin^2(theta)/L^2`, so the tradeoff between stronger D/N overlap and higher axial-gradient cost is derived, not asserted.
- The substitution back into the Stage 8/9 module is faithful: the full branch really is the Stage-9 branch with `kappa_0 -> kappa(theta)` plus the axial-gradient penalty in `K_geo`.

**Script Review:**

- The audit script computes the basis integrals, overlap law, special-angle data, wall-stiffness expectation, branch substitution, Stage 9 recovery, and blind-angle no-go directly with SymPy. The saved output matches the note’s formulas and all checks pass.
- The script is somewhat repetitive because helper sections are re-run inside later sections, but that is only a presentation issue; I did not find a correctness bug or a tautological check.

**Issues Found:**

None.

**Questions:**

None.

---
