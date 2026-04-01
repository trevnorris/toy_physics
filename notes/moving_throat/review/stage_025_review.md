# Review: Stage 025 — Rank2 selected mode normalization

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage025_rank2_selected_mode_normalization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage025_rank2_selected_mode_sympy_audit.py`

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

1. **Equation-level correctness.** Eigenvector ratio e_1/e_0 verified from both rows of null-space condition (script confirms agreement). Overlap formulas for outgoing/mixed and source channels correct — numerator structures verified by expanding (1 + q*e_1/e_0) type expressions. Generalized normalization F_{q,r,t} four-factor rational form correct. Tracking collapse r=q kills m-dependence, recovers Stage-23. Source-tied specialization with q=sqrt(lam0)*R_U, r=t=sqrt(lam0) verified both by construction and substitution. Flat-U recovery at R_U=1 matches Stage-18/19. First-order deformation H_n^{src} and H_F^{src} verified by series expansion.

2. **Logical flow.** Clean arc from eigenvector → overlaps → normalization → specializations → deformation.

3. **Notation consistent** with Stages 22-24.

**Script Review:**

**B.1-B.7.** All 6 sub-checks faithful. Eigenvector from both matrix rows. Overlaps from eigenvector formula vs closed forms. Source-tied by both construction and substitution. Linearization by sp.series. No tautologies. All pass (exit code 0, 109s runtime from high-degree rational simplification). Complete coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The selected-mode geometry is correct. Using the exact Stage-24 support load, both null-space rows give the same eigenvector ratio, and the closed form `e_1/e_0 = [m(q-r) + r xi] / [delta + xi - m q(q-r)]` matches the symbolic reduction.
2. The generalized normalization function is also correct. The overlapping factors for the loading direction and the source direction combine into the stated three-direction function `F_(q,r,t)(xi,delta;m)`, with the exact denominator `D_(q,r)` controlling the deformation.
3. The stage’s two special reductions are faithful. Tracking support (`r=q`) collapses back to the Stage-23 geometry, and the source-tied split-`U` specialization gives the exact `F_src` with the expected flat-`U` recovery at `R_U = 1`.
4. The first-order deformation formulas in `R_U - 1` are consistent with the exact rational expressions and the natural constructive branch sign choice.

**Script Review:**

1. The script checks both row equations for the eigenvector ratio, builds the generalized normalization law from the exact overlap formulas, verifies the tracking collapse, and then specializes to the source-tied split-`U` case.
2. The saved output matches the note claims, and the linearized deformation checks are nontrivial.

**Issues Found:**

None.

**Questions:**

None.

---
