# Review: Stage 027 — Continuum selected rank2 closure

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage027_continuum_selected_rank2_closure.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage027_continuum_selected_rank2_sympy_audit.py`

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

1. **Equation-level correctness.** Continuum-selected branch equation correctly substitutes Stage-24 n_req with Stage-22/26 data. Quadratic rearrangement xi^2 + B_cont xi + C_cont = 0 verified. Physical root selection: at M_mix=M_supp=0, roots are 0 and -delta, + branch gives xi=0. Normalization function F_cont correctly specializes Stage-25. Source-tied surface (R_phi=1) and tracking surface (R_phi=R_U) correctly identified from Stage 26. Mismatch penalty lambda_0 M_mix M_supp (R_U-R_phi)^2 in C_cont correct.

2. **Logical flow.** Clean: data collection → quadratic theorem → normalization → special surfaces → mismatch penalty.

**Script Review:** Verifies quadratic form by expanding n_req-M_supp numerator. Special surface collapses verified by substitution. Mismatch penalty by R_phi → R_U-Delta_R. No tautologies. All pass (exit code 0). Complete coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The continuum-selected branch equation is correct. Substituting the actual continuum `M_supp`, `R_phi`, `M_mix`, and `R_U` data into the Stage-24 support theorem produces the exact quadratic
   `xi^2 + B_cont xi + C_cont = 0`,
   and the physical root `xi_phys` is the `+` branch that returns to zero when the loads vanish.
2. I independently checked the two most failure-prone reductions: the quadratic branch equation and the normalization function. Both agree with the note and the saved output, including the source-tied limit `R_phi = 1` and the interference-matched tracking limit `R_phi = R_U`.
3. The mismatch penalty is correctly positive and quadratic in `R_U - R_phi`, so the generic first extended kernel really does define an intermediate closure rather than collapsing to either special surface.
4. The stage is a good checkpoint: it turns the abstract rank-2 geometry from Stages 24-25 into explicit continuum-selected branch equations without changing the meaning of the selected-branch data.

**Script Review:**

1. The audit script faithfully checks the quadratic theorem, the zero-load root, the continuum-selected normalization function, the source-tied collapse, the tracking collapse, and the exact mismatch penalty. The saved output matches all of those claims.
2. I did not find a coding error or a tautological assertion. The script’s factor-9 convention in the quadratic expansion is just a rescaling of the numerator, not a mismatch.

**Issues Found:**

None.

**Questions:**

None.

---
