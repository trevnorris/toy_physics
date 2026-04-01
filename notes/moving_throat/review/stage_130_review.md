# Review: Stage 130 — First order rigidity kernel

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (2× PASS, 2026-04-03)

**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage130_first_order_rigidity_kernel.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage130_first_order_rigidity_kernel_sympy_audit.py`

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
**Notes Derivation Review:** Traction formula and branch law from Stages 125/128. Partial derivatives of Sigma_0 w.r.t. Pi and S_q correct. A_T ~ -4.273 and B_T ~ 0.135 correctly assembled. Kernel representation with centered W_*(x) valid (centering justified by source normalization preservation). |A_T|/B_T ~ 31.68 quantifies overlap dominance.
**Script Review:** All quantities from symbolic formulas. Pi_* via nsolve at 40 digits. A_T, B_T computed. Output matches notes. No hardcoded values, no tautologies. All pass (exit code 0).
**Cross-Stage Consistency:** Pi_*, T_hat_*, Sigma_*, S_* all agree across Stages 125-130 to 17-30 digits. Logical chain from deformation family (129) through traction variation (130) complete. R_q = 1/4 correctly used. No issues papered over.
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The traction variation is assembled correctly from the canonical branch laws `Sigma_0 = Pi / (1 - S_q/4)` and `T_hat_m = sqrt(9 Sigma_0 / 20)`. Differentiating with respect to `Pi` and `S_q` gives the stated linear form.
2. Substituting the Stage-129 retuning relations produces the coefficients `A_T` and `B_T` with the right signs and relative magnitudes. The large ratio `|A_T| / B_T` really does show that overlap perturbations dominate mixed-kernel perturbations at first order.
3. The centered kernel representation `W_*(x)` is justified by preservation of source normalization, so the first-order mouth-shape sensitivity genuinely collapses to a single weighted overlap.

**Script Review:**

The script recomputes the canonical point, the traction coefficients, and the centered rigidity kernel directly from the earlier formulas. The saved audit output matches the note, and my independent check reproduced the canonical `Pi_*` and positive derivative `g'_*` used in the retuning law.

**Issues Found:** None.

---
