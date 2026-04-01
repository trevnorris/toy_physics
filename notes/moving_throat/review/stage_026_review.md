# Review: Stage 026 — Support direction extraction

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage026_support_direction_extraction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage026_support_direction_sympy_audit.py`

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

1. **Equation-level correctness.** Continuum extension with c_{Uphi} coupling correctly projected. Effective support vector y = g_B v + g_U g_S D_U v from standard tree-level elimination. R_phi = [1+sigma_0/(1+delta_U)]/(1+sigma_0) correct. Splitting invariant D_phi verified. Overlap contraction v.D_U.v = (sigma/K_U)[1-(2/11)delta_U/(1+delta_U)] verified by expanding kappa_0^2/K_U + kappa_1^2/(K_U(1+delta_U)). Tracking condition g_B g_R = g_W g_S (codimension-1) correct. Mismatch R_phi - R_U = delta_U(rho_0-sigma_0)/[(1+delta_U)(1+rho_0)(1+sigma_0)] verified.

2. **Logical flow.** Clean: continuum coupling → effective vector → direction factor → splitting invariant → tracking condition → mismatch.

3. **Physical interpretation.** Trichotomy (source-tied / intermediate / tracking) cleanly stated and logically follows.

**Script Review:**

**B.1-B.7.** All 7 claims verified: effective vector, R_phi, D_phi, overlap contraction, support baseline, tracking condition (substitutes g_S = g_B g_R/g_W and verifies R_phi-R_U = 0), mismatch formula. No tautologies. Coupling constants correctly real (either sign). All pass (exit code 0, 2s). Complete coverage.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The support-direction extraction is mathematically consistent. The effective support vector `y = g_B v + g_U g_S D_U v` is the correct tree-level result after eliminating the split `U` doublet, and the resulting direction factor
   `R_phi = [1 + sigma_0/(1 + delta_U)] / (1 + sigma_0)`
   matches the note and the saved audit output exactly.
2. I independently rechecked the two core identities with SymPy: the support-direction factor and the splitting invariant `D_phi`. Both reduce to zero against the stated formulas, so the source-tied limit `sigma_0 = 0` and the `delta_U = 0` degenerate limit are both correctly captured.
3. The support pole shift and the `2/11` contraction factor are consistent with the frozen overlap data from Stage 20. The physical support baseline `M_supp` and the tracking condition `g_B g_R = g_W g_S` are therefore derived, not assumed.
4. Physically, the trichotomy is clean: minimal kernel is source-tied, the generic first extension is intermediate, and exact tracking is a codimension-one interference match.

**Script Review:**

1. The audit script faithfully computes the effective support vector, the direction factor, the splitting invariant, the overlap contraction, the support baseline, the tracking condition, and the mismatch formula. The saved output matches those claims.
2. I did not find a bug, hardcoded result, or tautological check that would undermine the verification.

**Issues Found:**

None.

**Questions:**

None.

---
