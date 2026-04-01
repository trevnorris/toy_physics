# Review: Stage 030 — Coherent kernel map

**Batch:** 5 — Coherent Kernel Map [CP]
**Status:** Verified (2× PASS, 2026-04-03)

**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage030_coherent_kernel_map.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage030_coherent_kernel_map_sympy_audit.py`

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

1. **Equation-level correctness.** Coherent-kernel coupling pattern (shared gamma ratio) correctly forces rho_0 = sigma_0 = chi_0, eps_phi = zeta eps_W, Z_phi = zeta Z_W. All verified by hand-substituting coherent couplings into Stage-021/026 definitions. Support-enhancement factor S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps) algebraically correct — verified M_supp/M_mix ratio. Monotonicity dS/dzeta = (1-eps)/(1-zeta eps)^2 > 0 correct. Product law R_target * M_tr = 8 Lambda(1-eps)/pi^2 * S verified — Z_W, (1+chi_0)^2, (1-eps_eta) cancel. R_target independent of zeta correct.

2. **Logical flow.** Clean: coherent coupling → ratio collapse → split data → enhancement factor → product law → physical interpretation.

3. **Assumptions.** All explicit: coherent local kernel hypothesis, split-U structure from Stage 022, tracking from Stage 028.

4. **Notation consistent** with Stages 021-029. Evolution from rho → chi_0 clearly motivated.

**Script Review:**

**B.1-B.7.** Faithful implementation. Coherent identities from microscopic couplings (not assumed). Proportionality checks use microscopic zeta_def, not abstract zeta. M_tr = M_mix*S check is non-tautological (assembled from different starting points). Product law independently constructed. All pass (exit code 0). Coverage: coherent collapse, proportionality, S factor, product law, R_target independence, dS/dzeta, S(0)=1.

**Cross-Stage Consistency:**

- **Stage 021:** Product law correctly generalized with support enhancement factor S.
- **Stage 022:** delta and eps_W^{split} formulas match exactly.
- **Stage 026:** Tracking condition (sigma_0 = rho_0) explicitly derived as automatic for coherent kernel, not silently assumed.
- **Stage 028:** R_tr formula, M_mix, M_supp baselines all match exactly.
- **Stage 029:** Split-U deficit not contradicted — Stage 030 correctly frames support lane as addressing the follow-up question.
- **Lambda and all effective stiffness definitions consistent across all stages.**
- **No silently papered-over issues.**

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The coherent-kernel compression is algebraically correct. Imposing the shared local source density forces the coherent identities `rho_0 = sigma_0 = chi_0`, `eps_phi = zeta eps_W`, and `Z_phi = zeta Z_W` exactly as stated.
2. The mixed/support baseline split is clean. Rewriting the support contribution in terms of `zeta` and the split blocking `eps` gives
   `M_tr = M_mix + M_supp = M_mix S(zeta;eps)`
   with
   `S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)`.
   I independently checked both the ratio identity and the derivative formula `dS/dzeta = (1-eps)/(1-zeta eps)^2`.
3. The checkpoint product law is also correct:
   `R_target M_tr = 8 Lambda (1-eps)/pi^2 * S(zeta;eps)`.
   This preserves the Stage-21/22 scalar placement structure while isolating the support lane as a pure baseline enhancer.
4. Physically, the stage does the right consolidation work: it shows that coherent support does not alter `R_target` directly, only the total baseline that must overcome the tracking-branch deficit.

**Script Review:**

1. The script faithfully derives the coherent identities from microscopic couplings, checks the support/mixed proportionality, verifies `M_tr = M_mix S`, the product law, `dR_target/dzeta = 0`, the monotonicity of `S`, and `S(0)=1`.
2. The saved output matches the note and I did not find a tautological check or code bug.

**Cross-Stage Consistency (Checkpoint):**

1. **Stage 021:** The original product law survives with the expected coherent support multiplier `S(zeta;eps)`.
2. **Stage 022:** The split quantities `eps = eps_W^(split)` and `delta` are carried over unchanged.
3. **Stages 028-029:** The tracking factor `R_tr` and the coherent local D/N framing are used consistently; Stage 30 is a support-lane continuation of that branch rather than a redefinition of it.
4. I did not find any place where Stage 30 silently drops a prior condition or changes the meaning of the placement variables.

**Issues Found:** None.

**Questions:** None.

---
