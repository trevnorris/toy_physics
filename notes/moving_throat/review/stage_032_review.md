# Review: Stage 032 — Dn overlap zeta

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage032_dn_overlap_zeta.md`
- **Script:** None (status/consolidation stage)

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

**Notes Derivation Review (Status Stage — No Script):**

All claims verified against prior stages and by manual computation:

1. **zeta_n^{phys} = (K_W^{eff}/K_{phi,n}^{eff})(I_n/I_0)^2** correctly derived from Stage 30 definition by substituting coherent coupling pattern.

2. **Overlap integral I_n = sqrt(2L)/((n+1/2)pi)** for uniform source verified by direct integration. cos((n+1/2)pi) = 0 provides the key simplification. Ratio I_n/I_0 = 1/(2n+1) confirmed.

3. **Same-operator twin: zeta_n^{twin} = 1/[(2n+1)^2(1+x n(n+1))]** with x = pi^2 T_X/(L^2 K_W^{eff}) verified. zeta_0^{twin} = 1 confirmed. Higher harmonics suppressed: zeta_n < 1/(2n+1)^2 for n >= 1, x > 0 (denominator exceeds (2n+1)^2).

4. **Notation consistent** with Stages 30-31. D/N basis functions standard.

5. **No claims exceed what is established.** D/N family is new physical input, honestly presented.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review (Status Stage — No Script):**

1. The microscopic replacement for the abstract support ratio is derived correctly. Starting from the coherent-source amplitudes, the stage gets
   `zeta_n^(phys) = (K_W^(eff)/K_(phi,n)^(eff)) (I_n/I_0)^2`
   directly from the Stage 30 definition of `zeta`.
2. For the uniform local source `sigma(s)=1`, I independently checked the D/N overlap integral and obtained
   `I_n = 2 sqrt(2L)/(pi (2n+1))`,
   hence `I_n / I_0 = 1/(2n+1)`, exactly as claimed.
3. The same-operator twin specialization is also consistent:
   `zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ]`.
   This gives `zeta_0^(twin)=1`, while higher harmonics are suppressed both by overlap and by the extra stiffness term.
4. The comparison with the Stage 31 threshold is logically sound. The stage does not overclaim; it turns the support-feasibility question into a branch-by-branch microscopic inequality rather than asserting success automatically.

**Issues Found:** None.

**Questions:** None.

---
