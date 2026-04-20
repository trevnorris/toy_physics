# Review: Stage 034 — Lowest twin criterion

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage034_lowest_twin_criterion.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage034_lowest_twin_criterion_sympy_audit.py`

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

1. **Equation-level correctness.** Pi_tr = F_tr * G_tr closed form verified by algebraic cancellation. Endpoints Pi_tr(0)=0, Pi_tr→∞ as xi→1 correct. Elimination: zeta_req <= 1 ⟺ Pi_tr <= 2 C_mix via S_req <= 2, factor of 2 maps correctly. Three threshold scales (Lambda_twin,req, M_mix^{twin,req}, Z_W^{twin,req}) all correct rearrangements. Twin-saturation depth xi_{2x} from quadratic 9xi^2 + [...] xi - 18 M_mix delta = 0, discriminant 4×9×18=648 verified.

2. **Logical flow.** Clean: product formula → elimination of zeta_req → threshold scales → saturation depth.

3. **Physical interpretation.** "Doubling criterion" — twin sufficient iff branch product <= 2× mixed product scale. Clear.

**Script Review:** Product from sp.factor(Ftr*Gtr) vs closed form. Anchor values zeta_req(C_mix)=0, zeta_req(2C_mix)=1. Z_W round-trip through Stage-30 map. Quadratic root substituted back into G_tr. All genuine, all pass (exit code 0).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The branch-product algebra is correct. I independently checked that `Pi_tr = F_tr * G_tr` reduces to the closed form shown in the note, and that the stable-branch endpoints are `Pi_tr(0)=0` and `Pi_tr -> +infinity` as `xi -> 1^-`.
2. The elimination of `zeta_req` is also exact. Using `S_req = Pi_tr / C_mix`, the boundary values `zeta_req(C_mix)=0` and `zeta_req(2 C_mix)=1` are correct, and the lowest-twin criterion really is equivalent to `Pi_tr <= 2 C_mix`.
3. The threshold scales `Lambda_twin,req`, `M_mix^(twin,req)`, and `Z_W^(twin,req)` are just clean rearrangements of the same inequality, not new assumptions. The quadratic root `xi_(2x)` also satisfies `G_tr(xi_(2x)) = 2 M_mix`.
4. Physically, the stage is doing the intended consolidation work: it turns the twin-sufficiency question into one branch-product inequality on the physical tracking branch.

**Script Review:**

1. The script faithfully verifies the product formula, endpoint limits, `zeta_req` elimination, threshold scales, and the quadratic saturation depth.
2. The saved output matches the note, and the checks are genuinely computed rather than hardcoded.

**Issues Found:**

None.

**Questions:**

None.

---
