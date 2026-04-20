# Review: Stage 031 — Support compensation theorem

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage031_support_compensation_theorem.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage031_support_compensation_sympy_audit.py`

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

1. **Equation-level correctness.** G_tr formula matches Stage 30. dG_tr/dxi numerator manifestly positive (sum of positive terms). M_crit = G_tr(xi=1) = 9(1+delta)/(9 delta+9+2R^2) correct. M_crit - G_tr gap formula confirmed positive for 0 < xi < 1. S(zeta;eps) and dS/dzeta standard. Inverse map zeta_req = (S_req-1)/(1+eps(S_req-2)) verified by substitution. Pole margin and branch margin formulas confirmed.

2. **Logical flow.** Clean: tracking branch critical load → enhancement factor → IVT existence on F_tr → compensation theorem (two cases) → monotone softening response.

3. **Assumptions.** All inherited from Stage 30: 0 < eps < 1, delta > 0, R > 0, xi in (0,1).

4. **Physical interpretation.** Sound: every finite target can be reached before softening. Correctly avoids overclaiming about whether physical PDE produces large enough zeta.

**Script Review:**

**B.1-B.7.** 11 nontrivial symbolic checks: G_tr derivative, F_tr endpoints, M_crit gap, S properties, inverse map, pole/branch margins, implicit derivative. No tautologies. All pass (exit code 0).

**Issues Found:**

1. **(MINOR)** F_tr contains a factor (9 delta+(9+2R)xi) with bare R instead of R^2 — consistent between notes and script, endpoints verified, but asymmetry with R^2 elsewhere worth confirming as intentional.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The compensation theorem is mathematically coherent. `G_tr` is strictly increasing on `0 < xi < 1`, the critical load `M_crit = G_tr(1^-)` is correct, and the exact gap formula `M_crit - G_tr > 0` shows every stable selected-branch point lies below softening.
- The support-enhancement factor behaves exactly as claimed. I checked the inverse map
  `zeta_req = (S_req - 1)/(1 + eps (S_req - 2))`
  and the margin formulas against direct algebra, and they reduce cleanly.
- The existence argument is also sound at the reduced level: `F_tr(0,delta;R)=1`, `F_tr -> +infinity` as `xi -> 1^-`, so every finite `R_target > 1` has at least one stable-side root `xi_req`.
- The final implicit derivative
  `dxi_phys/dzeta = M_mix (dS/dzeta)/(dG_tr/dxi) > 0`
  is correct, so coherent support enhancement always pushes the physical branch deeper into the tracking family.

**Script Review:**

- The script faithfully checks the tracking derivative, the critical-load gap, the enhancement inverse map, the pole and branch margins, and the implicit derivative formula. The saved output matches the note.
- I did find one small assumption gap: the note states `lim_{zeta -> 1/eps^-} S = +infinity`, but the script leaves `eps` only as `positive=True`, so SymPy prints the more generic `-oo*sign(eps - 1)` rather than proving the physical-branch `+infinity` limit directly. The theorem is still correct on the stated window `0 < eps < 1`; the script just does not encode that interval sharply enough.

**Issues Found:**

- **(MINOR)** The audit does not encode the full physical assumption `0 < eps < 1`, so the endpoint divergence of `S(zeta;eps)` is not asserted in the exact form used in the note.

**Questions:**

None.

---
