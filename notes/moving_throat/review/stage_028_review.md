# Review: Stage 028 — Coherent local tracking

**Batch:** 4 — Kernel Continuation
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage028_coherent_local_tracking.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage028_coherent_local_tracking_sympy_audit.py`

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

1. **Equation-level correctness.** Coherent local kernel hypothesis cleanly stated. Coupling amplitudes from modal projection correct. Exact tracking theorem: g_B g_R = g_W g_S manifestly from identical wall/U ratio gamma — algebraically immediate. Common R_tr = [1+chi_0/(1+delta_U)]/(1+chi_0) correct. Range identities 1-R_tr and R_tr-1/(1+delta_U) verified by direct algebra. Total M_tr = M_mix + M_supp correctly factored. Quadratic collapse: mismatch penalty vanishes, C_cont = -delta M_tr, branch equation reduces correctly. G_tr = 9 xi(xi+delta)/[9 delta+(9+2R^2)xi] from lambda_0=2/9 substitution.

2. **Physical interpretation.** Sound: coherent local kernel automatically lands on tracking surface.

**Script Review:** Coupling amplitudes from scratch, g_B g_R - g_W g_S = 0 verified. Range identities, total baseline, branch collapse all checked. No tautologies. All pass (exit code 0). Adequate coverage (normalization collapse already verified in Stage 027 check 27.4).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The coherent-kernel tracking identity is correct. I independently rechecked the local-kernel amplitude relations and confirmed `g_B g_R = g_W g_S`, which forces `rho_0 = sigma_0` and therefore lands the first concrete kernel exactly on the tracking surface.
- The collapse of the Stage-27 quadratic branch equation is also correct. Substituting `R_phi = R_U` removes the mismatch penalty and reduces the two-parameter rank-2 closure to the one-parameter tracking law `M_tr = G_tr(xi,delta;R_tr)`.
- The physical interpretation is sound: this is genuinely more specific than the generic Stage-27 intermediate closure, and the exact common factor `R_tr` really is trapped in the interval `1/(1+delta_U) < R_tr < 1` on the constructive branch.

**Script Review:**

- The script faithfully verifies the tracking identity, the ratio coincidence, the range identities, the total baseline, and the collapse of the Stage-27 branch equation.
- I did not find a code bug or a tautology. The output matches the printed formulas and the symbolic reductions are nontrivial.

**Issues Found:**

- **[MINOR] The stage’s normalization-law claim is not independently checked by the script.** The notes state that the coherent local kernel also gives the exact tracking normalization law `R_target = F_tr(xi_phys, delta; R_tr)`, but the audit script stops after the branch-equation collapse and does not programmatically verify the `F_tr` formula on the coherent branch. The derivation is likely correct, but the coverage is incomplete relative to the note.

---
