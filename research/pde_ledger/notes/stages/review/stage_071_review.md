# Review: Stage 071 — Continuum kernel extraction

**Batch:** 3 — Continuum Kernel Extraction [CP]
**Status:** Verified (2× PASS, 2026-04-03)

**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage037_continuum_kernel_extraction.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`

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

1. **Equation-level correctness.** All formulas verified by substitution chain: A = [K_U K_eta^{eff} - c_{etaU}^2]/(mu_eta K_U) from A = K_0 - g_U^2/Omega_U^2 with mass-normalized definitions. DeltaK_ax = pi^2 T_w/(mu_eta L^2) from u_1 eigenvalue. Delta_0, Chi, beta_0, alpha_mix, M_mix all verified through mass-factor cancellations. delta = DeltaK_ax/A verified. Stability inequalities correctly translated to continuum level.

2. **Logical flow.** Clean: continuum Lagrangian → mode basis → mass normalization → Schur complement → closed substitution → stability → summary.

3. **Assumptions.** All explicit: N/N and D/N BCs, brane-like flat-doublet limit for U, conservative static branch, local bilinear couplings.

4. **Completeness.** Both stability surfaces identified as necessary conditions. No missing branches.

5. **Notation consistency.** Continuum couplings c_{etaU} etc. introduced alongside mass-normalized g_* from stages 033-036. Mapping explicit in section 3.

6. **Physical interpretation.** Sound: grounds abstract reduced quantities in explicit continuum operator.

**Script Review:**

**B.1 Faithful.** Mode functions, mass-normalized kernels, 4×4 Schur complement, continuum substitutions all match notes.
**B.2 No bugs.** Sign on g_R in B matrix consistent with Lagrangian coupling convention.
**B.3 No hardcoded values.** Overlap constants verified by integration in script, not assumed.
**B.4 No tautologies.** Schur complement C*B^{-1}*C^T computed by SymPy matrix inversion vs closed-form factorization — cannot pass trivially. Continuum formulas built from abstract definitions vs independently constructed expected forms.
**B.5 Symbol assumptions correct.** Mass/stiffness positive, coupling constants real (either sign).
**B.6 All pass.** Exit code 0. 7 continuum formula checks + 2 stability checks + Schur complement matrix check.
**B.7 Coverage complete.** Mode normalization (9 checks), overlaps (3), Schur factorization (1 matrix), all 7 continuum formulas, both stability inequalities.

**Cross-Stage Consistency (Checkpoint):**

- **Stage 033:** Definitions of A, Delta_0, Chi, beta_0, alpha_0 reproduced identically before continuum substitution.
- **Stage 034:** Softening depth parameterization preserved. Decomposition alpha_0 = g_B^2/varpi^2 + alpha_mix maintained.
- **Stage 035:** delta = DeltaK_ax/A, R_target, M_mix definitions all consistent. M_mix = 8 alpha_mix/(pi^2 A) verified.
- **Stage 036:** Support-feasibility function G and condition M_mix <= G(xi_req, delta) not re-derived but M_mix now given continuum expression. No conflict.
- **Stage 029:** Schur complement factorization re-derived with g_* notation, structurally identical. Sign conventions consistent.
- **No silently papered-over issues.** Stability gates correctly translated without introducing new conditions or dropping old ones.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The continuum reduction is algebraically consistent. Projecting onto `(u_0,u_1,f_0)` and mass-normalizing gives exactly the reduced Stage-12 structure, with `K_0 = (K_eta + 6 T_Omega)/mu_eta`, `DeltaK_ax = pi^2 T_w/(mu_eta L^2)`, `Omega_U^2 = K_U/mu_U`, and `Omega_W^2 = K_W^(eff)/mu_W`.
2. I independently rechecked the two most failure-prone mass-cancellation steps. Both reduce cleanly to zero:
   `A = [K_U (K_eta + 6 T_Omega) - c_(etaU)^2]/(mu_eta K_U)`
   and
   `beta_0 = (mu_W/mu_eta) (K_U c_(etaW) + c_(UW) c_(etaU))^2 / (K_U K_W^(eff) - c_(UW)^2 sigma)^2`.
3. The Schur-complement formulas for `Delta_0`, `Chi`, `alpha_mix`, `M_mix`, and `delta` all follow from the same reduction without sign or factor mismatch. The continuum stability conditions are translated correctly to `K_U K_eta^(eff) > c_(etaU)^2` and `K_U K_W^(eff) > c_(UW)^2 sigma`.
4. Physically, the stage does the right checkpoint work: it turns the abstract reduced inputs from Stages 16-19 into explicit functionals of one concrete finite-throat continuum kernel rather than introducing new phenomenological data.

**Script Review:**

1. The script faithfully implements the note: exact basis integrals, mass-normalized kernels, the `4 x 4` Schur complement, the closed continuum substitutions, and the two stability inequalities.
2. I independently recomputed the continuum Schur factorization and got zero difference from `Sigma_wall = Xi I_2 + alpha v v^T`. The saved output matches the note and all assertions pass.
3. I did not find a coding bug, hidden tautology, or mismatched assumption set.

**Cross-Stage Consistency (Checkpoint):**

1. **Stage 029:** The continuum wall self-energy reduces to the same `Xi I + alpha v v^T` structure, so the stage is genuinely a derivation of the earlier reduced loading law rather than a redefinition.
2. **Stage 033:** The reduced objects `A`, `Delta_0`, `Chi`, `beta_0`, and `alpha_mix` are the same ones used in the microscopic normalization equation, now written in continuum-kernel variables.
3. **Stages 034-036:** The normal-form variables `delta = DeltaK_ax/A` and `M_mix = 8 alpha_mix/(pi^2 A)` are carried over exactly, so the Stage-18/19 admissibility geometry is preserved without any silent modification.
4. I did not find any place where Stage 20 papers over an earlier issue or changes the meaning of the selected-branch gates.

**Issues Found:** None.

**Questions:** None.

---
