# Review: Stage 022 — Split u sector

**Batch:** 4 — Kernel Continuation
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage022_split_u_sector.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage022_split_u_sector_sympy_audit.py`

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

1. **Equation-level correctness.** U-mode splitting K_{U0}=K_U, K_{U1}=K_U(1+delta_U) correct. Split softening delta_split = (delta_0 + eps_eta delta_U/(1+delta_U))/(1-eps_eta) verified by expanding A_1/A_0 - 1. Mixed blocking eps_W_split = eps_W(9 delta_U+11)/(11(1+delta_U)) verified via S_U computation — alternate form 1-(2/11)delta_U/(1+delta_U) confirmed. Direction-splitting invariant D_dir = -kappa_0 kappa_1 g_W rho_0 delta_U/(1+delta_U) verified by direct algebra; collinearity iff delta_U=0 or rho_0=0. Split placement map maintains same structure with eps_W → eps_W_split.

2. **Logical flow.** Clean: U splitting → direct softening → mixed blocking → direction splitting → placement map.

3. **Physical interpretation.** Sound: scalar placement survives but directional collinearity breaks.

**Script Review:** Constructs all quantities from mode data. No tautologies. All pass (exit code 0). Comprehensive coverage including series expansions.

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The split-`U` continuation is consistent with Stage 21. Turning on `T_U` gives the expected internal-mode split `K_(U1) = K_U (1 + delta_U)`, and the direct wall softening then yields the shifted anisotropy
  `delta_split = [delta_0 + eps_eta delta_U/(1 + delta_U)]/(1 - eps_eta)`.
- The mixed blocking correction is also correct. Recomputing the overlap-weighted inverse kernel `S_U` gives the stated
  `eps_W_split = eps_W [1 - (2/11) delta_U/(1 + delta_U)]`,
  so the first axial `U` structure lowers the effective mixed blocking relative to the flat-doublet case.
- I independently checked the load-bearing directional statement:
  `D_dir = kappa_0 z_1 - kappa_1 z_0 = - kappa_0 kappa_1 g_W rho_0 delta_U/(1 + delta_U)`,
  and the difference reduces to zero. That confirms the key theorem that scalar placement survives while source/loading collinearity generically fails once both `delta_U` and `rho_0` are nonzero.
- The split placement map and product law remain algebraically consistent, so the stage cleanly separates “scalar placement” from “directional geometry.”

**Script Review:**

- The script faithfully checks the split wall softening, the mixed blocking ratio, the rotated loading vector, the direction-splitting invariant, the split placement formulas, and the small-`delta_U` expansions.
- The saved output matches the note and I did not find a coding error or a tautological assertion.

**Issues Found:**

None.

**Questions:**

None.

---
