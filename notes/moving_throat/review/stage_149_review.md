# Review: Stage 149 — Bundle inversion four drifts

**Batch:** 18 — Linear Defect Transport & Final
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage149_bundle_inversion_four_drifts.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage149_bundle_inversion_four_drifts_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** 4×4 forward transport laws verified: dTheta=2drho, dKs=2da+drho, dKq=dZ+2dcs-2da, dP=5(dcs-da). Inversion by forward substitution: drho=dTheta/2, da=dKs/2-dTheta/4, dcs=dKs/2-dTheta/4+dP/5, dZ=dKq-2dP/5. All correct. Bundle form with dP=dN0-dD0 correct. Frozen-wall (dTheta=0) simplifies correctly. rho_w^{chi} = sqrt(4.069/25) ~ 0.4034 verified.
**Script Review:** sp.solve for 4×4 system. Forward verification, bundle form, frozen-wall all compared against independently typed targets. Numerical rho_w from Stage 60 datum. All pass (exit code 0).
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The bundle inversion is correct. The four observable laws invert exactly to the four remaining drifts, the bundle form with `P_0 = N_0/D_0` matches, and the frozen-wall corollary reduces the remaining freedom exactly as claimed.

**Script Review:**

The script solves the 4×4 linear system, verifies the forward substitutions, checks the bundle identities, and confirms the frozen-wall simplification. The saved output matches the note.

**Issues Found:**

None.

---
