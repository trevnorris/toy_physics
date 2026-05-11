# Review: Stage 043 — Concrete axial overlaps

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (2× MINOR, 2026-04-02)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage026_concrete_axial_overlaps.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`

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
**Verdict:** MINOR

**Notes Derivation Review:**

**1. Equation-level correctness.** All formulas verified:

- u_0(s) = 1/sqrt(L) with int u_0^2 = 1. Correct.
- D/N modes f_n = sqrt(2/L) sin((n+1/2)pi s/L), k_n = (n+1/2)pi/L. Normalization verified.
- Overlap kappa_n = sqrt(2)/((n+1/2)pi). kappa_0 = 2 sqrt(2)/pi ≈ 0.900316. Verified by integration.
- All four overlap assignments (I_{eta,phi}, I_{eta,u}, I_{eta,w}, I_{u,w}) correctly traced through mode definitions.
- All Stage 8 substitutions (C, G_U, G_W, R, Delta, Q, P, N_0, B_0, Z_0, D_0, P_0) correctly computed with kappa factors.
- K_req solution from target equation algebraically correct.

**2. Logical flow.** Clean progression: concrete modes → overlap computation → Stage 8 substitution → K_req solution.

**3. Assumptions.** All stated: flat axial measure, N/N zero mode for wall, D/N half-wave for support/mixed, minimal branch n=0.

**4. Completeness.** Eight undetermined branch parameters clearly enumerated. Hierarchy kappa_n ~ 1/(n+1/2) correct.

**5. Notation consistency.** Consistent with Stage 8 and Stage 7. Minor dual use of u_0 (Neumann zero mode function and brane-like profile) ultimately unambiguous in context.

**6. Physical interpretation.** Sound. Higher D/N levels axially suppressed. K_req as balance condition well-motivated.

**Script Review:**

**B.1 Faithful.** All mode definitions, overlaps, substitutions, and K_req solution implemented.

**B.2 No bugs.** Integration limits, mode expressions, algebraic substitutions all correct.

**B.3 Hardcoded values.** Factor 6 in `K = K_eta + 6 T_Omega` (line 234) used in print statement but not verified by script. Minor gap.

**B.4 No tautologies.** Normalization and overlap checks are genuine SymPy integrals. Branch substitution checks compare SymPy-computed overlaps against hand-derived compact forms.

**B.5 Symbol assumptions correct.** s,L positive; coupling constants real; frequencies positive.

**B.6 Output agreement.** All expect_zero pass (exit code 0). kappa = 2*sqrt(2)/pi matches notes.

**B.7 Coverage.** Mode normalization, overlap law, all substituted quantities, K_req verified. Gap: K = K_eta + 6 T_Omega decomposition not verified.

**Issues Found:**

1. **(MINOR)** The `K = K_eta + 6 T_Omega` wall stiffness decomposition (notes Sec. 5, script line 234-235) is stated without derivation or cross-reference to the specific earlier stage where this is established. The script prints it but does not verify the factor of 6. This feeds into the interpretation of K_req. A cross-reference or brief derivation would improve traceability.

---

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. **Equation-level correctness.** The concrete finite-throat mode choice is consistent: `u_0 = 1/sqrt(L)` is correctly normalized, the D/N ladder `f_n = sqrt(2/L) sin((n+1/2) pi s/L)` satisfies the stated boundary conditions, and the overlap law `kappa_n = sqrt(2)/((n+1/2) pi)` is exact. Substituting the branch overlaps into the Stage 8 scalars produces the quoted `Delta`, `Q`, `P`, `B_0`, `Z_0`, `N_0`, and the branch-level normalization equation without algebraic mismatch.

2. **Logical flow.** The stage is a clean continuation of Stage 8: it replaces formal overlap integrals with a concrete finite-throat family, then propagates those constants into the final normalization test. The hierarchy statement that higher D/N levels are axially suppressed follows directly from the `kappa_n ~ 1/(n+1/2)` scaling.

3. **Assumptions.** The branch assumptions are explicit and sufficient for the calculation: flat axial measure, constant N/N zero mode for the wall/brane-like sector, D/N half-wave support/mixed sector, and the `n = 0` truncation on the minimal branch.

4. **Notation consistency.** The reuse of `u_0` for both the wall profile and the brane-like zero mode is acceptable here because the stage states they are intentionally the same constant mode. The reduced couplings and frequencies carry through Stage 8 consistently.

5. **Physical interpretation.** The branch interpretation is sound. The mode selection makes the constant overlap `kappa = 2 sqrt(2) / pi` the only nontrivial axial constant, and the final `K_req` equation is a sensible balance among support softening, Maxwell/mixed self-load, and outgoing normalization.

**Script Review:**

The script faithfully mirrors the notes and the saved output agrees with every checked identity. The normalization integrals, the general `kappa_n` formula, the branch substitutions, and the exact `K_req` solve are all computed symbolically rather than hardcoded. My independent spot-check of the overlap integral also reproduces `kappa_n` exactly when the integer assumption on `n` is used.

**Issues Found:**

1. The final wall-stiffness identification `K = K_eta + 6 T_Omega` is carried through from prior reduction but is only printed here, not rederived or checked by the script. That does not affect the algebra in this stage, but it leaves a small traceability gap in the branch interpretation.

**Questions:** None.

---
