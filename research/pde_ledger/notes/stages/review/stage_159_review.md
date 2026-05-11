# Review: Stage 159 — Full grouped bundle

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (current PASS after one-port replay, 2026-04-20)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage023_full_grouped_bundle.md`
- **SymPy:** `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
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

**1. Equation-level correctness.**

- **Projector calculus (Sec. 1):** Weighted grouped metric Ggrp = diag(1,2,2) consistent with (1,2,2) degeneracy. Basis vectors e_bar=(1,1,1), e_a=(4,-1,-1), e_b=(0,1,-1) verified orthogonal under Ggrp: e_bar^T Ggrp e_a = 4-2-2 = 0. Projectors are correct rank-1 outer products, idempotent, and complete (P_bar + P_a + P_b = I_3). All verified.

- **Lane Lagrangian (Sec. 2):** Direct multi-port generalization of Stage 021 single-port model. Wall backbone (M_A, K_A), BdG support modes with bilinear coupling, Maxwell/mixed pairs with frequencies, internal mixing R, and wall couplings g_U, g_W. No sign errors.

- **Low-frequency coefficients (Sec. 3):** B_{A,n} moments match Stage 003 convention. Z_{A,n} and N_{A,n} formulas are multi-port generalizations of Stage 021. Notation correspondence verified: Stage 023's (Delta, S, Q, G, P) = Stage 021's (D0, S2, N0, G2, P0_proto).

- **Isotropic branch formulas (Sec. 5):** u_2 = -D_2/D_0 and u_4 = (D_2^2 - D_0 D_4)/D_0^2 match Stage 022. Prefactor formulas P_0, P_2, P_4 match Stage 022 Section 4.

- **Constant-prefactor conditions (Sec. 6):** N_2 and N_4 conditions verified. Substituting N_2 = 2 D_2 N_0/D_0 into N_4 formula yields N_0(2 D_0 D_4 + D_2^2)/D_0^2, matching script output.

- **First-order anisotropy (Sec. 7):** delta u_2 = -(delta D_2 + u_2 delta D_0)/D_0 and delta P_0 = (delta N_0 - P_0 delta D_0)/D_0 are standard quotient-rule expansions. Correct.

- **Monotonicity (Sec. 8):** dP_0/dN_0 > 0 and dP_0/dB_0 > 0 on stable branch (D_0 > 0). Straightforward calculus. Correct.

**2. Logical flow.** Clean progression from Stage 022. Normalization target correctly cited.

**3. Assumptions.** Main assumption (grouped lanes don't mix at linear order on isotropic background) explicitly stated and consistent with Stage 003 angular selection rules.

**4. Completeness.** All six anisotropy channels handled. Edge case D_0 = 0 acknowledged via "stable branch" qualifier.

**5. Notation consistency.** Consistent with Stages 003-022. G_{A,r} notation noted to avoid clash with Newton's G.

**6. Physical interpretation.** Sound: bundle organizes all microscopic inputs, projectors separate isotropy from anisotropy, theorem gap sharpened to one ratio.

**Script Review:**

**B.1 Faithful.** Covers projector calculus (I), bundle assembly (II), isotropic formulas and constant-prefactor conditions (III), anisotropy transport (IV), monotonicity (V).

**B.2 No bugs.** Series expansions correct. `expect_zero` and `grouped_parts` helpers correct.

**B.3 Hardcoded values.** Denominators 5, 10, 2 fixed by (1,2,2) weighting. Target carried symbolically.

**B.4 Tautological checks.** Section II.2 (grouped parts round-trip) is mildly tautological (tests helper function consistency). Critical checks in Sections I (projector properties), III (series identities), IV (epsilon expansion) are substantive and non-tautological.

**B.5 Symbol assumptions.** D0, D2, D4 nonzero+real (appropriate for denominators). K, B0, Z0, N0 positive for stable branch.

**B.6 Output agreement.** All expect_zero pass. N_2 target = 2*D2*N0/D0. N_4 target = N0*(2*D0*D4 + D2^2)/D0^2. Consistent.

**B.7 Coverage.** Does not re-derive Z_{A,n} and N_{A,n} from Lagrangian (inherited from Stage 021 audit). Acceptable.

**Issues Found:**

1. **(MINOR, coverage)** Script does not independently re-derive Z_{A,n} and N_{A,n} low-frequency formulas from the Lagrangian (tested in Stage 021). A single-port spot-check within Stage 023's script would strengthen confidence.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The weighted projector calculus is correct: the `diag(1,2,2)` metric, the `e_bar/e_a/e_b` basis, and the three projectors are all mutually orthogonal and complete, and the grouped vector decomposition is exact.
2. The full coupled bundle formulas, isotropic-branch response moments, constant-prefactor conditions, and first-order anisotropy transport laws all match the reduced algebra in the note. The normalization target is carried through consistently to `mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

**Script Review:**

1. The script is faithful to the note structure and the cached output confirms the projector identities, the full-bundle coefficient decomposition, the isotropic formulas, and the monotonicity derivatives.
2. A direct SymPy spot-check of the `N4` prototype formula agrees with the note and output.

**Issues Found:**

1. **(MINOR)** The stage treats `Z_(A,n)` and `N_(A,n)` as assembled reduced coefficients, but the script does not independently reconstruct them from the underlying Lagrangian in this stage. That is a coverage gap, not a correctness problem.

**Questions:**

None.

---

### Agent: GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

1. The checkpoint claim is still the exact full grouped-bundle packet:
   weighted projector calculus for the real grouped `P_2` lanes, exact
   grouped decomposition of the conservative/output coefficients, isotropic
   branch formulas for `u_2`, `u_4`, `P_0`, `P_2`, `P_4`, the
   constant-prefactor constraints, and first-order anisotropy transport.
2. I did not find a mathematical mismatch between the current note and the
   live CAS outputs. The grouped `(1,2,2)` metric, projector normalization,
   decomposition weights, and the carried Stage-022 normalization target are all
   used consistently.

**Script Review:**

1. The earlier coverage caveat is stale. The current SymPy audit now includes
   a representative one-port reconstruction of both `Z_n` and `N_n` directly
   from the underlying `(\Delta, S, Q, H, P)` lane data before assembling the
   grouped packet.
2. The Mathematica mirror independently replays the same structure: weighted
   projector calculus, representative one-port `Z_n/N_n` reconstruction,
   grouped bundle assembly, isotropic branch formulas, constant-prefactor
   conditions, anisotropy transport laws, and monotonicity derivatives.
3. I did not find a remaining trust defect at the checkpoint level. The open
   gap is only the absence of a dedicated numerical-stress layer for this
   checkpoint family, which is a coverage priority but not a symbolic-fidelity
   failure.

**Issues Found:**

- None. The old “assembled inputs only” caveat is closed in the live audit
  surface.

**Questions:**

- None.

---

### Agent: GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

1. The grouped-bundle theorem is unchanged, but the last remaining carry-forward
   literal at the normalization-product seam is now gone.
2. `Gamma5_port` is rebuilt through the same exact Stage-021/022 outgoing `l=2`
   route used in Stage `022`, so the bundle-level invariant product no longer
   silently depends on a retyped `a^5/(27 c_s^5)` coefficient.

**Script Review:**

1. The SymPy audit now derives `Gamma5_port` from the spherical Hankel
   `l=2` DtN branch before solving the Stage-6 normalization target.
2. The Mathematica mirror now does the same rather than carrying the `27`
   denominator as a local literal.

**Issues Found:**

- None. The Stage-022 carry-forward seam is now anchored in both CAS layers.

**Questions:**

- None.
