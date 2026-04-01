# Review: Stage 005 — Grouped p2 normalization bridge

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (2× PASS, 2026-04-02)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage005_grouped_p2_normalization_bridge.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage005_grouped_p2_sympy_audit.py`

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

**1. Equation-level correctness.** All formulas verified:

- **Response coefficients (Sec. 2):** u2 = -D_{A2}/D_{A0} and u4 = (D_{A2}^2 - D_{A0} D_{A4})/D_{A0}^2 from geometric series expansion of D_{A0}/(D_{A0} + D_{A2} omega^2 + D_{A4} omega^4). Verified.

- **Grouped inverse map (Sec. 3):** Standard 5-dimensional real P2 trace/anomaly decomposition with multiplicities (1,2,2). Round-trip identities verified by direct substitution. Anisotropy norm A^2 = 4a_x^2 + (4/5)b_x^2 consistent. Correct.

- **Prefactor Taylor coefficients (Sec. 4):** P0, P2, P4 from series expansion of ratio of polynomials. P4 formula verified by tracing series arithmetic. Correct.

- **Branch coefficients (Sec. 5):** K0=P0, K2=P2+A*P0, K4=P4+A*P2+B*P0, Gamma5=G5*P0 from truncated power series product. Only P0 enters leading odd coefficient. Physically significant and correctly derived.

- **Prototype transfer factor (Sec. 6):** N(omega) = (P0_proto - gW omega^2)^2 / (Delta0 - S2 omega^2 + omega^4)^2 matches Stage 4 after substitutions. N0, N2 verified by direct differentiation. N4 verified through SymPy audit.

- **Normalization algebra (Sec. 8):** mhat^2 P0 = 54 G c_s^5 / (5 a^5 c^5). K2 and K4 target values at mhat=1 verified arithmetically: 54/(9×5) = 6/5 and 54×4/(81×5) = 8/15. All correct.

**2. Logical flow.** Cleanly builds on Stages 3-4: conservative BdG self-energy, Maxwell/mixed sector, then this stage converts operator moments to response moments, defines internal prefactor, lifts to grouped P2 bundle, isolates normalization product.

**3. Assumptions.** All explicit: conservative front-end frozen from Stages 3-4, first-order outgoing expansion, lane-by-lane grouped P2, isotropic branch requires channel-independent coefficients.

**4. Completeness.** Does not overclaim. Correctly identifies remaining task (computing mhat^2 P0 on true moving-throat branch).

**5. Notation consistency.** Consistent with Stages 3-4 and Phase 1 scaffold.

**6. Physical interpretation.** Sound. Leading 2.5PN odd coefficient depends only on static prefactor P0 — clean consequence, not additional assumption.

**Script Review:**

**B.1 Faithful implementation.** All five note sections covered: response moments, prefactor/branch coefficients, grouped P2 inverse map, Stage-4 prototype N0/N2/N4, quadrupole normalization.

**B.2 Correctness of code.** No bugs. `sp.series` with `.removeO()` and `.coeff()` used correctly throughout.

**B.3 Hardcoded values.** Numeric constants (27, 5, 54, 9, 81, 6, 8, 15) are exact rational coefficients verified by the script's own expect_zero checks.

**B.4 Tautological checks.** None found. Series expansions computed from generating functions, not hand-built. Inverse-map round-trip substitutes forward definitions into inverse formulas. All checks have genuine content.

**B.5 Symbol assumptions.** All correct: omega real, D0/D2/D4 nonzero+real, N0/N2/N4 real, Delta0/S2/P0_proto nonzero+real, G/c/c_s/a/mhat positive.

**B.6 Output agreement.** All 18 expect_zero checks pass (exit code 0). Formulas match notes.

**B.7 Coverage.** Comprehensive. Minor gap: Sec. IV.2 dictionary printed but not round-trip tested via expect_zero.

**Issues Found:**

1. **(MINOR)** Section IV.2 dictionary `{Delta0: OA^2*OW^2 - R^2, ...}` is printed but not substituted back into N0/N2/N4 to confirm consistency with Stage 4 forms. Low risk (definitional mapping) but a round-trip check would strengthen the audit.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The algebraic bridge is correct end to end. `u2`/`u4` match the normalized-response expansion, the grouped `P2` inverse map is consistent with the `(1,2,2)` weighting, the outgoing prefactor coefficients `P0/P2/P4` are correct, and the Stage-4 prototype transfer series reproduces the stated `N0/N2/N4` coefficients.
2. The final normalization relation `mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)` and its constant-prefactor consequences for `K2` and `K4` are consistent with the stated `Gamma_5 = P_0 a^5/(27 c_s^5)` bridge.

**Script Review:**

1. The script faithfully covers the five note sections and the cached output matches every `expect_zero` check.
2. The Stage-4 prototype expansion is independently computed from a generating function, not hardcoded, and the `N4` coefficient matches a separate SymPy spot-check.

**Issues Found:**

1. **(MINOR)** The `Omega_A/Omega_W/R/g_A` back-substitution is only printed, not round-tripped through the prototype coefficients. That leaves one small coverage gap in an otherwise clean audit.

**Questions:**

None.

---
