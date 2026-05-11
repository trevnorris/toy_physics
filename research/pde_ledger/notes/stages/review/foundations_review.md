# Review: Foundations — Roadmap, Phase 0 Spec, Phase 1 Scaffold, Unit Test

**Batch:** 0 — Foundations  
**Status:** Verified (2× PASS, 2026-04-02)

## Files Under Review

| Doc | Path |
|-----|------|
| Roadmap | `notes/moving_throat/moving_throat_pde_roadmap.md` |
| Phase 0 spec | `notes/moving_throat/moving_throat_pde_phase0_spec.md` |
| Phase 1 scaffold | `notes/moving_throat/moving_throat_pde_phase1_linearized_scaffold.md` |
| Unit test | `notes/moving_throat/moving_throat_pde_unit_test_dn_branch.md` |

## Review Focus

- Are the frozen parent equations in phase0_spec self-consistent?
- Does the linearized scaffold correctly specialize the parent equations?
- Is the unit-test benchmark correctly derived from the finite-throat limit?
- Do these documents provide a sound foundation for the 169-stage derivation?

## Agent Reviews

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Findings:**

**Document 1: `moving_throat_pde_roadmap.md`**

A.1 *Equation-level correctness.* The roadmap contains no displayed equations beyond the response operator definition j_i(omega) = sum_j Z_eff,ij(omega) u_j(omega) and the quadrupole target gamma_quad^eff = 2G/(5c^5). Both match the parent hierarchy (confirmed against the 2.5PN summary, which states hat_m_0^2 * Gamma_5 = 2G/(5c^5)).

A.2 *Logical flow.* The six-phase ordering is sound: freeze target (Phase 0), lift geometry (Phase 1), write linearized problem (Phase 2), compute response (Phase 3), match PN ledger (Phase 4), close side conditions (Phase 5), close quadrupole (Phase 6). Each phase depends on its predecessor. The document correctly identifies that the reduced hierarchy leaves the dynamic pole data and quadrupole normalization open, and that these must be derived from a PDE rather than postulated.

A.3 *Assumptions.* The key assumption — that a linearized finite-throat response problem is the correct first target rather than a fully nonlinear PDE — is explicitly stated and well-justified by the existing stack structure. The assumption that the far-field brane reduction must not be imposed at the PDE level is stated clearly.

A.4 *Completeness.* The roadmap deliberately defers dipole and scalar side conditions to Phase 5, which is reasonable for a staged program. No critical branch is missing from the stated scope.

A.5 *Notation consistency.* The acceptance outputs (Omega_{1perp}, Omega_{10}, Omega_0, Omega_{20}, Omega_{21}, Omega_{22}, Omega_g) match exactly the open dynamic pole list in the 2PN summary. The grouped P2 coefficients u_2^(20,21,22) and quadrupole pair (Kbar_0, Kbar_2) are introduced here and carried forward consistently into later documents.

A.6 *Physical interpretation.* The strategy of using the PDE as a generator of response data that must pass back through frozen PN gates is physically sound. The falsification logic in Phase 4 is sharp: if the extracted data miss the gates, the PDE/branch choice is wrong.

**Document 2: `moving_throat_pde_phase0_spec.md`**

A.1 *Equation-level correctness.* The GNLS equation (i hbar D_t psi = [-(hbar^2/2m) D_i D_i + V_conf + h(rho)] psi) matches the parent 4D summary exactly (cf. 4d_summary.md, Section 5.1). The Maxwell equation (partial_M[Z(w) F^{MN}] + (1/xi) partial^N(partial . A) = mu_0(J_psi^N + J_ext^N)) matches the parent (cf. 4d_summary.md, Section 5.3). The geometry equations (M_a a_ddot + Gamma_a a_dot = -partial H_tot/partial a, and similarly for L) match the parent (cf. 4d_summary.md, Section 4.3). The equation of state P(rho) = K_EOS rho^5 is stated; the enthalpy h(rho) = dU/drho is consistent with U(rho) = K rho^5/4 giving h = 5K rho^4/4, which requires K_EOS = K/4 or a similar identification depending on conventions. This is not contradicted by anything in the document but the mapping from K_EOS to K is not made explicit. This is a minor notation point, not an error.

A.2 *Logical flow.* The document progresses cleanly from frozen parent equations (Section 1) to restrictions on what must not be applied (Section 2), to the geometry lift proposal (Section 3), linearization target (Section 4), observables (Section 5), drive/measure protocol (Section 6), acceptance tests (Section 7), and working theorem statement (Section 8). No logical gaps.

A.3 *Assumptions.* All major assumptions are explicit: (i) the parent action is exact at the bulk level, (ii) far-field reductions must not be imposed at the PDE level, (iii) the geometry lift is the first new ingredient not yet frozen. The design requirements for the geometry lift (Section 3) are clearly stated.

A.4 *Completeness.* The document identifies the geometry lift as "not yet frozen" and marks it as the first new ingredient. The continuity equation (partial_t rho + partial_i j^i = 0) is stated alongside the GNLS. The 2PN, 3PN, 2.5PN, and 4PN acceptance gates are listed. No critical items missing.

A.5 *Notation consistency.* Consistent with the roadmap and with the parent 4D summary. The notation D_t, D_i for gauge-covariant derivatives is inherited from the parent summary (cf. 4d_summary.md: D_t psi = partial_t psi + i q_star/hbar A_0 psi, D_i psi = partial_i psi - i q_star/hbar A_i psi). The use of partial_M rather than nabla_M for the Maxwell divergence is consistent with flat-space 4+1 conventions.

A.6 *Physical interpretation.* Correct. The distinction between microscopic ontology (keeping A_w, J^w, F_{mu w}) and brane-effective suppressions is physically important and clearly stated.

**Document 3: `moving_throat_pde_phase1_linearized_scaffold.md`**

A.1 *Equation-level correctness.*

- The level-set parameterization Sigma(X,t) = Sigma_0(X) - eta(Omega,s,t) is clean. The sign convention (eta > 0 means inward wall displacement / throat shrinkage) is consistent with Sigma < 0 being interior.

- The confinement potential linearization (Section 4): V_conf(X;Sigma) = V_conf(X;Sigma_0) + (dV/dSigma)|_0 * (-eta) + O(eta^2), giving delta V_conf = -(dV/dSigma)|_0 * eta. Verified: since Sigma = Sigma_0 - eta, the chain rule gives the correct sign.

- The Bogoliubov coupling C_Sigma (Section 5.1): the upper component +delta_V_conf * psi_0 and lower component -delta_V_conf * psi_0* are correct. The GNLS for delta_psi has +delta_V_conf * psi_0 on the RHS. The conjugate equation, after the standard BdG sign flip, gives -delta_V_conf * psi_0* in the lower component.

- The linearized Maxwell equation (Section 5.2) is the straightforward linearization of the parent equation. Correct.

- The distributed geometry force (Section 5.3): f_Sigma^(psi) = -delta H_psi / delta Sigma |_0. This is consistent with the parent Hellmann-Feynman form F_a^(psi) = -partial H_psi / partial a = -integral rho * partial_a V_conf d^4X.

- The isotropic decomposition formulas (Section 8.2): verified by hand that (u_bar_2, a_2, b_2) form an invertible parameterization of (u_2^(20), u_2^(21), u_2^(22)), and that isotropy (a_2 = b_2 = 0) correctly forces all three to equal u_bar_2. The weight factors (1:2:2)/5 in u_bar_2 correctly reflect the five independent l=2 real harmonics with the degeneracy pattern (1, 2, 2) for (m=0, |m|=1, |m|=2).

- The benchmark mouth operator Z_00^DN = -(omega/c_s)*tan(omega L_0/c_s) in Section 9.1 is correct (verified in detail under Document 4 below).

- The pole ladder omega_n = (pi c_s/L_0)(n + 1/2) is correct.

- The response expansion Z^eff = Z^(0) + i omega Z^(1) + omega^2 Z^(2) + ... (Section 7) uses alternating real/imaginary coefficient convention. The benchmark has only even powers of omega (purely real coefficients), consistent with a lossless system having no odd/dissipative terms. The "danger" terms at i*omega (scalar), i*omega^3 (dipole), i*omega^5 (quadrupole) correspond to the standard radiation-reaction multipole scaling. Physically correct.

A.2 *Logical flow.* The progression from frozen input (Section 1) through geometry lift (Section 2), mode content (Section 3), potential promotion (Section 4), linearized system (Section 5), branch conditions (Section 6), response operator (Section 7), extraction formulas (Section 8), benchmark (Section 9), acceptance tests (Section 10) is clean and each step builds on the previous one.

A.3 *Assumptions.* (i) The level-set representation is presented as a choice, not derived from first principles. This is appropriate — the document acknowledges it is "the minimal honest lift." (ii) The restriction to a stationary finite-throat background is stated. (iii) The compact/passive/outgoing branch conditions are explicitly listed in Section 6.

A.4 *Completeness.* The l=1 (dipole) mode is deliberately omitted from the spherical harmonic expansion of eta. This is justified because the roadmap defers dipole to Phase 5 and the scaffold states "higher multipoles: deferred unless the first extraction fails." The geometry lane eta_g is included but its angular structure is not specified; this is a gap in explicitness (see Issues below).

A.5 *Notation consistency.* Consistent across documents. One point: the geometry sector operator is called G_Sigma (using Sigma subscript) but acts on eta. This is a naming convention, not an inconsistency, since the operator is derived from the Sigma equation of motion.

A.6 *Physical interpretation.* The recovery of collective variables a(t), L(t) as lowest moments of the distributed field eta is physically well-motivated and explicitly demonstrated (Section 2, equations for a(t) and L(t) as integrals of eta). The no-phenomenological-dissipation design choice for G_Sigma is physically sound: dissipation should emerge from the outgoing branch conditions, not be put in by hand.

**Document 4: `moving_throat_pde_unit_test_dn_branch.md`**

A.1 *Equation-level correctness.* All equations verified by hand:

- General solution hat_phi(s) = A cos(ks) + B sin(ks) to hat_phi'' + k^2 hat_phi = 0. Correct.
- Neumann condition hat_phi'(L_0) = 0: gives B = A tan(kL_0). Correct.
- Simplification to hat_phi(s) = hat_phi_m * cos(k(L_0 - s))/cos(kL_0): verified using the addition formula. Correct.
- Mouth derivative: partial_s hat_phi(0) = hat_phi_m * k tan(kL_0). Correct.
- Mouth operator: Z_00^DN = -k tan(kL_0) = -(omega/c_s) tan(omega L_0/c_s). Correct.
- Pole ladder: cos(kL_0) = 0 gives omega_n = (pi c_s/L_0)(n + 1/2). Correct.
- Low-frequency expansion using tan(x) = x + x^3/3 + 2x^5/15 + O(x^7): Z_2 = -L_0/c_s^2, Z_4 = -L_0^3/(3 c_s^4). Correct.

A.2 *Logical flow.* Clean: problem statement, exact solution, operator extraction, pole structure, low-frequency expansion. Each step follows immediately.

A.3 *Assumptions.* The benchmark assumes: (i) frozen wall, (ii) 1D longitudinal wave equation with sound speed c_s, (iii) Dirichlet at mouth and Neumann at cap. All stated explicitly.

A.4 *Completeness.* The benchmark is the simplest possible unit test. The document correctly identifies the next extensions (breathing wall mode, grouped P2 modes, isotropy test).

A.5 *Notation consistency.* k = omega/c_s is used consistently. Z_2, Z_4 naming as coefficients of omega^2, omega^4 is clear.

A.6 *Physical interpretation.* The half-shifted pole ladder is the standard result for a D/N resonator. The absence of odd-power omega terms reflects the lossless, time-reversal-invariant nature of the frozen-wall benchmark. Correct.

**Cross-document consistency checks:**

- The benchmark operator Z_00^DN in Document 4 matches exactly the formula stated in Document 3 (Section 9.1). Consistent.
- The observable list in Document 2 (Section 5) matches Document 1 (Phase 0). Consistent.
- The acceptance gates in Documents 1, 2, and 3 are compatible and progressively more specific.
- The quadrupole target gamma_quad^eff = 2G/(5c^5) in Documents 1 and 3 matches the 2.5PN summary.
- The parent GNLS, Maxwell, and geometry equations in Document 2 match the 4D summary.

**Issues:**

1. **[MINOR] Unspecified angular content of eta_g** — In Document 3, Section 3, the spherical harmonic expansion of eta includes eta_g(s,t) without an accompanying angular basis function. The document says eta_g is "not reducible to the grouped real mouth harmonics alone," but does not specify what angular or functional form it carries. Stage 001 later clarifies this. *Suggested fix:* add a brief clarifying sentence specifying that eta_g is the longitudinal/length-like mode with trivial angular dependence, or specifying its angular kernel.

2. **[MINOR] K_EOS vs K convention not pinned** — In Document 2, the equation of state uses K_EOS while the parent 4D summary uses plain K. The mapping between K_EOS and K is not stated. Not an error but a notation gap. *Suggested fix:* either pin K_EOS = K/4 explicitly or use the same letter K throughout.

3. **[MINOR] G_Sigma operator content not specified** — In Document 3, Section 5.3, G_Sigma appears but its structure (inertia, stiffness, boundary conditions) is left open. Appropriate at the scaffold stage, but the document should note that its specification is among the first tasks when moving beyond the frozen-wall benchmark.

---

### Agent: Claude Opus 4.6 (B) — 2026-04-02
**Verdict:** PASS

**Findings:**

**Document 1: `moving_throat_pde_roadmap.md`**

A.1 *Equation-level correctness* — The roadmap contains no numbered equations per se, only the response operator schema `j_i(omega) = sum_j Z_eff,ij(omega) u_j(omega)` and the low-frequency target expansion `Z_eff ~ Z^(0) + i*omega*Z^(1) + ...`. These are standard definitions, correctly stated. The quadrupole normalization target `gamma_quad^eff = 2G/(5c^5)` matches the 2.5PN summary's canonical value.

A.2 *Logical flow* — The roadmap correctly sequences the work: freeze theorem target (Phase 0) -> geometry lift (Phase 1) -> linearized PDE (Phase 2) -> response operator (Phase 3) -> match to PN ledger (Phase 4) -> scalar/dipole side conditions (Phase 5) -> quadrupole normalization closure (Phase 6). Each phase depends only on its predecessors. The rationale for deferring scalar and dipole channels to Phase 5 (after the quadrupole closure) is sound since those channels have known subtleties ("doom stories") that should not contaminate the main derivation path.

A.3 *Assumptions* — The roadmap explicitly states: (i) the parent 4D action is exact at the bulk level, (ii) the open problems at each PN order, (iii) the restriction to not impose far-field brane-Maxwell reductions prematurely. All assumptions are justified by referencing the existing stack.

A.4 *Completeness* — The roadmap covers all necessary phases. One minor gap: it does not mention what happens if the linearized problem reveals instabilities on the compact branch. The acceptance criterion (Phase 4) is binary (pass/fail), which is correct for a falsification stage, but no contingency is discussed for what alternative branch to try if the first candidate fails.

A.5 *Notation* — Consistent throughout. Uses `Z_eff`, `gamma_quad^eff`, `Omega` for dynamic pole scales, and standard PN nomenclature.

A.6 *Physical interpretation* — Sound. The strategy of targeting the linearized finite-throat response before the full nonlinear problem is physically well-motivated: the PN hierarchy only requires low-frequency data, which the linearized problem can provide.

**Document 2: `moving_throat_pde_phase0_spec.md`**

A.1 *Equation-level correctness* — The GNLS equation `i*hbar*D_t*psi = [-(hbar^2/(2m))*D_i*D_i + V_conf + h(rho)]*psi` matches the parent action's Euler-Lagrange equation (confirmed against `4d_summary.md` Section 5.1). The Maxwell equation matches the parent (confirmed against `4d_summary.md` Section 5.3). The geometry evolution laws match the parent (confirmed against `4d_summary.md` lines 256-260). All equations are correctly transcribed.

A.2 *Logical flow* — Clean progression: frozen parent equations (Sec 1) -> prohibited premature reductions (Sec 2) -> proposed geometry lift (Sec 3) -> linearization target (Sec 4) -> observables list (Sec 5) -> drive/measure protocol (Sec 6) -> acceptance tests (Sec 7) -> working theorem statement (Sec 8). No logical gaps.

A.3 *Assumptions* — The key assumption that `J_ext` is stationary (so `delta J_ext = 0` under linearization) is not explicitly stated. This matters because the phase0 spec includes `J_ext^N` in the parent Maxwell equation, but later the scaffold's linearization drops it. While this is correct (the external current is a fixed background), the assumption should be stated somewhere in this document for completeness.

A.4 *Completeness* — The observables list (Sec 5) covers all seven dynamic pole scales, the three grouped P2 low-frequency coefficients, and the quadrupole normalization pair. Complete for the stated gates. The acceptance tests correctly cover 2PN, 3PN, 2.5PN, and 4PN.

A.5 *Notation* — One minor inconsistency: the EOS is written as `P(rho) = K_EOS*rho^5`, using `K_EOS` rather than the plain `K` used throughout the parent summaries. Cosmetic issue — the subscript "EOS" is more descriptive but breaks notational uniformity with the parent files.

A.6 *Physical interpretation* — The restriction against premature brane-Maxwell reduction (Sec 2) is physically well-motivated: the mixed-sector components `A_w, J^w, F_{mu w}` are needed for the microscopic PDE even though they vanish in the far-field effective theory.

**Document 3: `moving_throat_pde_phase1_linearized_scaffold.md`**

A.1 *Equation-level correctness* —

*Section 2 (Geometry lift):* The level-set parameterization `Sigma = Sigma_0 - eta` with `eta` as normal wall displacement is standard. The moment-recovery formulas for `a(t)` and `L(t)` from eta are structurally correct.

*Section 4 (Confinement linearization):* The linearized confinement perturbation `delta V_conf = -(dV_conf/dSigma)_0 * eta` is correct. The sign arises because `Sigma = Sigma_0 - eta`, so `delta Sigma = -eta`, and by chain rule `delta V_conf = (dV_conf/dSigma)_0 * delta Sigma = -(dV_conf/dSigma)_0 * eta`. Verified.

*Section 5.1 (Matter sector):* The Bogoliubov/Nambu doubling and the wall-drive coupling `C_Sigma[eta] ~ (delta_V_conf*psi_0, -delta_V_conf*psi_0*)^T` have the correct sign structure.

*Section 5.2 (Gauge sector):* The linearized Maxwell equation correctly drops `J_ext^N` (since the external current is a fixed background). The linearized field strength `delta F_{MN} = partial_M delta A_N - partial_N delta A_M` is correct.

*Section 5.3 (Geometry sector):* The distributed Hellmann-Feynman forces `f_Sigma^(psi) = -delta H_psi / delta Sigma |_0` correctly generalize the parent's `F_a = -integral(rho * partial_a V_conf)` to functional-derivative form.

*Section 8.2 (Isotropic decomposition):* Verified algebraically. Given `u_bar = (u_2^(20) + 2*u_2^(21) + 2*u_2^(22))/5`, the transformation is invertible: `u_2^(20) = u_bar + 4*a_2`, `u_2^(21) = u_bar - a_2 + b_2`, `u_2^(22) = u_bar - a_2 - b_2`. Setting `a_2 = b_2 = 0` correctly forces all three to equal `u_bar`. The degeneracy weights (1,2,2) correctly reflect the m-degeneracy of l=2 real spherical harmonics.

*Section 9.1 (D/N benchmark):* The quoted mouth operator `Z_00^DN = -(omega/c_s)*tan(omega*L_0/c_s)` and pole ladder match the unit-test document. Verified.

A.2 *Logical flow* — Each section builds on its predecessors. Clean and complete.

A.3 *Assumptions* — All assumptions are explicitly stated: stationary background, compact/passive/outgoing branch, rotational invariance on isotropic branch.

A.4 *Completeness* — The scaffold correctly identifies that the l=1 channel appears in the "danger" list but is deferred. The eta_g geometry lane is introduced but its angular structure is left under-specified at the scaffold level. Stage 001 later clarifies it as the axisymmetric residual orthogonal to the `a` and `L` moments, which resolves this.

A.5 *Notation* — Consistent with Phase 0 and the roadmap. The scaffold switches to LaTeX notation while Phase 0 uses plaintext, but all symbols carry through correctly.

A.6 *Physical interpretation* — The low-frequency danger channels (scalar at i*omega, dipole at i*omega^3, quadrupole at i*omega^5) correctly identify the leading dissipative powers for each multipole. The benchmark's purely even expansion correctly reflects its conservative character.

**Document 4: `moving_throat_pde_unit_test_dn_branch.md`**

A.1 *Equation-level correctness* — Verified every step:

1. ODE `phi'' + k^2*phi = 0` with `k = omega/c_s` from the wave equation after Fourier transform. Correct.
2. General solution `phi(s) = A*cos(ks) + B*sin(ks)`. Correct.
3. Neumann BC gives `B = A*tan(kL_0)`. Verified.
4. Closed-form `phi(s) = phi_m * cos(k(L_0-s)) / cos(kL_0)`. Verified by expanding with addition formula.
5. Mouth derivative `phi'(0) = phi_m * k * tan(kL_0)`. Verified by differentiation.
6. Mouth operator `Z_00^DN = -k*tan(kL_0) = -(omega/c_s)*tan(omega*L_0/c_s)`. Correct.
7. Pole ladder `omega_n = pi*c_s/L_0 * (n+1/2)`. Correct (D/N resonator).
8. Low-frequency expansion using `tan(x) = x + x^3/3 + 2x^5/15 + O(x^7)`: `Z_2 = -L_0/c_s^2`, `Z_4 = -L_0^3/(3*c_s^4)`. Verified term by term. Correct.

A.2–A.6 — Clean logical flow, all assumptions explicit, complete for its scope, notation consistent, physically correct (half-shifted pole ladder = quarter-wave resonator).

**Cross-Document Consistency Checks:**

1. Phase0 parent equations -> Phase1 linearization: correct perturbative expansion. The GNLS linearization to BdG form is standard. The geometry lift from `(a,L)` to `Sigma` is correctly flagged as the first new ingredient.
2. Phase1 benchmark -> Unit-test document: exact match, no discrepancy.
3. Observables lists: identical across roadmap, phase0_spec, and scaffold.
4. Quadrupole target: `gamma_quad^eff = 2G/(5c^5)` consistent everywhere.
5. Structural restriction (retaining `A_w, J^w, F_{mu w}`): consistently stated across all documents.

**Issues:**

1. **[MINOR] Missing explicit assumption about delta J_ext = 0** — Phase 0 includes `J_ext^N` in the parent Maxwell equation, but Phase 1's linearized Maxwell equation writes only `mu_0 * delta J_psi^N` on the RHS without explicitly noting that `delta J_ext^N = 0` because the external/background current is stationary. Correct physics but the assumption should be stated. *Location:* `phase1_linearized_scaffold.md`, Section 5.2. *Fix:* Add a one-line note: "The external source `J_ext^N` is a fixed background, so `delta J_ext^N = 0`."

2. **[MINOR] K_EOS vs K notation** — Phase 0 writes `P(rho) = K_EOS rho^5`, while all parent summaries use plain `K`. Trivial naming difference but could cause confusion. *Location:* `phase0_spec.md`, Section 1. *Fix:* Either use `K` consistently or add a note identifying `K_EOS = K`.

3. **[MINOR] eta_g angular structure under-specified** — The scaffold introduces `eta_g(s,t)` as "the geometry lane not reducible to the grouped real mouth harmonics alone" but does not specify its angular dependence. Stage 001 later clarifies it as axisymmetric and orthogonal to the `a` and `L` moments. *Location:* `phase1_linearized_scaffold.md`, Section 3. *Fix:* Add a note indicating that `eta_g` is axisymmetric (l=0 residual orthogonal to the breathing and length moments).

4. **[MINOR] Potential subscript clash for Z_2, Z_4** — The unit-test document defines `Z_2` and `Z_4` as the omega^2 and omega^4 coefficients of `Z_00^DN`, while the scaffold uses `Z^(n)_{AB}` for the n-th coefficient of the response matrix. These are in different contexts and unlikely to cause confusion in practice, but could clash in a future stage. *Location:* `unit_test_dn_branch.md`, Section 5. *Suggestion:* Consider labeling as `Z_{00,2}^{DN}` and `Z_{00,4}^{DN}` for explicit disambiguation.

---
