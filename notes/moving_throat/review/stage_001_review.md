# Review: Stage 001 — Geometry lift

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage001_geometry_lift.md`
- **Script:** None (status/consolidation stage)

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

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Notes Derivation Review:**

**A.1 Equation-level correctness**

All displayed equations check out:

- **Confinement linearization (Section 3):** The chain rule `delta V_conf = -(V_wall'(Sigma_0/ell_c)/ell_c) delta R` is correct. Since `Sigma = r - R`, we have `delta Sigma = -delta R` (at fixed spatial point), so `delta V = V_wall'(Sigma_0/ell_c) * delta Sigma / ell_c = -V_wall'(Sigma_0/ell_c) * delta R / ell_c`. Signs and factors confirmed.

- **Spherical harmonic normalization (Section 4):** `Y_00 = 1/(2 sqrt(pi))` is the standard real normalized harmonic. The claimed integral `(1/4pi) int_{S^2} Y_00 dOmega = 1/(2 sqrt(pi))` is correct. The extraction `q_00 = 2 sqrt(pi) delta_a` follows correctly. Confirmed by SymPy audit Section I.2.

- **Spherical Laplacian eigenvalue (Section 5):** `-Delta_{S^2} Y_{lm} = l(l+1) Y_{lm}` is standard, used correctly to get modal effective stiffness `K_eta + l(l+1) T_Omega`. For l=2 this gives `K_eta + 6 T_Omega`. Confirmed by SymPy audit Section III.

- **Modal geometry PDE (Section 5 and Section 9):** The Euler-Lagrange equation from the quadratic action correctly yields `mu_eta q_{lm,tt} - d_w(T_w d_w q_{lm}) + [K_eta + l(l+1) T_Omega] q_{lm} = source`. Verified by SymPy audit (Section III, exit code 0). The note uses the product-rule form `d_w(T_w d_w q)` which is correct for varying `T_w`; the audit uses constant `T_w` as a special case.

- **Linearized BdG matter sector (Section 7):** The two-component form is standard Bogoliubov-de Gennes. The relative signs between the `delta psi` and `delta psi*` rows are correct and consistent with Phase 1 scaffold Section 5.1.

- **Linearized Maxwell (Section 8):** Directly transcribed from Phase 0 spec with `delta` perturbation prefix. No index or sign issues.

- **Geometry back-reaction source (Section 9):** `S_eta^(psi) ~ -(V_wall'/ell_c) delta rho` is the correct Hellmann-Feynman structure. This closes the two-way coupling loop correctly.

**A.2 Logical flow** — Clean progression from motivation through reference throat, confinement promotion, mode decomposition, new geometry action, background fields, linearized sectors, mode content, response operator, and observable extraction. No skipped intermediate steps.

**A.3 Assumptions** — All explicitly stated: hybrid level-set/shape-field representation as a choice, quadratic wall action flagged as "new and not yet frozen," brane-Maxwell non-suppression carried from Phase 0, axisymmetric reference throat as a background choice.

**A.4 Completeness** — Higher multipoles (l >= 3) correctly deferred. Dissipation/odd parts correctly deferred to emerge from coupled problem. Bottom closure condition left open, appropriate for skeleton stage.

**A.5 Notation consistency** — Phase 1 scaffold uses `s`, Stage 001 uses `w` (bulk axial direction); identification made explicit. Level-set conventions match: `delta Sigma = -delta R = -eta`. Port basis notation matches Phase 1 scaffold.

**A.6 Physical interpretation** — Sound. Geometry lift from `(a,L)` to distributed field is minimal upgrade to expose grouped P2 modes. Isotropic reference gives degenerate P2 modes. Dynamic pole scales emerge as poles of mode-resolved DtN operator.

**Issues Found:**

None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The stage does the right structural job: the level-set lift `Sigma = r - R` preserves the bulk coupling, the `Y_00` normalization and `q_00 = 2 sqrt(pi) delta a` bridge are correct, and the grouped real `l=2` split is the cleanest way to make the conservative `P2` lane literal.
- I independently re-checked two claims with focused SymPy computations. First, the confinement-chain-rule sign is correct. Second, the variable-coefficient Euler-Lagrange equation from the wall action gives the stated `-partial_w(T_w partial_w q)` form when no extra surface-weight factor is present.
- The linearized BdG and Maxwell schematic equations are consistent with the Phase 1 scaffold: the wall perturbation enters the matter sector through `delta V_conf`, and the mixed `w`-channel gauge data are correctly kept alive rather than projected out too early.
- The mode-level interpretation is sound: isotropic background implies an unsplit geometry-only `l=2` block, while later matter/gauge couplings are the natural source of splitting and odd/passive pieces.

**Issues Found:**

- **[MINOR] Surface-measure handling is under-specified.** Section 5 writes the action with an explicit `sqrt(gamma_0)` factor, but Sections 5 and 9 then use the unweighted modal operator `mu_eta q_tt - partial_w(T_w partial_w q) + ...`. An independent SymPy check with a generic weight `g(w)` gives an extra term `T_w(w) g'(w) q_w / g(w)`, so the stated PDE is correct only if `sqrt(gamma_0)` has already been absorbed into the coefficient functions or the `w` coordinate has been chosen so that the weight is constant. The note should say that explicitly.

---
