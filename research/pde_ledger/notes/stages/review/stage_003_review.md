# Review: Stage 003 — Bdg coupling

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (current PASS with dual-CAS and shared-stress coverage, 2026-04-20)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage003_bdg_coupling.md`
- **SymPy:** `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- **Numerical stress:** `scripts/numerical/stage003_004_foundational_stress.py`
- **Numerical stress (Mathematica):** `mathematica/numerical/stage003_004_foundational_stress.wl`

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
-->

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Notes Derivation Review:**

**1. Equation-level correctness.** All displayed equations are algebraically correct:

- **Reduced Lagrangian (Sec. 3.1):** Standard quadratic form for wall coordinates Q^A coupled bilinearly to harmonic oscillator matter modes X_alpha. Signs, factors, and conventions consistent with Stage 002.
- **Frequency-space elimination (Sec. 4.1):** `D^eff_0(omega) = K_0 - omega^2 M_0 - C(Omega_0^2 - omega^2 I)^{-1} C^T`. Signs correct: self-energy enters with minus sign. Low-frequency expansion verified: `K_eff = K - C Omega^{-2} C^T` (stiffness lowered), `M_eff = M + C Omega^{-4} C^T` (inertia increased), `N_eff = C Omega^{-6} C^T`. All correct.
- **Grouped P2 kernels (Sec. 4.2):** `D_A^eff(omega) = K_A - M_A omega^2 - sum g_{A,alpha}^2/(varpi_{A,alpha}^2 - omega^2)` with correct low-frequency coefficients d_{0A}, d_{2A}, d_{4A}.
- **Two-pole formula (Sec. 5):** Dispersion `(K - M omega^2)(varpi^2 - omega^2) - g^2 = 0` gives discriminant `(Omega_eta^2 - varpi^2)^2 + 4g^2/M`. Verified. Perturbative shift formulas are standard level-repulsion results.
- **Trace/anomaly decomposition (Sec. 6):** Weights (1, 2, 2) for real l=2 harmonics. Invertible: d_20 = d2bar + 4a_2, d_21 = d2bar - a_2 + b_2, d_22 = d2bar - a_2 - b_2. Consistent with Phase 1 scaffold Sec. 8.2.

**2. Logical flow.** Follows cleanly from Stages 001-002. Stage 002 Sec. 5 outlined exactly this next step (coupling to BdG matter modes).

**3. Assumptions.** All explicit: stable normal-mode reduction (flagged as "controlled reduced-sector move"), isotropic reference throat, weak coupling for perturbative shift (varpi^2 > Omega_eta^2).

**4. Completeness.** Degenerate case varpi^2 = Omega_eta^2 not discussed (shift becomes O(g) rather than O(g^2)), but the perturbative formula explicitly assumes separation. Minor observation, not a deficiency.

**5. Notation consistency.** Q^A = (delta a, delta L) matches Stage 002. M_0, K_0 match axisymmetric inertia/stiffness. Y_{2m}^{real} and grouped labeling consistent.

**6. Physical interpretation.** Sound: matter support lowers stiffness and increases inertia (level repulsion). Isotropy theorem (a_2 = b_2 = 0) correct for identical channels.

**Script Review:**

**B.1 Faithful implementation.** All five sections of the notes covered: E-L reduction, effective kernel, low-frequency renormalization, two-pole spectrum, perturbative shift, channelwise P2, trace/anomaly, isotropy, selection rule.

**B.2 Correctness of code.** No bugs. `euler_equations` used correctly, `expect_zero` helper properly simplifies.

**B.3 Hardcoded values.** None problematic. All parameters symbolic. Spherical harmonics use standard analytic forms.

**B.4 Tautological checks.** None found. E-L checks test Lagrangian construction against independently written expected equations. D_eff check tests matrix formula against component form. Series match tests expansion formulas. Exact roots use sp.solve. Perturbative shift tests perturbative formula against exact solution. Selection rule performs actual angular integrals.

**B.5 Symbol assumptions.** All physically appropriate: wa, wb positive (stable frequencies), coupling constants real, delta positive (enforces assumed ordering), g20/g21/g22 real.

**B.6 Output agreement.** All expect_zero checks pass (exit code 0). All output matches notes.

**B.7 Coverage.** Complete: E-L reduction, effective kernel, low-frequency renormalization, two-pole spectrum, perturbative shift, channelwise P2, trace/anomaly, isotropy, selection rule.

**Issues Found:**

1. **(Cosmetic)** Notes Sec. 9 lists filenames `moving_throat_pde_stage3_bdg_*` but actual files are `moving_throat_pde_stage003_bdg_*`.

---

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The axisymmetric reduction, exact Schur-complement kernel, one-mode pole shift, grouped `P2` low-frequency coefficients, and isotropy diagnostics are mathematically consistent with stage 002 and with the audit output.
2. The logical flow is clean: reduced wall sector, exact elimination of stable matter modes, then the single-mode pole-shift example and grouped `P2` channel bookkeeping.
3. Assumptions are stated clearly enough for a reduced-sector note: stable normal modes, weak coupling for the perturbative shift, and an isotropic reference throat for the degeneracy statement.
4. The only issue I found is cosmetic: Section 9 still cites `moving_throat_pde_stage3_bdg_sympy_output.txt`, but the actual audit artifact is the staged `003` filename.

**Script Review:**

1. The script faithfully checks the claims it advertises: the Euler-Lagrange reduction, the effective matrix kernel, the low-frequency series, the exact two-pole formula, the perturbative wall-like shift, the grouped `P2` coefficients, and the isotropy/selection-rule checks.
2. I did not find a code bug or a tautological identity. The SymPy checks are independent enough to support the note's reduced algebra.

**Issues Found:**

- Cosmetic stale filename in the note's script-backed-status section.

**Questions:**

- None.

---

### Agent: GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

1. The checkpoint claim remains the exact reduced Schur-complement replay: axisymmetric kernel elimination, grouped real `P_2` channel kernels, the one-mode avoided-crossing formula, and the isotropic grouped-degeneracy criterion. I did not find a mathematical mismatch against the canonical stage card or the current audit outputs.
2. The earlier defect was bookkeeping, not derivational weakness. The note’s audit-status section had stale filenames and understated the actual audit surface.
3. The shared numerical layer is properly scoped. Its JSON sample values are probe-only test cases for root ordering, perturbative-validity breakdown, and Stage `004` outgoing-scaling checks; they are not used to derive the Stage `003` formulas.

**Script Review:**

1. The SymPy audit still gives a full symbolic replay of the checkpoint claim: Euler-Lagrange reduction, exact Schur complement, low-frequency series, exact two-pole spectrum, perturbative wall-like shift, grouped trace / anomaly bookkeeping, isotropy, and the harmonic selection rule.
2. The Mathematica mirror independently reproduces the same symbolic checks and passes cleanly.
3. The actual runnable defect was in the shared numerical-stress harness path resolution after the repo move into `research/pde_ledger/`. That pathing is now corrected in both the Python and Mathematica stress scripts, and both scripts run successfully against the repo-local sample JSON.

**Issues Found:**

- None. The stale filename issue is resolved, and the shared-stress harness now resolves the current repo-local sample path correctly.

**Questions:**

- None.
