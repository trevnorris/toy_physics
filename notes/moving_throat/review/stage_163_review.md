# Review: Stage 163 — Effective Transfer-Shape Collapse

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage163_effective_transfer_shape_collapse.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage163_effective_transfer_shape_collapse_sympy_audit.py`
- **Output:** `scripts/moving_throat/output/moving_throat_pde_stage163_effective_transfer_shape_collapse_sympy_audit.txt`
- **Prior:** `notes/moving_throat/moving_throat_pde_stage162_transfer_shape_theorem.md`

## Review Checklist

- [x] Equation-level correctness (signs, factors, indices, limits)
- [x] Logical flow from prior stage(s)
- [x] Assumptions stated and justified
- [x] Notation consistent with prior stages
- [x] Physical interpretation sensible
- [x] SymPy script faithfully implements notes
- [x] Script runs without error
- [x] Script output matches notes claims
- [x] No missing edge cases or branches

## Agent Reviews

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS

---

**A. Notes Derivation Review**

**A.1 Equation-level correctness**

All displayed equations verified:

- **Section 1 (multi-port collapse):** Starting from the Stage-162 factorization N_0^{(r)} = K * T_r^2, the sum N_0 = K * sum_r T_r^2 is immediate. The definition T_eff^2 := sum_r T_r^2 = N_0/K is consistent. The logarithmic perturbation delta ln T_eff^2 / (eps*lambda) = 2 * sum_r rho_r^{(N)} * tau_r follows correctly: perturbing T_r -> T_r * exp(eps*lam*tau_r) gives T_eff^2 -> sum_r T_r^2 * exp(2*eps*lam*tau_r), and differentiating ln of this at eps=0 yields 2*lam * sum_r (T_r^2/sum_s T_s^2) * tau_r, recovering the Stage-162 weights rho_r^{(N)} exactly.

- **Section 2 (one-port continuum formula):** On the one-port branch, T^2 = beta_0/K_0. Substituting K_0 = K_eta^{eff}/mu_eta and beta_0 = (mu_W/mu_eta)(K_eta^{eff}/K_W^{eff}) * Z_W(1+rho)^2/(1-eps_W)^2, the K_eta^{eff} and mu_eta factors cancel exactly, giving T^2 = (mu_W/K_W^{eff}) * Z_W(1+rho)^2/(1-eps_W)^2. Using Omega_W^2 = K_W^{eff}/mu_W yields T^2 = Z_W(1+rho)^2/(Omega_W^2(1-eps_W)^2). All signs and factors confirmed.

- **Section 3 (direct slope law):** From ln T^2 = ln Z_W + 2*ln(1+rho) - ln Omega_W^2 - 2*ln(1-eps_W), linearizing with delta ln Z_W = eps*lam*zeta_W, delta ln Omega_W^2 = eps*lam*omega_W, delta rho = eps*lam*rho_1, and delta eps_W = eps*lam*eps_W_1 yields Xi_1 = zeta_W - omega_W + 2*rho_1/(1+rho) + 2*eps_W_1/(1-eps_W). Signs verified individually: zeta_W enters with +1 (numerator), omega_W with -1 (denominator), rho_1 gets +2/(1+rho) from d/d(rho)[2*ln(1+rho)], and eps_W_1 gets +2/(1-eps_W) from d/d(eps_W)[-2*ln(1-eps_W)].

- **Section 4 (selected-branch reformulation):** Using R_target = Lambda*(1-eps_eta)*(1-eps_W)^2/(Z_W*(1+rho)^2) with Lambda = (27*pi^2*G*c_s^5/(20*a^5*c^5))*Omega_W^2, solving for T^2 gives (27*pi^2*G*c_s^5/(20*a^5*c^5))*(1-eps_eta)/R_target. The slope law Xi_1 = -eta_1/(1-eps_eta) - R_1 follows from perturbing ln(const*(1-eps_eta)/R_target) with the front factor inert at linear grouped P_2 order on the isotropic branch.

- **Section 5 (zero-defect theorem):** Setting Xi_1 = 0 in the direct-port form gives delta ln[Z_W(1+rho)^2/(Omega_W^2(1-eps_W)^2)] = 0. In the selected-branch form, delta ln R_target = -delta(eps_eta)/(1-eps_eta). Corollaries 5.1 (target-rigid implies eps_eta-rigid) and 5.2 (eps_eta-rigid implies target-rigid) are both immediate. The bidirectionality is genuine since the selected-branch Xi_1 is linear in both eta_1 and R_1.

- **Section 6 (quadrupole defect):** Substituting Xi_1 = 2*tau into the Stage-156 pattern: Delta_Q^{(20)} = eps*2*tau = 2*eps*tau, Delta_Q^{(21)} = (eps/2)*2*tau = eps*tau, Delta_Q^{(22)} = -eps*2*tau = -2*eps*tau. All three correct.

**A.2 Logical flow from Stage 162**

Stage 162 established: (1) N_0^{(r)} = K * T_r^2 with an explicit transfer shape, (2) Xi_1 = 2*sum_r rho_r^{(N)} * tau_r where tau_r is the transfer-shape slope. Stage 163 picks up exactly these results and performs two clean operations: collapses the multi-port sum to a single effective shape T_eff^2, then evaluates that shape on the one-port continuum branch using the Stage-21 formulas for K_0 and beta_0. The logical chain is tight with no gaps.

**A.3 Assumptions**

The key assumptions are: (1) the Stage-162 portwise factorization, (2) the Stage-21 continuum-kernel formulas for K_0 and beta_0, (3) the minimal one-port isotropic branch specialization, (4) scalar observables (a, c_s) are inert at linear grouped P_2 order. Assumptions (1)-(3) are clearly sourced from prior stages. Assumption (4) is stated and physically justified on the isotropic branch.

**A.4 Notation**

Consistent with Stage 162. New notation introduced: T_eff (effective transfer shape), zeta_W, omega_W, rho_1, eps_W_1 for the four microscopic drift channels, eta_1 and R_1 for the selected-branch perturbation parameters, tau_eff for the effective slope. All clearly defined at point of introduction.

**A.5 Physical interpretation**

The observation that eps_eta drops out of the direct continuum transfer-shape formula (Section 2) is a genuine structural result -- the wall-U dressing lane affects port transfer only indirectly. The identification of four microscopic drift channels is clean and physically transparent. The dual interpretation via direct-port and selected-branch variables is well-motivated and provides a useful cross-check.

---

**B. Script Review**

**B.1 Faithful implementation**

The script tests four core claims from the notes:

1. **Multi-port collapse (lines 38-48):** Uses two ports with independent T_1, T_2, tau_1, tau_2. Constructs T_eff^2 = T_1^2*exp(2*eps*lam*tau_1) + T_2^2*exp(2*eps*lam*tau_2), differentiates ln(T_eff^2) at eps=0, divides by lam, and checks against 2*(rho_1*tau_1 + rho_2*tau_2). Two ports suffice since the algebra extends by linearity.

2. **One-port continuum formula (lines 54-72):** Constructs K_0 and beta_0 from Stage-21 definitions, computes T^2 = beta_0/K_0, and verifies it equals mu_W/K_W * Z_W*(1+rho)^2/(1-eps_W)^2. Also checks the Omega_W^2 form.

3. **Selected-branch reformulation (lines 78-89):** Defines Lambda and R_target from Stage-21 formulas, constructs T^2_selected, and verifies both expressions agree when R_target is expanded.

4. **Slope laws (lines 94-119):** Constructs perturbed T^2 with exponential perturbations in the four drift variables, differentiates at e=0, and checks against the four-term slope formula. Also checks the selected-branch slope law.

**B.2 Code correctness**

- The `expect_zero` function applies `sp.simplify(sp.expand(...))` before comparing to 0, appropriate for these algebraic identities.
- **Lines 70-72 precedence:** The `.subs({OmegaW2: KW/muW})` method call binds to the parenthesized sub-expression `(OmegaW2 * (1 - epsW)**2)` due to Python method-call precedence being higher than arithmetic operators. Since `OmegaW2` appears inside that parenthesized group, the substitution replaces it with `KW/muW`, making the denominator `KW*(1-epsW)^2/muW`. The test is correct despite the visually misleading line break.
- **Lines 86-89 selected-branch check:** `T2_direct` does not contain `Rtarget`, so `.subs({Rtarget: Rtarget_def})` on it is a harmless no-op. The substitution is effective only on `T2_selected`, expanding it to the direct-port form. This is algebraically valid and genuinely tests the claimed equivalence.
- The perturbation technique (exponential perturbation + differentiation at e=0) is standard and correctly implemented.

**B.3 Hardcoded values**

No hardcoded numerical values. The factor 27*pi^2/20 appears in both script and notes, sourced from Stage-21 formulas.

**B.4 Tautological checks**

None found. Each check computes from one representation and compares to a distinct one:
- Check 1: perturbation-based Xi vs. weight-based formula.
- Check 2: ratio beta_0/K_0 vs. simplified algebraic form.
- Check 3: direct-port T^2 vs. selected-branch T^2 (after expansion).
- Checks 4-5: perturbation-based slopes vs. closed-form slope formulas.

**B.5 Symbol assumptions**

Appropriate throughout. `rho` and `epsW` are `real=True` without positivity (correct -- these are dimensionless ratios). Positive-definite quantities (mu_W, mu_eta, K_eta_eff, K_W_eff, G, c_s, a, c) are correctly `positive=True`. The perturbation parameter `e` is correctly `real=True`.

**B.6 Output agreement**

All six checks return 0. Exit code 0. Output matches all notes claims that are tested.

**B.7 Coverage**

Sections 1-4 of the notes are directly tested by the script. Section 5 (zero-defect conditions and corollaries) and Section 6 (quadrupole defect substitution) are not independently scripted, but these are trivial corollaries -- setting Xi_1=0 and substituting Xi_1=2*tau into already-established Stage-156 results, respectively. Coverage is adequate.

---

## Issues Found

None. The derivation is clean, the script faithfully tests the core algebraic claims, and all checks pass. The logical flow from Stage 162 is tight, the notation is consistent, and no tautological, hardcoded-value, or symbol-assumption issues are present.

## Questions

None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The effective transfer-shape collapse is correct. The multi-port reduction to `T_eff^2 = sum_r T_r^2`, the one-port continuum form `T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-eps_W)^2]`, the selected-branch reformulation, and the direct/selected slope laws all match the derivation chain.

**Script Review:** The audit script independently checks the multi-port collapse, the one-port continuum factorization, the selected-branch identity, and the weak-axisymmetric slope laws. The saved output is clean and agrees with the note.

**Issues Found:** None.
