# Review: Stage 162 — Transfer shape theorem

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage162_transfer_shape_theorem.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage162_transfer_shape_sympy_audit.py`

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

**Notes Derivation Review:**

1. **Equation-level correctness.** All algebraic steps verified by hand against the primitive definitions and confirmed by SymPy.

   - *Wall-normalized variable definitions (Section 1).* GW_hat = GW/(OW^2 sqrt(K)), GU_hat = GU/(OU OW sqrt(K)), R_hat = R/(OU OW). Substituting back into P_r = OU^2 GW + R GU reproduces P_r = sqrt(K) OU^2 OW^2 (GW_hat + R_hat GU_hat). Similarly Delta_r = OU^2 OW^2 (1 - R_hat^2). The ratio N0/K = T^2 follows immediately. All signs and factors correct.

   - *Weak-axisymmetric slopes (Section 2).* Taking log-derivatives of each wall-normalized variable: w = gW - oW - kappa1/2, u = gU - oU/2 - oW/2 - kappa1/2, c = rr - oU/2 - oW/2. Each correctly accounts for all denominator factors in the respective definitions. Verified by direct differentiation.

   - *Transfer-shape slope tau (Section 3).* The logarithmic derivative of T = numerator/denominator splits into an alpha-weighted wall-leg slope, a beta-weighted mixed-leg slope, and the 2R_hat^2/(1-R_hat^2) coupling term from the denominator. The identity nu_r = kappa1 + 2 tau_r follows from N0 = K T^2 by the chain rule. This is the central identity, and it is correct.

   - *Defect collapse (Section 4).* Substituting nu_r = kappa1 + 2 tau_r into Xi_1 = sum rho_r (nu_r - kappa1) with sum rho_r = 1 gives Xi_1 = 2 sum rho_r tau_r. Pure algebra, correct.

   - *Equivalence to Stage 159/160/161 (Section 5).* Verified: M_r = GW_hat, I_r = R_hat GU_hat / GW_hat, H_r = R_hat^2. Slope identities m = w, i = (u+c) - w, h = 2c all checked against Stage 161 definitions. The slippage form tau = m + I/(1+I) i + H/(1-H) h recovers exactly from the alpha/beta form. sigma_r = 2 tau_r is consistent with Stage 160's sigma_r = nu_r - kappa1.

   - *Corollaries (Section 7).* Dominant-port limit, square-root mixed-leg recovery (u=c=0 implies tau=w), and common-leg co-scaling (c=0, u=w implies tau=w) are all straightforward specializations of the general formula. No errors.

2. **Logical flow.** Clean linear progression: recall Stage 161 data, introduce wall normalization, derive factorization, compute slopes, collapse defect, prove equivalence to prior stages, state corollaries. Each section depends only on what precedes it plus the Stage 161 results.

3. **Assumptions.** Same conservative-shape-preserving branch as Stages 158-161. The factorization N0 = K T^2 is exact (no approximation); the weak-axisymmetric expansion is first-order in epsilon as throughout the chain. No new assumptions introduced.

4. **Notation.** Consistent with prior stages. The hat notation (GW_hat, etc.) for wall-normalized variables is new to this stage and clearly defined. The fraktur slopes (w, u, c) are new compact names for the wall-normalized counterparts of the primitive slopes; their relations to the prior fraktur slopes (m, i, h, gW, etc.) are explicitly given.

5. **Physical interpretation.** Sound. The wall-normalized transfer shape T_r is a dimensionless object measuring how the outgoing port deviates from what the wall baseline alone would predict. The remaining theorem gate is cleanly recast: "are the transfer shapes rigid?" This is a genuine sharpening of the Stage 161 gate.

**Script Review:**

1. **Faithful implementation.** The script builds the primitive port data (P, Delta, N0) from generic symbolic variables, constructs the wall-normalized variables (GWh, GUh, Rh, T) from the same definitions as the notes, and verifies N0/K = T^2 symbolically. The weak-axisymmetric perturbation is implemented via multiplicative (1 + eps lam slope) factors applied to each primitive, matching the notes' convention. The slopes w, u, c are constructed from the explicit formulas in Section 2. The tau formula uses the alpha/beta weights built from the hatted variables. All definitions match the notes.

2. **No hardcoded values or tautologies.** The central check (line 91) compares nu_direct -- obtained independently by Taylor-expanding log(N0A) in epsilon to first order -- against kappa1 + 2 tau, where tau is assembled from the wall-normalized slope formulas. These two expressions are built by independent code paths (series expansion of a ratio of polynomials vs. algebraic combination of named slopes), so the check is non-tautological. The slippage equivalence check (lines 101-103) similarly builds tau_slippage from a separate formula involving I, H, m, i, h and compares it to the alpha/beta form of tau. The weighted defect check (lines 107-112) uses three generic ports with symbolic rho values summing to 1.

3. **Symbol assumptions.** K, OU2, OW2, GW, GU are declared positive+real. R is real (allowing negative coupling, appropriate). eps, lam, kappa1, and all slopes are real. No inappropriate constraints. The positivity of K, OU2, OW2 is physically motivated and avoids sign ambiguities in sqrt.

4. **Coverage.** The script tests four distinct claims from the notes: (a) exact factorization N0/K = T^2; (b) slope identity nu = kappa1 + 2 tau via independent series expansion; (c) equivalence to Stage 159/160/161 slippage form; (d) weighted defect identity Xi_1 = 2 sum rho tau. This covers Sections 1, 3, 4, and 5 of the notes. Section 2 (slope definitions) is implicitly verified through the slope identity check. The corollaries in Section 7 are simple specializations and do not need separate script coverage.

5. **Output agreement.** All four expect_zero checks print "= 0" and the script exits with code 0. No warnings or errors.

6. **Minor observation.** The script does not independently verify the per-port rigidity sufficient condition (Section 6, tau_r = 0 for all r implies Xi_1 = 0), but this follows trivially from the weighted defect identity already tested and does not need a separate check.

**Issues Found:** None.

---
-->

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:** The wall-normalized transfer-shape theorem is correct. The factorization `N0/K = T^2`, the weak-axisymmetric slope identity `nu_r = kappa_1 + 2 tau_r`, the slippage-form equivalence, and the weighted defect collapse `Xi_1 = 2 sum_r rho_r tau_r` all line up with the note.

**Script Review:** The audit script independently rebuilds the wall-normalized variables, checks the `N0/K = T^2` factorization, differentiates the perturbed `N0` directly, verifies the slippage-form identity, and confirms the weighted defect law. The saved output matches the note.

**Issues Found:** None.
