---
unit_id: 204
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T09:55:17-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 204

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose documents — this directive only adds a script.

After writing the new script, RUN it (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting it to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.wl`

**Issue:** Stage 204 is non-status-only and non-checkpoint and computes a large body of independently verifiable symbolic math, but only a SymPy script exists. The card states "Mathematica audit: none yet." (`paper/stages/stage_204.tex:11`). The dual-engine contract requires a Mathematica `.wl` wherever Mathematica can independently verify the stage; it can verify every claim below natively. Add the `.wl`.

**Anti-transliteration guard:** The new script MUST derive each result independently using native Mathematica primitives (`Solve`/`Reduce`, `Series` + `SeriesCoefficient`/`Coefficient`, `D[]`, `Integrate`, `FindRoot`, matrix operations) via a DIFFERENT decomposition than the `.py`. In particular: do NOT mirror the SymPy script's "build `const_0` and `const_0*exp(sigma*tau)` then subtract" choreography. Instead, recover each dependent exponent `sigma` as the τ-logarithmic-derivative of the closed-form graph quantity — e.g. `sigmaDelta = D[Log[deltaGraph[tau]], tau]` simplified to a τ-independent constant, then compare that recovered constant against the paper's formula — and recover the predictor result by extracting series coefficients with `SeriesCoefficient` rather than subtracting a guessed tail. A line-by-line port of the SymPy algebra (same variable choreography, same direct-minus-exp subtraction) is rejected as `mathematica_transliteration`.

**Required change:** Create the `.wl` at the Target path. It must declare the symbols with the same domains the SymPy script uses (positivity on `chi0_star, deltaU_star, E_star, F_star`, the base-point coordinates `lambda0, c0, gamma0, KU0, KW0`, the target constants, `L, sigma`; reals on `tau` and the five log-slopes) and verify the claim manifest below with `expectZero`-style checks (`FullSimplify[lhs - rhs]` must be `0`, with `Exit[1]` on any nonzero). Use `expectZero` patched to strip `ConditionalExpression[0, ...]` per project idiom.

**Claim manifest** (each must be independently verified; symbols match the paper card / notes / appendix):

- **M1 — constant free log-slopes.** With `lambdaW[tau] = lambda0 Exp[s_lambda tau]` (and analogously for `c_etaU, gamma, K_U, K_W`), verify `D[Log[lambdaW[tau]], tau] - s_lambda == 0` and the four analogues. (paper notes lines 124-135; appendix line 533-543.)

- **M2 — dependent exponent σ_δ.** Define `aStar = (1 + deltaU_star)/(1 + chi0_star)` and the closed-form graph quantity `deltaGraph[tau] = (Ctr_target / ((gamma[tau] c_etaU[tau] / K_U[tau])^(1 + deltaU_star)))^(1/(1 + chi0_star))`. Verify that `D[Log[deltaGraph[tau]], tau]` is τ-independent and equals `sigmaDelta = -aStar (s_gamma + s_c - s_U)`. (notes lines 172-178; appendix line 609-612.)

- **M3 — dependent exponents σ_T, σ_{K_η}, σ_μ.** With `TGraph[tau] = L^2 K_U[tau] deltaGraph[tau]/Pi^2`, `KetaGraph[tau] = c_etaU[tau]^2/(K_U[tau] epsEta_target)`, and `muGraph[tau] = Cnt_target c_etaU[tau]^2 K_W[tau]^2/(epsEta_target K_U[tau] lambdaW[tau]^2) * ((gamma[tau]^2 lambdaW[tau]^2 sigma)/(K_U[tau] K_W[tau]))^(-E_star) * deltaGraph[tau]^F_star`, verify each `D[Log[·], tau]` is τ-independent and equals respectively `sigmaT = s_U + sigmaDelta`, `sigmaKeta = 2 s_c - s_U`, and `sigmaMu = 2 s_c - s_U + 2 s_W - 2 s_lambda - E_star (2 s_gamma + 2 s_lambda - s_U - s_W) + F_star sigmaDelta`. (notes lines 180-194; appendix lines 613-625.)

- **M4 — finite monomial invariance.** Substitute the full graph-lifted ray `x_s^graph(tau)` (eight components, the four dependent ones carrying `e^{sigma tau}` factors) into the three target monomials and verify each is τ-independent and equal to its target constant: `Ctr(tau) == Ctr_target`, `Cnt(tau) == Cnt_target`, `epsEta(tau) == epsEta_target`, for all `tau`. The monomial definitions match the SymPy Sec. III forms (lines 141-147 of the `.py`); reconstruct them from the paper §4 monomials independently rather than copying the SymPy expressions verbatim. (notes lines 257-296.)

- **M5 — quotient-map kernel.** With the carried Stage 192/243 quotient matrix `Mstar` (3×8, rows as in the SymPy Sec. IV) and the constant ray tangent `dxRay = {s_lambda, s_c, s_gamma, s_U, sigmaKeta, s_W, sigmaMu, sigmaT}`, verify `Mstar . dxRay == {0,0,0}`. (notes lines 299-304; appendix line 571-573.)

- **M6 — primitive direction table.** For each of the five coordinate directions `e_lambda, e_c, e_gamma, e_U, e_W` (unit log-slope vectors), verify the four exponents `(sigmaDelta, sigmaT, sigmaKeta, sigmaMu)` equal the table literals:
  - `e_lambda`: `(0, 0, 0, -2 - 2 E_star)`
  - `e_c`: `(-aStar, -aStar, 2, 2 - F_star aStar)`
  - `e_gamma`: `(-aStar, -aStar, 0, -2 E_star - F_star aStar)`
  - `e_U`: `(aStar, 1 + aStar, -1, -1 + E_star + F_star aStar)`
  - `e_W`: `(0, 0, 0, 2 + E_star)`
  (notes lines 332-338; appendix table.)

- **M7 — first-order predictor agreement.** With `Phi0 = 1 + eps`, `Phi1 = L0 (1 + eps)`, `tauAff = (1 - Phi0)/Phi1`, `tauLog = -Log[Phi0]/L0`, verify via `Series`/`SeriesCoefficient` that the difference `tauLog - tauAff` has vanishing constant and linear coefficients in `eps` and leading (eps^2) coefficient `-1/(2 L0)`, i.e. `SeriesCoefficient[tauLog - tauAff, {eps, 0, 2}] == -1/(2 L0)` and `SeriesCoefficient[..., {eps, 0, 0}] == 0`, `SeriesCoefficient[..., {eps, 0, 1}] == 0`. Optionally also confirm the eps^3 coefficient `== 2/(3 L0)`. (notes lines 472-487; appendix lines 668-674.) Also verify the first-order log-ray completeness M0-style check `Series[Y0 Exp[(Y1/Y0) epsStep], {epsStep, 0, 1}]` reproduces `Y0 + Y1 epsStep`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 204` and confirm the new `.wl` appears at the Target path, contains substantive checks for M1-M7, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.wl`
- summary: Created the independent Mathematica audit for Stage 204 covering M1-M7 with log-derivative exponent recovery, graph-ray monomial invariance, quotient-kernel, primitive-table, and predictor coefficient checks.
- deviation: none
