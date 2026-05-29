---
unit_id: 151
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T19:25:17-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md]
  paper_appendix: present
---

# Audit unit 151 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_151.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_151}` at line 1336; no separate status row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.txt`

## What the paper claims

The stage card is a terse "finite mouth-profile corrections ledger step" that points to the two audit transcripts; the load-bearing math is in the notes. The notes state four boxed deliverables. (1) Expanding the exact self-consistent positive source `Σ_full(x) ∝ e^{-Φ_*(x)} = e^{-Π_* x - R_*(x)}` about the canonical exponential source `Σ_*(x) = Π_* e^{-Π_* x}/(1-e^{-Π_*})` gives the normalized first-order correction `Σ_act = Σ_*[1 - R̃_*] + O(R_*²)` with `R̃_* := R_* - <R_*>_*` and `<f>_* := ∫_0^1 Σ_* f dx`. (2) The first-order moment shifts are `δg_act = ∫_0^1 c δΣ_* dx = -Cov_*(c, R_*)` and `δS_act = ∫_0^1 K_q δΣ_* dx = -Cov_*(K_q, R_*)`, with kernels `c(x) = cos(πx/2)`, `K_q(x) = cosh((π/2)(1-x))/cosh(π/2)` and `Cov_*(f,h) = <fh>_* - <f>_*<h>_*`. (3) The bias retunes as `δΠ_act = -δg/g_*' = Cov_*(c,R_*)/g_*'` and the normalized mouth traction shift is `δT_{m,act} = A_T δg + B_T δS = -A_T Cov_*(c,R_*) - B_T Cov_*(K_q,R_*)`. (4) The whole mouth-side ambiguity reduces to the single object `R_*(x)`. The notes keep `Π_*`, `R_*`, `g_*'`, `A_T`, `B_T` symbolic; stage 151 asserts no numeric value of `Π_*`.

## What the script claims to verify

The SymPy docstring states it verifies the moment-shift identities `δg = -Cov_*(c,R)`, `δS = -Cov_*(K_q,R)` by (1) symbolically defining `δΣ = -Σ_*(R - <R>_*)` as a hand form and (2) numerically integrating `∫c δΣ` / `∫K_q δΣ` and checking they match the abstract covariance form to high precision. The Mathematica script builds `Σ_full = e^{-Φ}/Z` with `Φ = Π_* x + ε(r1 x + r2 x²)`, takes a genuine `Series` expansion in `ε` to first order, extracts `SigmaStar` and `deltaSigma`, and symbolically checks: the centering `<δΣ>_* = 0`; that the series-derived `deltaSigma` equals `-Σ_*(R - <R>)` (this is paper deliverable #1); the two moment-shift identities; and the bias/traction retunings.

## Paper ↔ script cross-check

| Paper deliverable | SymPy | Mathematica | Status |
|---|---|---|---|
| (1) `Σ_act = Σ_*[1-R̃_*]` derived from expansion of `e^{-Φ_*}` | hand-coded `δΣ = -Σ_*(R-Rbar)` at line 99-100; expansion NOT derived | `Series` of `e^{-Φ}/Z` → `deltaSigma`; checked `deltaSigma + Σ_*(R-<R>) == 0` (line 54) | `partial` (SymPy assumes it; Mathematica derives it) |
| (2) `δg = -Cov_*(c,R)`, `δS = -Cov_*(K_q,R)` | numeric `expect_close` at one anchor (lines 110-111) | symbolic `expectZero` (lines 58-59) | `match` (Mathematica exact; SymPy single-point numeric) |
| (3) `δΠ = Cov_*(c,R)/g_*'`, `δT_m = -A_T Cov(c,R) - B_T Cov(K,R)` | numeric `expect_close` (lines 116-122) | symbolic `expectZero` (lines 63-64) | `match` (Mathematica exact; SymPy single-point numeric) |
| (4) reduction to single object `R_*(x)` | print only (lines 127-129) | print only (lines 67-69) | `match` (narrative, no check needed) |

Centering `<δΣ>_* = 0` (mass preservation) is also checked by both engines (SymPy line 105, Mathematica line 53). `paper_alignment: aligned` — every load-bearing identity is exercised by at least the Mathematica engine, and the SymPy values use `Π_*` only as an arbitrary input anchor, not as a paper-asserted constant.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 81 | `expect_close(<1>_*, 1)` numeric | normalization (supports #2) | partial (one point) |
| A2 | sympy | 105 | `expect_close(<δΣ>_*, 0)` numeric | centering (supports #1/#2) | partial (one point) |
| A3 | sympy | 110 | `expect_close(δg_int, -CovcR)` numeric | #2 (δg) | partial (one point) |
| A4 | sympy | 111 | `expect_close(δS_int, -CovKR)` numeric | #2 (δS) | partial (one point) |
| A5 | sympy | 116 | `expect_close(δΠ, CovcR/gprime)` numeric | #3 (δΠ) | partial (one point) |
| A6 | sympy | 117-122 | `expect_close(δT, -AT·CovcR-BT·CovKR)` numeric | #3 (δT) | partial (one point) |
| A7 | math | 53 | `expectZero(∫δΣ dx)` symbolic | centering (#1/#2) | yes |
| A8 | math | 54 | `expectZero(deltaSigma + Σ_*(R-<R>))` symbolic | #1 (expansion form) | yes |
| A9 | math | 58 | `expectZero(δGInt + CovCR)` symbolic | #2 (δg) | yes |
| A10 | math | 59 | `expectZero(δSInt + CovKR)` symbolic | #2 (δS) | yes |
| A11 | math | 63 | `expectZero(δΠ - CovCR/gPrime)` symbolic | #3 (δΠ) | yes |
| A12 | math | 64 | `expectZero(δT + aT·CovCR + bT·CovKR)` symbolic | #3 (δT) | yes |

The Mathematica rows (A7-A12) are genuine exact-zero `FullSimplify === 0` checks over fully symbolic `piStar, r1, r2, gPrime, aT, bT`, so each non-tautologically exercises the paper claim across the whole parameter space. The SymPy rows (A1-A6) are tolerance-accepted numeric checks at a single fixed anchor (`Pi_star=1.50882951349316`, `r1=1.7`, `r2=-0.9`), with nonzero accepted diffs (~1e-42); they cannot distinguish the exact identity from a value that happens to coincide at that point, and A8's paper deliverable #1 (the exponential-expansion derivation) has no SymPy counterpart at all (SymPy hand-codes the result it should derive).

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:24` (mpmath import)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:50-55` (fixed numeric anchors)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:99-100` (hand-coded δΣ, not derived)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:108-122` (tolerance-accepted numeric `expect_close`)

**What's wrong:**
The paper's deliverables #2 and #3 are *exact symbolic identities* that must hold for all `Π_* > 0` and all residual coefficients (`δg = -Cov_*(c,R)`, `δS = -Cov_*(K_q,R)`, `δΠ = Cov_*(c,R)/g_*'`, `δT_m = -A_T Cov_*(c,R) - B_T Cov_*(K_q,R)`). The SymPy script verifies them only as a single-point floating-point spot-check: it imports `mpmath` (line 24), fixes `Pi_star = mp.mpf("1.50882951349316")`, `r1 = 1.7`, `r2 = -0.9`, `gprime`, `AT`, `BT` (lines 50-55), integrates with `mp.quad`, and accepts via `expect_close` with nonzero residuals. The saved transcript confirms the checks are numeric, not exact: line 11 shows `<delta_Sigma>_* ... diff 1.375847592069469627811980620665879893659e-42`, line 12 `delta_g_int ... diff 2.152394441202919020938848639933311177646e-42`, line 14 `deltaPi ... diff 2.29588740394978028900143854926219858949e-41`, line 15 `deltaT ... diff 5.739718509874450722503596373155496473724e-42`. A single numeric anchor cannot distinguish the true identity from one that merely coincides at that point.

Compounding this, the SymPy script does not even *cover* paper deliverable #1: the linearized correction is hand-asserted at lines 99-100 (`delta_Sigma(x) = -Sigma_star(x) * (R_residual(x) - Rbar)`) rather than derived by expanding `e^{-Φ_*} = e^{-Π_* x - R_*(x)}` about the canonical source. The SymPy script assumes the very form whose origin the paper's Section 1 derives. (The Mathematica script derives it correctly via `Series`, lines 35-37, and checks it at line 54.)

**Why this matters:**
The two-engine policy requires each engine to independently establish the claim. As it stands only the Mathematica engine performs an exact symbolic verification; the SymPy engine provides one numeric data point and assumes deliverable #1. If the paper identity were subtly wrong off the chosen anchor (e.g., a `Π_*`-dependent factor that vanishes only near `Π_* ≈ 1.509`), the SymPy check would still PASS while the identity is false in general. SymPy is fully capable of the closed-form integrals here (`∫_0^1 cos(πx/2) e^{-Π x} dx`, `∫_0^1 cosh((π/2)(1-x)) e^{-Π x} dx`, and `series(exp(-Φ)/Z, ε, 0, 2)` all have closed forms for symbolic `Π`), so the numeric shortcut is not forced by any computational limitation.

**Required change:**
Rewrite the SymPy script to perform a genuine symbolic derivation matching the Mathematica engine: declare `Pi_star` (positive), `r1, r2, gprime, AT, BT` (real, `gprime != 0`) as `sp.Symbol`s; build `Phi = Pi_star*x + eps*(r1*x + r2*x**2)`, `Z = sp.integrate(sp.exp(-Phi), (x,0,1))`, `Sigma_full = sp.exp(-Phi)/Z`, take `sp.series(Sigma_full, eps, 0, 2).removeO()` and extract the `eps^0` term as `Sigma_star` and the `eps^1` coefficient as `delta_Sigma`; then assert via `sp.simplify(...) == 0` (or `.equals(0)`): (i) `delta_Sigma + Sigma_star*(R - Rbar) == 0` (deliverable #1, currently missing on the SymPy side), (ii) `<delta_Sigma>_* == 0` (centering), (iii) `delta_g_int + CovcR == 0` and `delta_S_int + CovKR == 0` (deliverable #2), (iv) `deltaPi - CovcR/gprime == 0` and `deltaT + AT*CovcR + BT*CovKR == 0` (deliverable #3), where every integral is `sp.integrate(..., (x,0,1))` over symbolic `Pi_star, r1, r2`. Remove the `mpmath`/`expect_close`/fixed-anchor scaffolding (lines 24, 34-38, 50-55, and the numeric `mp.quad` usages). See the directive for the full claim manifest.

**Verification:**
The verifier reruns `redteam exec-sympy 151`. The refreshed transcript must show exact `= 0` (or `True`) results — not tolerance-accepted nonzero diffs — for at least the six identities above, and a new check corresponding to deliverable #1 (`delta_Sigma + Sigma_star*(R - <R>) == 0`) must appear. The script must exit 0 with no `mpmath` import remaining.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. The two engines take materially different derivation paths for deliverable #1. The SymPy script hand-codes the correction: `def delta_Sigma(x): return -Sigma_star(x) * (R_residual(x) - Rbar)` (lines 99-100). The Mathematica script instead derives it from first principles: `Phi[x_] := piStar*x + epsilon*(r1*x + r2*x^2); unnorm[x_] := Exp[-Phi[x]]; Z = Integrate[unnorm[x], {x,0,1}]; SigmaFull[x_] := unnorm[x]/Z; SigmaSeries = Normal[Series[SigmaFull[x], {epsilon,0,1}]]; deltaSigma = Coefficient[SigmaSeries, epsilon, 1]` (lines 31-37), then independently checks `deltaSigma + SigmaStar*(RResidual[x] - rBar) === 0` (line 54). The shared content (kernel definitions, covariance algebra) is dictated by the physics, not copied choreography. The Mathematica engine is the stronger of the two and is a valid independent derivation.

## Engine cross-check

Both engines pass and agree on the verified structure. Mathematica returns exact symbolic `= 0` for centering, the deliverable-#1 form, both moment shifts, and both retunings (output lines 5-16, all `PASS`). SymPy returns the same identities numerically at its anchor: `Cov(c,R) = -0.0596749...`, `Cov(K,R) = -0.0462935...`, `delta_g_int = -Cov(c,R)`, `delta_S_int = -Cov(K,R)`, `deltaPi = -0.835252...`, `deltaT = -0.248725...` (output lines 9-18). No sign, factor, or structural disagreement; the only difference is rigor (exact vs. single-point numeric), captured in F1. No `engine_disagreement` finding.

## Verdict justification

Verdict is `findings` (one finding, medium). I attacked the identities directly: re-derived `δg = ∫c δΣ = -(<cR>_* - <R>_*<c>_*) = -Cov_*(c,R)` and confirmed it requires `<Σ_*>_* = 1` (checked at SymPy line 81 / built into Mathematica's normalization) — the algebra is correct and the Mathematica `expectZero` checks exercise it symbolically across the full parameter space, so the paper's claims hold. The kernels (`c`, `K_q`), the canonical source `Σ_*`, and the covariance definition all match the notes exactly; no `paper_misalignment`, no `tautological_check`, no `hardcoded_result` (the numeric `Pi_star` is an arbitrary input anchor, not a checked-against-itself answer), no `symbol_assumption_error`, no `engine_disagreement`, no `mathematica_transliteration` (the engines derive deliverable #1 differently), and outputs are fresh (SymPy 11:27/11:30, Mathematica 10:03/11:31). The single defect is that the SymPy engine verifies exact symbolic identities only as a tolerance-accepted single-point numeric spot-check (nonzero diffs ~1e-42) and additionally never derives paper deliverable #1 (it hand-codes the correction form), so it does not independently establish the claim — `insufficient_verification`.

## Self-test notes

I walked the proposed SymPy symbolic rewrite through the traps. (1) Variable independence: the fix uses `sp.series`/`sp.integrate`, no `sp.diff(EXPR, VAR)` with a phantom variable, so no identically-zero-derivative trap. (2) Symmetry/parity: all integrals are over the bounded domain `[0,1]`, not an unbounded symmetric domain, so parity-cancellation traps do not apply; closed forms exist for `∫_0^1 cos(πx/2)e^{-Π x}dx`, `∫_0^1 cosh((π/2)(1-x))e^{-Π x}dx`, and the `ε`-series of `e^{-Φ}/Z`. (3) Trivial-case pre-check: with `Π_*` symbolic, `δg + Cov(c,R)` reduces to `-(<cR> - <c><R>) + (<cR> - <c><R>) = 0` exactly given `<Σ_*> = 1`, and deliverable #1 `δΣ + Σ_*(R-<R>)` is exactly the `ε^1` series coefficient by construction, reducing to literal 0. (5) Paper round-trip: the rewrite introduces no new constant — it keeps `Π_*, r1, r2, g_*', A_T, B_T` symbolic, matching the notes which assert no numeric value, so no new `paper_misalignment` is created.
