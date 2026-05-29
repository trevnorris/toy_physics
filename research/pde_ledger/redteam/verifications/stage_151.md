---
unit_id: 151
batch: IV.6
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T23:45:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 151

> Note: this file supersedes a prior stale verification that described the
> earlier mpmath single-point implementation (and mislabeled the unit as having
> 3 findings). The original auditor report for unit 151 records exactly ONE
> finding (F1, insufficient_verification). The resolution evolved across
> directive iterations: delta1's fully-symbolic SymPy target was intractable
> (SymPy hangs on symbolic Pi_star), so the agreed final target (delta2 + delta3)
> is: Mathematica = full all-Pi_star symbolic verifier; SymPy = an EXACT
> multi-point cross-check at 5 rational Pi_star samples. Judged against that
> final target below.

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py` was fully rewritten (diff reworks ~271 lines). Removed: `import mpmath`, the `expect_close` tolerance helper, the fixed numeric anchors (`Pi_star = mp.mpf("1.50882951349316")`, `r1=1.7`, `r2=-0.9`, `gprime`, `AT`, `BT`), and all `mp.quad` numeric integration. Added:
- Exact-zero helper `expect_zero` = `sp.simplify(sp.cancel(sp.together(expr)))` then `assert res == 0` — no float, no tolerance, no `.is_zero` (`.py:37-40`).
- All of `r1, r2` (real), `gprime` (real, nonzero), `AT, BT` (real), `eps`, `x` kept SYMBOLIC; `k = sp.pi/2` (`.py:45-50`).
- `PI_SAMPLES = [1/2, 1, 3/2, 2, 5/3]` as exact rationals (`.py:121-127`).
- Per-sample fresh build: `Sigma_star`/`delta_Sigma` derived from an eps-first `sp.series` of `e^{-Pi*x}·e^{-eps·pert}/(Z0+eps·Z1)` (`.py:133-141`); M1–M7 asserted via `expect_zero` with `r1,r2,A_T,B_T,gprime` symbolic (`.py:153-186`).
- delta3 anti-footgun comment block present verbatim right after the docstring (`.py:13-23`).
- Mathematica `.wl` NOT touched (confirmed: absent from `git status` modified list; mtime 10:03, before the 23:18 `.py` edit).

**Assessment:**
Confirmed against the four `verified` criteria:

1. SymPy is now an exact multi-point check: 5 rational `Pi_star` samples, `r1/r2/A_T/B_T/gprime` symbolic, exact `expect_zero`. No mpmath, no `expect_close`, no numeric tolerance, no single fixed anchor.

2. M1–M7 are genuine, non-tautological identities. M1 (`delta_Sigma + Sigma_star*(pert - Rbar)`) verifies the `sp.series`-derived first-order coefficient equals the centered hand form — this is the deliverable-#1 derivation that previously had NO SymPy counterpart. M4/M5 integrate the series-derived `delta_Sigma` against the `cos`/`cosh` kernels and compare to independently-formed covariances `CovcR`/`CovKR`; the two sides have different integrands and agree only if the paper identity holds. M6/M7 are the bias/traction retunings. None is a definition checked against itself. The exec log shows exactly 35 `= 0` lines (5 samples × 7 M-items) and exit 0.

3. Mathematica script unchanged and still does the full all-`piStar` symbolic proof: `$Assumptions = piStar > 0`, symbolic `Phi`, `Series` expansion, `Integrate` over symbolic `piStar`, `expectZero` via `FullSimplify[...] === 0`; all 6 checks PASS, exit 0.

4. delta3 anti-footgun comment block present near the top of the `.py`.

Collateral edit (non-blocking, see Side Observations): Codex introduced a custom monkeypatched `sp.integrate` (`_exact_unit_integrate`, `.py:52-119`) for the `(x,0,1)` case rather than relying on stock `sp.integrate`. I verified its math is correct: `_poly_exp_moment(rate,deg)` computes ∫₀¹ xⁿe^{ax}dx by the correct integration-by-parts recursion (`e^a/a − (n/a)·I_{n-1}`; base `(e^a−1)/a`; `a=0 → 1/(n+1)`), and `_expand_linear_exponentials` rewrites the trig/hyperbolic kernels to linear-in-x exponentials handled exactly by that formula. It is a mathematically sound exact evaluator (a speed workaround for the 600s cap — delta3 reported 126s runtime), not a tautology shim; it cannot manufacture a false zero because the compared sides have genuinely different integrands.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `[Pi=1/2] M1 delta_Sigma + Sigma_star*(R - <R>) = 0`
- `[Pi=1] M4 dg + CovcR = 0`
- `[Pi=5/3] M6 deltaPi - CovcR/gprime = 0`
- `[Pi=5/3] M7 deltaT + A_T*CovcR + B_T*CovKR = 0`
All 35 M-lines print `= 0` — exact, with NO `diff ...e-NN` tolerance residuals anywhere. The old ~1e-42 numeric diffs are gone.

**Mathematica:** exit=0. Notable lines:
- `deltaSigma + SigmaStar*(R - <R>) = 0` / `PASS`
- `deltaGInt + Cov(c,R) = 0` / `PASS`
- `deltaPi - Cov(c,R)/gPrime = 0` / `PASS`
- `deltaT + aT*Cov(c,R) + bT*Cov(K,R) = 0` / `PASS`
All 6 `expectZero` checks pass over fully symbolic `piStar>0, r1,r2,gPrime,aT,bT`; `Stage 151 Mathematica audit passed.`

**Output freshness:** confirmed. `.py` mtime 23:18:31 < SymPy `.txt` 23:29:23. `.wl` mtime 10:03:12 < Mathematica `.txt` 23:29:23. Both transcripts regenerated post-fix.

## Material-change assessment

`material_change`: false. No derived constant or downstream-consumed numeric result changed — the fix only strengthened the SymPy engine's rigor (single mpmath point → exact 5-sample multi-point) and updated the docstring/comment. The verified identities and the Mathematica symbolic proof are identical to before. (The prior stale verification's `material_change: true` was an artifact of the now-removed hardcoded `deltaPi`/`deltaT` numeric values; the current script asserts no such numbers.) No downstream unit is affected.

## Side observations (non-blocking)

- The custom monkeypatched integrator (`.py:52-119`) is a substantial collateral implementation beyond delta3's literal text (which envisioned stock `sp.integrate` at ~5-10s/sample). It is mathematically correct and is the apparent reason the run fits the 600s cap. Not a finding — flagged for visibility so a future maintainer knows this helper, not stock SymPy, performs the unit-interval integrals.
- All 5 `PI_SAMPLES` rationals {1/2, 1, 3/2, 2, 5/3} match delta3's spec exactly; none dropped (delta3 floor was 4).

## Verdict justification

The sole finding F1 — that the SymPy engine verified exact symbolic identities only as a single-point mpmath numeric spot-check and never derived deliverable #1 — is resolved against the agreed final target. The SymPy engine is now an exact, tolerance-free, 5-rational-sample cross-check symbolic in `r1,r2,A_T,B_T,gprime`, with M1 supplying the previously-missing deliverable-#1 derivation, while the unchanged Mathematica engine carries the full all-`Pi_star` symbolic proof. Both engines exit 0 with all checks `= 0`/PASS, outputs are fresh, and the diff's one collateral edit (a custom exact integrator) is mathematically sound and non-tautological. Verdict: verified.
