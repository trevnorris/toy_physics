# Path-A C0g Diagnostic Battery Final Report (Steps 0-8)

Step-8 verdict: **MIXED/INCONCLUSIVE**

Recommended C0g build branch: Gauge-fix plus LM/PTC first, then rerun the battery before pseudo-arclength.

The original `patha_closed_branch_residual` remains the sole convergence/progress arbiter. Steps 4-7 are diagnostic only; no C0g fix is implemented here.

## Provenance

- C0f2 stalled-state provenance: `GENUINE_WARM_START_NOT_COLD_LOADED`
- prefer_existing_b2c_background_predictor: `False` from `crawl_config.prefer_existing_b2c_background_predictor`
- source: `previous_c0_converged_state`; used_existing_b2c: `False`

## Step 0 Premise Gate

Call: `PREMISE_HOLDS` with `f_nn=9.999999e-01`.

| tau | f_nn | gauge_fraction | mode2_fraction | gauge_image_F_fraction | alpha1_actual_ratio | best_alpha | best_reduction_percent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2.906250e-02 | 9.999999e-01 | 4.987886e-27 | 4.691096e-06 | 1.612607e-19 | 4.283280e+00 | 7.812500e-03 | 4.902981e-02 |

## Step 0b Best-Alpha Trend

Trend call: `DEGRADING_TOWARD_STALL`.

| tau | best_alpha | best_actual_l2_ratio | best_percent_reduction | alpha1_actual_l2_ratio | alpha1_predicted_l2_ratio |
| --- | --- | --- | --- | --- | --- |
| 2.912500e-02 | 1.562500e-02 | 9.448496e-01 | 5.515044e+00 | 1.733119e+01 | 6.300907e-12 |
| 2.925000e-02 | 5.000000e-01 | 5.688369e-01 | 4.311631e+01 | 8.698908e-01 | 3.122009e-13 |
| 2.950000e-02 | 1.000000e+00 | 3.437082e-01 | 6.562918e+01 | 3.437082e-01 | 1.756710e-13 |
| 3.000000e-02 | 1.000000e+00 | 2.469664e-01 | 7.530336e+01 | 2.469664e-01 | 9.305143e-14 |

## Step 1 Gauge-Complement SVD

| tau | status | sigma_min | sigma_min/sigma_max | gauge_rank | mode1_lane |
| --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | MEASURED | 2.245655e-03 | 2.286600e-04 | 511 | r0 |
| 2.950000e-02 | MEASURED | 1.502323e-03 | 1.529715e-04 | 511 | r0 |
| 2.925000e-02 | MEASURED | 8.301718e-04 | 8.453086e-05 | 511 | r0 |
| 2.912500e-02 | MEASURED | 4.267758e-05 | 4.345574e-06 | 511 | r0 |
| 2.906250e-02 | MEASURED | 2.448979e-04 | 2.493632e-05 | 511 | r0 |

## Step 2 wT F_tau

| tau | status | cos_theta | call | stability |
| --- | --- | --- | --- | --- |
| 3.000000e-02 | MEASURED | 5.138510e-01 | FOLD_SUPPORT | STABLE |
| 2.950000e-02 | MEASURED | 5.572308e-01 | FOLD_SUPPORT | STABLE |
| 2.925000e-02 | MEASURED | 5.902699e-01 | FOLD_SUPPORT | STABLE |
| 2.912500e-02 | MEASURED | 6.235099e-01 | FOLD_SUPPORT | STABLE |
| 2.906250e-02 | MEASURED | 6.348990e-01 | FOLD_SUPPORT | STABLE |

## Step 3 Sigma-Min Squared Fit

- call: `LINEAR_MONOTONE_FOLD_SUPPORT`
- linear slope: `5.783863e-03`
- linear R2: `9.994113e-01`
- tau_fold zero crossing: `2.912335e-02`
- quadratic R2: `9.997300e-01`
- monotone: `True`

## Step 4 Mach Context

Context label: `SONIC_HYPOTHESIS_NOT_REPRESENTED`; represented: `False`.
Monotone M_max approach as tau decreases: `True`.

| tau | M_full_max | M_w_max | rho_at_max | j_at_max | r_star | w_star | r_star_over_R0 | current_represented | floor_stability |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 0.000000e+00 | 0.000000e+00 | 1.229132e-02 | 0.000000e+00 | 5.000000e-02 | 5.781250e-02 | 5.076365e-02 | false | STABLE |
| 2.950000e-02 | 0.000000e+00 | 0.000000e+00 | 1.221981e-02 | 0.000000e+00 | 5.000000e-02 | 5.781250e-02 | 5.072874e-02 | false | STABLE |
| 2.925000e-02 | 0.000000e+00 | 0.000000e+00 | 1.216943e-02 | 0.000000e+00 | 5.000000e-02 | 5.781250e-02 | 5.070026e-02 | false | STABLE |
| 2.912500e-02 | 0.000000e+00 | 0.000000e+00 | 1.212110e-02 | 0.000000e+00 | 5.000000e-02 | 5.781250e-02 | 5.066921e-02 | false | STABLE |
| 2.906250e-02 | 0.000000e+00 | 0.000000e+00 | 1.210311e-02 | 0.000000e+00 | 5.000000e-02 | 5.781250e-02 | 5.065797e-02 | false | STABLE |

## Step 5 Capped Scipy Probe

Call: `NO_PROGRESS_NO_EVIDENCE`; status: `NOT_MEASURED`.
Caps used: maxfev `5192` (`4*(n+1)`), method wall `285.0` s, global wall `600.0` s.

| method | status | nfev | elapsed_seconds | best_l2_ratio | best_linf | progress_call | abort_reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lm | NOT_MEASURED | 5192 | 4.699241e+01 | 5.397528e-01 | 2.434638e-05 | NO_PROGRESS_NO_EVIDENCE | maxfev_cap |
| hybr | NOT_MEASURED | 5192 | 4.588796e+01 | 9.930337e-01 | 2.474210e-05 | NO_PROGRESS_NO_EVIDENCE | maxfev_cap |

Branch-overlap aggregate: `NOT_COMPUTED_NO_PROGRESS`.

| lane | scaled_overlap | call | candidate_delta_scaled_norm | prestall_delta_scaled_norm |
| --- | --- | --- | --- | --- |

A non-progressing capped probe is interpreted as no evidence, not fold proof.

## Step 6 Bordered Reduced Conditioning

Call: `NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD`.

| tau | status | cond_JQ_perp | cond_Jb | sigma_min_JQ_perp | sigma_min_Jb | ratio | call |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | MEASURED | 4.373305e+03 | 5.699423e+02 | 2.245655e-03 | 1.723145e-02 | 7.673240e+00 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.950000e-02 | MEASURED | 6.537164e+03 | 5.699420e+02 | 1.502323e-03 | 1.723146e-02 | 1.146988e+01 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.925000e-02 | MEASURED | 1.183000e+04 | 5.699417e+02 | 8.301718e-04 | 1.723146e-02 | 2.075650e+01 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |
| 2.912500e-02 | MEASURED | 2.301192e+05 | 5.699416e+02 | 4.267758e-05 | 1.723147e-02 | 4.037593e+02 | NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD |

## Step 7 Commutator And P_G

C0g gauge-complement tracked modes; mode 2 here is the a0 lane with sigma about 1.7e-2, not the r0 stall-driving minimum mode.

Commutator rows:

| tau | norm_CG_over_G | boundary_fraction | mixed_control_one_minus_P_G | mixed_control_curl_over_A |
| --- | --- | --- | --- | --- |
| 3.000000e-02 | 1.001686e-15 | 8.831918e-01 | 5.829003e-01 | 1.009363e+01 |
| 2.950000e-02 | 1.001686e-15 | 8.831918e-01 | 5.713494e-01 | 9.397264e+00 |
| 2.925000e-02 | 1.001686e-15 | 8.831918e-01 | 5.880325e-01 | 9.608546e+00 |
| 2.912500e-02 | 1.001686e-15 | 8.831918e-01 | 5.776887e-01 | 1.012324e+01 |
| 2.906250e-02 | 1.001686e-15 | 8.831918e-01 | 5.977188e-01 | 1.085252e+01 |

Per-mode `1-P_G` rows:

| tau | tracked_lane | ascending_rank | sigma | one_minus_P_G | P_G | curl_over_A | boundary_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000e-02 | 1 | 1 | 2.245655e-03 | 1.000000e+00 | 6.216658e-26 | 8.349897e+00 | 7.392434e-01 |
| 3.000000e-02 | 2 | 2 | 1.723152e-02 | 1.000000e+00 | 8.074042e-29 | 1.020920e+01 | 7.746370e-01 |
| 3.000000e-02 | 3 | 3 | 3.362335e-02 | 1.000000e+00 | 1.734813e-24 | 9.888601e+00 | 8.480159e-01 |
| 3.000000e-02 | 4 | 4 | 6.195948e-02 | 1.000000e+00 | 3.352760e-24 | 9.232153e+00 | 8.473621e-01 |
| 3.000000e-02 | 5 | 5 | 6.807263e-02 | 1.000000e+00 | 2.719566e-24 | 9.644617e+00 | 8.840342e-01 |
| 2.950000e-02 | 1 | 1 | 1.502323e-03 | 1.000000e+00 | 4.567242e-26 | 1.057973e+01 | 8.609067e-01 |
| 2.950000e-02 | 2 | 2 | 1.723152e-02 | 1.000000e+00 | 7.164336e-29 | 7.140771e+00 | 7.189182e-01 |
| 2.950000e-02 | 3 | 3 | 3.236622e-02 | 1.000000e+00 | 2.688817e-24 | 5.831982e+00 | 7.280060e-01 |
| 2.950000e-02 | 4 | 4 | 6.209394e-02 | 1.000000e+00 | 1.441133e-24 | 1.132127e+01 | 9.155540e-01 |
| 2.950000e-02 | 5 | 5 | 6.658556e-02 | 1.000000e+00 | 3.408841e-24 | 6.195849e+00 | 7.564511e-01 |
| 2.925000e-02 | 1 | 1 | 8.301718e-04 | 1.000000e+00 | 1.990203e-26 | 7.206434e+00 | 7.955257e-01 |
| 2.925000e-02 | 2 | 2 | 1.723152e-02 | 1.000000e+00 | 5.941138e-29 | 1.110514e+01 | 8.176556e-01 |
| 2.925000e-02 | 3 | 3 | 3.142533e-02 | 1.000000e+00 | 2.051503e-24 | 7.791855e+00 | 8.079037e-01 |
| 2.925000e-02 | 4 | 4 | 6.206744e-02 | 1.000000e+00 | 3.062406e-24 | 5.187191e+00 | 7.352392e-01 |
| 2.925000e-02 | 5 | 5 | 6.585785e-02 | 1.000000e+00 | 2.521891e-24 | 1.061374e+01 | 8.480389e-01 |
| 2.912500e-02 | 1 | 1 | 4.267758e-05 | 1.000000e+00 | 1.093785e-26 | 9.612163e+00 | 8.159216e-01 |
| 2.912500e-02 | 2 | 2 | 1.723151e-02 | 1.000000e+00 | 5.659374e-29 | 8.797871e+00 | 7.213283e-01 |
| 2.912500e-02 | 3 | 3 | 3.049417e-02 | 1.000000e+00 | 2.243021e-24 | 6.844563e+00 | 8.007722e-01 |
| 2.912500e-02 | 4 | 4 | 6.186087e-02 | 1.000000e+00 | 2.910725e-24 | 1.021346e+01 | 8.904766e-01 |
| 2.912500e-02 | 5 | 5 | 6.553625e-02 | 1.000000e+00 | 1.287672e-24 | 1.182555e+01 | 8.794100e-01 |
| 2.906250e-02 | 1 | 1 | 2.448979e-04 | 1.000000e+00 | 1.143200e-26 | 8.314462e+00 | 7.684710e-01 |
| 2.906250e-02 | 2 | 2 | 1.723151e-02 | 1.000000e+00 | 6.461278e-29 | 7.601729e+00 | 7.531209e-01 |
| 2.906250e-02 | 3 | 3 | 3.016670e-02 | 1.000000e+00 | 1.722129e-24 | 5.317101e+00 | 7.711767e-01 |
| 2.906250e-02 | 4 | 4 | 6.174833e-02 | 1.000000e+00 | 2.574017e-24 | 9.742013e+00 | 8.869321e-01 |
| 2.906250e-02 | 5 | 5 | 6.549392e-02 | 1.000000e+00 | 4.860511e-25 | 5.204555e+00 | 7.192031e-01 |

## Step 8 Verdict Support

- primary fold support: `True`
- primary conditioning support: `False`
- gray items: `[]`
- not measured: `[]`
- disagreements: `[]`
- Step-5 overlap call: `NOT_COMPUTED_NO_PROGRESS`
- Step-6 call: `NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD`
- gray-zone/agreement rule honored: `True`

## Post-hoc adjudication (orchestrator + dual review, 2026-06-20)

The machine Step-8 verdict above is `MIXED/INCONCLUSIVE` — the faithful output of the directive's literal gray-zone rule —
and is retained as-is. However, both post-execution reviews (fidelity = `FAITHFUL_WITH_MINOR_NOTES`; adversarial =
`SHOULD_BE_FOLD_LEANING_NOT_FLAT_MIXED`) concluded the body of evidence is **asymmetrically fold-leaning, not a symmetric
"could not decide."** Recorded so a future reader does not misread `MIXED` as "conditioning is plausible":

- **The stall is a simple FOLD (turning point) in the throat-radius (`r0`) continuation at τ_fold ≈ 0.0291233 — NOT
  conditioning, NOT sonic.** Five independent lines agree and conditioning is affirmatively refuted
  (`primary_fold_support=True`, `primary_conditioning_support=False`, no disagreements, no gray items): Step-0 premise
  HOLDS (f_nn≈1.0; the Newton step rides the `r0` near-null mode; gauge fraction 5e-27 ⇒ the gauge framing is refuted as
  the stall *driver*); Step-2 cosθ 0.51–0.64 stable (wᵀF_τ≠0 ⇒ simple fold, not bifurcation); Step-3 σ_min² linear
  R²=0.9994 → τ_fold; Step-6 bordered-cond TREND (cond(Jb) flat ≈570 while cond(J·Q_perp) → 2.3e5, ratio → 404,
  accelerating); Step-4 Mach: zero current (ψ effectively real, static empty-flow background) ⇒ sonic NOT TESTABLE here
  (→ NON_SONIC).
- **Why a clean `FOLD_CONFIRMED` was (correctly) withheld:** only the Step-6 *absolute* threshold (cond(J·Q_perp) > 1e10)
  was not met. That threshold is **self-defeating for a fold diagnosis**: since σ_min ∝ √(τ − τ_fold), reaching cond > 1e10
  requires sampling ~1e-15 from τ_fold, but Newton stalls ~1.6e-6 short *by the very definition of a fold*. (The directive
  states "No absolute O(1) thresholds" yet sets an absolute 1e10 — an internal inconsistency.) The **bordered-cond TREND**
  — cond(Jb) flat while cond(J·Q_perp) and the ratio accelerate toward τ_fold — is the correct, scale-relative
  discriminant, and it supports a simple fold.
- **The 511-dim gauge near-null subspace is a SEPARATE, independently-real conditioning issue** (it is *not* the stall
  driver — the Newton step's gauge component is ~5e-27) that must be fixed regardless.

**Recommended C0g build (both reviews converge):** (a) gauge-fix FIRST (path-only; also lets the crawl sample closer to
τ_fold); (b) re-run cheap Steps 1–3 + 6 on the gauge-fixed crawl nearer τ_fold → promote `MIXED` → `FOLD_CONFIRMED` if
σ_min² stays linear-monotone and cond(J·Q_perp) crosses ~1e8–1e10 with cond(Jb) flat; (c) gauge-fixed pseudo-arclength
continuation to round the fold. Skip the LM/PTC conditioning detour (the evidence disfavors it). See `decisions/13` §0 (4).

## Scope Confirmation

- NO C0g fix implemented: `True`
- Single arbiter residual: `stage1_solver.coupled_branch.patha_closed_branch_residual`
- Depth continuation: `tau_only`
- Solver logic touched: `False`
- Faithful operators touched: `False`
- Frozen physics touched: `False`
- Physical export guard touched: `False`

Complete machine JSON: `software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/pathA_C0g_diag_fold_vs_conditioning.json`.
