# Directive pathA_C0g-diag — settle fold vs ill-conditioning vs sonic-horizon (cheap diagnostic battery, on existing states)

**Status:** READY (Claude-authored 2026-06-20 from the LOCKED spec `reports/pathA_throat_solver_literature_synthesis.md`
§6, produced by the 5-lane online literature research → Claude+Codex consult [§5, code-grounded corrections] → GLM 5.2
review [§6, the authoritative battery]). **Codex design-review of THIS directive = SOUND-WITH-FIXES → all 6 required fixes
+ anchor corrections + new fail conditions APPLIED** (1: project `F` via LEFT vectors `u`, not right-state subspaces;
2: gauge machinery is NOT turnkey — build a wrapper from the named C0c/C0d/C0e helpers; 3: gauge-projected spectrum =
`SVD(J·Q_perp)` NOT `SVD((I−P_G)J)`; 4: executable numeric thresholds + INCONCLUSIVE bands; 5: add the best-α trend across
converged τ; 6: cap scipy concretely; + anchor fixes + split-scripts + C0f2 provenance re-verify). **Codex CONFIRM-PASS
after fixes = SOUND-AS-IS** ("all six required fixes, anchor corrections, new fail conditions, premise halt behavior, and
scipy caps present and consistent; no new contradiction, dangling reference, or scope creep"). Logs
`_scratch/codex_c0g_diag_directive_{review,confirmpass}.log` (gitignored). **READY — awaiting the user EXECUTION gate.** Follows `pathA_C0f2` (which
reframed the τ≈0.029 "wall" as a crippled-defaults artifact + a REAL stall just past τ=0.029125). Step of option C, task
#78/#79.

**Date:** 2026-06-20
**Owner:** Codex (codes + iterates until exit 0; reuses the C0/C0c/C0e/C0f SVD + gauge-projector + JVP machinery; adds
bounded diagnostic + reporting code only). Claude reviews after (fidelity + adversarial agents).

**Trigger.** C0f2 established (fidelity-verified): the default-config crawl converges cleanly to τ=0.029125, then a REAL
stall remains just past it (best line-search step ~0.05%; α=1 overshoots ~13 orders; cond(J)≈1e13, σ_min≈1e-11). The
near-null subspace is gauge-dominated (C0c: U(1) phase mode exact-null; C0f/C0f2: modes 1/3/4 gauge by 1−P_G) with ONE
transverse candidate (mode 2, ~17% non-gradient). The OPEN question — **genuine fold/turning-point vs curable
ill-conditioning vs a physical sonic/horizon critical point** — gates which C0g machinery to build (pseudo-arclength for a
fold; PTC/LM + Sobolev preconditioner for conditioning). The literature converges on a CHEAP battery that settles this on
the EXISTING converged states with no production solver and no long run. **GLM's load-bearing catch:** the
"mode 2 DRIVES the stall" premise was only ever measured at the OLD crippled stall τ=0.02899 (C0e-0,
`near_null_component_fraction=1.8e-21`); the ACTUAL stall is a DIFFERENT state at τ=0.0290625, and the near-null direction
can rotate to align with the residual as τ↓ — so the premise MUST be re-validated FIRST (Step 0).

## Scope (this step DIAGNOSES; it does NOT build the fix)
C0g-diag READS the EXISTING converged states (τ = 0.03, 0.0295, 0.02925, 0.029125) + the C0f2 stalled state at
τ=0.0290625, and runs bounded analysis + ONE capped solver PROBE (scipy). It MUST NOT implement any C0g FIX —
**NO gauge-fixing/pinning in the solver, NO pseudo-arclength, NO trust-region/dogleg/LM build, NO PTC, NO Sobolev
preconditioner, NO analytic/sparse-assembly rewrite** (all DEFERRED to the gated C0g build, shaped by THIS verdict). It
MUST NOT modify the solver logic, the physical residual `patha_closed_branch_residual` (`coupled_branch.py:512`), the
faithful PDE operators, frozen physics, or `physical_export_permitted`; depth continuation stays **τ-only**.

**Allowed (read-only analysis, NOT a solver change):**
- **Gauge PROJECTION for the spectral tests — build a SMALL WRAPPER from named existing helpers (NOT turnkey).** There is
  NO `g_phase` symbol and NO ready-made gauge-projected-J API. Build the gauge generator matrix `G` from the EXISTING
  ANALYTIC helpers — C0c phase generator `_c0c_generators_for_state` (`patha_c0_conditioning_spike.py:3848`); C0d
  A-sector scalar-gradient basis `_c0d_scalar_gradient_matrix` / `_c0d_build_gauge_subspace` (`:2782` / `:2842`); C0e
  coupled local gauge basis `_c0e_coupled_gauge_matrix` (`:3118`); projection-fraction helper `_c0d_projection_fraction`
  (`:2891`) — then form the gauge-COMPLEMENT spectrum the CORRECT way (see the construction key-fact below), for ANALYSIS
  ONLY. This is new projection-aware spectral plumbing assembled from existing pieces, NOT a turnkey reuse. Projecting an
  analyzed matrix ≠ fixing the gauge in the solver.
- **ONE capped scipy PROBE** (Step 5) — `scipy.optimize.root(method='lm'|'hybr')` from the stalled state as a bounded
  globalization probe; judged ONLY by the ORIGINAL `patha_closed_branch_residual`; with a HARD `maxfev` + wall-time cap
  (Step 5 thresholds). It is a measurement, NOT the LM build, and MUST run in its OWN ≤600s script (not bundled with the
  spectral work).

**Single Arbiter Principle holds:** the ORIGINAL unmodified `patha_closed_branch_residual` is the SOLE arbiter of
convergence / progress / overlap; every metric is judged on the original ‖F‖, never a scaled/floored/least-squares
surrogate. **Genuine continuation:** any state recomputed here is warm-started from the prior converged state with
`prefer_existing_b2c_background_predictor=False` — NO cold-loading of pre-existing solutions (the C0-v2 sin). CPU;
`timeout 600` per script (split if needed; the capped scipy probe MUST fit well under the cap → else NOT_MEASURED, never
raise the cap); standalone `python3`; NO commentary `python3 -c`; YAML/markdown human output, JSON only for machine
artifacts; chunk-1a/1b/1c gates must still pass.

## Key facts (anchors VERIFIED in the Codex design-review; Codex re-confirm at execution)
- **State packing** `[psi_real, psi_imag, a0, ar, aw, r0, mu]` (pack `coupled_branch.py:137`, unpack `:166-174`); arbiter
  `patha_closed_branch_residual` `:512`; **gauge-covariant current** function starts ~`:327`, the `−(q/m)Aρ` terms at
  `:339` and `:342`; coupled matter uses the quintic enthalpy at `:385`. Quintic enthalpy `(5K/4)ρ⁴`, `ρ=|ψ|²`
  (`physics.py:46`) ⇒ residual is **C^∞** on the positive-R0 domain. The explicit sparse J
  (`assemble_closed_coupled_colored_sparse_jacobian`) is in **`preconditioners.py:476`** (NOT coupled_branch) ⇒
  **transpose-JVP is NOT needed** (use the explicit J for dense SVD / left-null / bordered work). `F_τ` = centered FD in
  τ at fixed x using the same residual + s_σ provider (with a stepsize-stability check — Step 2). grad closures
  `operators.py:254`/`:270`; the curl-like `F_rw = ∂_r a_w − ∂_w a_r` is `:431` (`:429` is the divergence).
- **Gauge-projected spectrum — the CORRECT construction (REQUIRED).** Do **NOT** dense-SVD `(I−P_G)J` in full
  coordinates: that injects ARTIFICIAL gauge-zero singular values and fabricates near-null modes. Instead: build the
  gauge generator matrix `G` from the ANALYTIC generators ONLY (the named C0c/C0d/C0e helpers); QR (or SVD) `G` to an
  orthonormal basis; form the complement basis `Q_perp` (columns spanning the gauge-orthogonal subspace); then dense-SVD
  the REDUCED operator `J·Q_perp` with DECLARED row/column scaling. Lift the right singular vector `v = Q_perp·v̂`; take
  the left-null vector `w` as the **left** singular vector of `J·Q_perp` at σ_min.
- **Right vs left subspaces (linear-algebra correctness).** The gauge generators span a RIGHT-state subspace: the
  state-space Newton step `δ` decomposes onto {gauge, near-null right-vectors `v`, complement}. The RESIDUAL `F` lives in
  the OUTPUT space and must **NOT** be projected onto right-state subspaces — decompose `F` using the **LEFT** singular
  vectors `u_i` (or project onto the IMAGES `J·subspace`).
- **The principal-symbol type-change test is a DEAD END (GLM, accepted):** in ψ-form the stationary operator is ALWAYS
  elliptic (principal symbol = the Schrödinger Laplacian ħ²|ξ|²/2m > 0; gauge terms lower-order; Maxwell
  curl-curl+grad-div positive-definite). The hydrodynamic elliptic→hyperbolic change at M=1 lives only in the
  Madelung frame and is SMEARED by quantum pressure. ⇒ Mach = physical CONTEXT (why a fold might exist); the Jacobian
  test (wᵀF_τ / σ_min²) = the DISCRIMINANT (whether one exists). Do NOT build a symbol/type-change test.
- **A simple fold ⟺ the bordered J_b is nonsingular ⟺ wᵀF_τ ≠ 0** (the transversality/regularity condition; a
  bifurcation has wᵀF_τ=0). This requires the gauge-projected J to have a 1-D null space — hence gauge projection FIRST.
- **No det/τ̇ SIGN-FLIP is observable from our states** — they are all on ONE side of the putative fold. Route on the
  scaled/gauge-projected σ_min(τ) trend, σ_min²(τ) scaling, tangent `|τ̇|→0`, and `‖x_τ‖` blow-up — NOT on raw det.

## Work items (cheapest-decisive-first; Step 0 is a HALT gate)

### Step 0 — framing check at the ACTUAL stall τ=0.0290625 (the PREMISE GATE)
Obtain the stalled iterate at τ=0.0290625 by GENUINE warm-start from the converged τ=0.029125 state + one τ-step
(`prefer_existing=False`); OR reuse the C0f2 saved stalled state ONLY after RE-VERIFYING its provenance from the C0f2 JSON
(`software/stage1_solver/runs/pathA_C0f2_timing_rerun/…`: it is attempt 004, `source: previous_c0_converged_state`,
`prefer_existing_b2c_background_predictor: false`, `used_existing_b2c: false`, saved merit matrix on the ORIGINAL/unscaled
residual ⇒ a genuine stalled iterate, NOT a re-verified old solution). Assemble the original-residual Jacobian + the clean
Newton step `Jδ=−F` (report the linear rel-resid). Build the gauge-complement basis `Q_perp` (Key-facts construction) and
the tracked near-null right-vectors `{v_i}` (= right singular vectors of `J·Q_perp` at the smallest σ).
- **Decompose the state-space step `δ`** onto: the **gauge** subspace (`G`), the **near-null** right-vectors `{v_i}`, and
  the **complement**. Report `near_null_component_fraction`, `gauge_component_fraction`,
  `transverse(mode-2)_component_fraction`.
- **Decompose the residual `F` via LEFT vectors** (NOT right-state subspaces): report `‖u_iᵀF‖/‖F‖` for the tracked left
  singular vectors `u_i` of `J·Q_perp`, and the residual's projection onto the gauge IMAGES `J·G`.
- Report the α=1 predicted-vs-actual residual gap and the **best-α merit reduction**.
**PREMISE GATE (executable):** with `f_nn = near_null_component_fraction` (Step δ-decomposition): **PREMISE_FAILED if
`f_nn < 1e-2`** (the stall is NOT near-null amplification, as at the stale τ=0.02899 where it was 1.8e-21) → HALT the
battery; **premise holds if `f_nn > 1e-1`**; **INCONCLUSIVE band `[1e-2, 1e-1]`** → verdict `DIAGNOSTIC_INCOMPLETE` unless
the other Step-0 evidence (α=1 overshoot + best-α near-zero) clearly corroborates near-null amplification (state which).
On PREMISE_FAILED the report/JSON MUST contain NO Step 1–7 measurements except the literal `SKIPPED_PREMISE_FAILED`.

### Step 0b — best-α descent TREND across the converged τ states (GLM §6.1)
Run the same best-α merit sweep at the converged states (τ=0.029125, 0.02925, 0.0295, 0.03), not just the stall. Report
the best-α and its %reduction per τ. **Degrading monotonically toward the stall** (best-α reduction → 0, direction →
tangent) ⇒ fold-support; **roughly constant across τ** ⇒ fixed-conditioning support. Context for the Step-8 verdict.

### Step 1 — gauge-complement SVD at the converged states
At τ = 0.03, 0.0295, 0.02925, 0.029125 (and the τ=0.0290625 stalled state): dense SVD of the REDUCED operator `J·Q_perp`
(Key-facts construction — NOT `(I−P_G)J`), with declared row/column scaling. Report `σ_min(τ)`, `σ_min/σ_max`, the
left-null vector `w(τ)` (left singular vector at σ_min), the lifted right-null vector `v(τ)=Q_perp·v̂`, and per-mode
tracking by vector overlap (indices swap — track, don't index). Report each tracked mode's lane fractions +
center-of-energy. (Dense SVD at n=1297 is ~14–15 s/state per the C0c/C0e validation; assemble each J fresh ~40 s — SPLIT
across ≤600 s scripts.)

### Step 2 — wᵀF_τ (the PRIMARY fold test)
At each state: `F_τ` = centered FD in τ, `F_τ ≈ (F(x,τ+h) − F(x,τ−h))/(2h)`, at fixed x. **Run a stepsize-stability
check** (e.g. `h ∈ {1e-4, 1e-5, 1e-6}·τ`; require `cosθ` stable to ≤10% across two h's, else report the h-sensitivity and
downgrade to NOT_MEASURED) — a clean `cosθ` off a single unstable h is a stepsize artifact, not a fold. Then the
SCALE-FREE transversality measure `cosθ = |wᵀF_τ| / (‖w‖·‖F_τ‖)` with `w` the gauge-complement left-null vector from
Step 1 (project `F_τ` consistently into the left/output space — Key-facts right-vs-left). **Executable gate:
`cosθ > 0.1` ⇒ FOLD-support; `cosθ < 1e-2` ⇒ BIFURCATION (F_τ ⟂ w, in the range of J); INCONCLUSIVE band `[1e-2, 0.1]`.**
Report `cosθ` (per h), `wᵀF_τ`, `‖F_τ‖` per state.

### Step 3 — σ_min²(τ) scaling fit (free, from Step 1)
Fit the gauge-projected `σ_min²(τ)` vs τ over the converged states. **Linear (nonzero slope) ⇒ simple fold** (σ_min ∝
√|τ−τ_fold|; the zero-crossing predicts `τ_fold`). **Quadratic / flat ⇒ bifurcation or conditioning** (σ_min ∝
|τ−τ|). Report the linear-fit slope + R², the quadratic-fit comparison, the predicted `τ_fold`, AND the monotonicity of
`σ_min(τ)` — flag the fit UNRELIABLE if σ_min is non-monotone across the 4 points.

### Step 4 — gauge-covariant Mach map (physical CONTEXT only)
Compute the Mach map on cell centers from the converged + stalled states using the gauge-covariant current:
`j_r=(ħ/m)(ψ_r ∂_r ψ_i − ψ_i ∂_r ψ_r) − (q/m)A_r ρ` (and `j_w` analogously), `v=j/max(ρ,floor)`,
`c_s = √(5K/m)·ρ²`, `M = |v|/c_s`. Report `M_full_max`, `M_w_max`, `rho_at_max`, **`j_at_max`**, `(r*,w*,r*/R0(w*))`,
stability under a density-floor sweep, and the `M_max(τ)` trend. **Context-only — NOT a verdict gate.** Executable
sonic-context labels: **sonic-support if `0.8 ≤ M_max ≤ 1.2` AND `M_max(τ)` increases monotonically toward that band as
τ↓ AND the current is non-vanishing at the argmax (`|j_at_max| > 1e-3 · median|j|` over non-vacuum cells)**; **no-sonic
if `M_max < 0.7`**; otherwise INCONCLUSIVE. If `|j_at_max|` is below the current threshold (empty core, ~zero current) the
sonic hypothesis is NOT represented in these stationary states regardless of M (report explicitly). Do NOT assume the
sonic point sits on R0 (R0 enters the confinement potential).

### Step 5 — capped scipy LM/hybr probe + branch-overlap (the both-sides signal; OWN ≤600s script)
From the τ=0.0290625 stalled state, run `scipy.optimize.root` in its OWN script (NOT bundled with the spectral work).
**Concrete caps (defaults are NOT acceptable at n=1297 — one dense FD Jacobian ≈ n+1 ≈ 1298 residual evals):** method
order = `'lm'` first, then `'hybr'` if wall-time remains; per-method `maxfev = 4·(n+1)` (≈ 4 Jacobian builds) AND a hard
wall timer of **≤ 300 s per method** (a Python-level elapsed check that aborts the call); on exceeding either cap →
`NOT_MEASURED` for that method (never raise the `timeout 600` cap). Judge ONLY by the ORIGINAL
`patha_closed_branch_residual`: progress = `Linf ≤ 1e-6` OR original-residual `L2 ↓ ≥ 10×`. If it converges, compute the
**branch-overlap** in the SAME scaled-coordinate metric as the spectral tests (NOT a raw normalized dot product across
mixed units/gauge), and report **lane-wise** overlaps (ψ / A / r0 / μ) with the pre-stall converged solution:
**overlap > 0.99 = SAME branch** (conditioning / fold deeper); **< 0.5 = POST-FOLD branch** (fold confirmed); the band
`[0.5, 0.99]` = INCONCLUSIVE. Report the caps used, iters/fev, the original-residual trajectory, the lane-wise overlaps,
and a "no evidence" (NOT "fold proven") result if it neither progresses nor finds a distinct branch.

### Step 6 — bordered cond on the gauge complement: cond(J·Q_perp) vs cond(Jb) (CONFIRMATORY)
Border the REDUCED operator `J·Q_perp` (Step 1) with the Step-2 tangent/null data and compare scaled
`cond(J·Q_perp)` vs `cond(Jb)`. **Fold-support if `cond(J·Q_perp) > 1e10` and `cond(Jb) < 1e8`, or
`σ_min(Jb)/σ_min(J·Q_perp) > 1e4`.** Confirms Step 2; report both conds + the ratio. (No absolute O(1) thresholds —
scale-relative only.)

### Step 7 — commutator ‖curl_h·G‖ + per-mode P_G overlap (MODE-2 characterization)
Using the SAME centered/one-sided closures (`operators.py:254`/`:270` grad, `:431` curl `F_rw=∂_r a_w−∂_w a_r`): report
`‖C·G‖/‖G‖`, the boundary-row fraction, per-tracked-mode `P_G` (re-confirm 1/3/4 gauge, mode 2 ≈ transverse), and the
residual curl of each mode. **Mode-2 characterization (stencil-artifact vs real transverse), NOT fold discrimination.**
`P_G` MUST be built from the ANALYTIC gradient generators `G` — NEVER from the SVD modes themselves (tautology) — and use
the unbiased dimensionless `1−P_G` energy fraction (NOT the k-biased curl ratio that sank C0e), with a mixed control.

**Run order + reporting gate:** run Step 0 (+ 0b) FIRST. **If PREMISE_FAILED, HALT** — the report/JSON contains NO Step
1–7 measurements except the literal `SKIPPED_PREMISE_FAILED`. Otherwise run the spectral/Mach work (Steps 1–4, 6, 7) and
the capped Step 5 in SEPARATE ≤600 s scripts (assemble ≈40 s/J × ~5 states + SVD ≈15 s/state will not fit one script;
the scipy probe MUST be its own script). Steps 0–3 are premise-deciding — surface them prominently so an early decisive
result (PREMISE_FAILED, or a clean FOLD/CONDITIONING agreement) is visible without wading through 4–7.

## Verdict (Step 8) — exactly one, with falsifiable numeric `verdict_support`; require AGREEMENT, do NOT force a clean call
**Gray-zone rule (binding):** a clean `FOLD_CONFIRMED` or `CONDITIONING` verdict may be declared ONLY if NO primary
evidence (Step-0 `f_nn`, Step-2 `cosθ`, Step-3 `σ_min²` fit, Step-5 overlap) sits in its explicit INCONCLUSIVE band. If
any primary item is in a gray band, or the items disagree, the verdict is **MIXED / INCONCLUSIVE** (or
`DIAGNOSTIC_INCOMPLETE` if items are NOT_MEASURED) — never a forced clean call.
  - **PREMISE_FAILED** — Step-0 `f_nn < 1e-2` at the real stall. ⇒ HALT; the mode-2/gauge/fold framing is wrong;
    recommend re-diagnosing what actually drives the stall before any C0g build.
  - **FOLD_CONFIRMED** — Step-2 `cosθ > 0.1` AND Step-3 `σ_min²` linear (good R², monotone σ_min) AND (Step-5 overlap
    `< 0.5` OR Step-6 `cond(Jb) < 1e8` with `cond(J·Q_perp) > 1e10`). Sub-label:
      - **SONIC_FOLD** if ALSO Step-4 sonic-support (`0.8 ≤ M_max ≤ 1.2`, monotone approach, `|j_at_max|` above the
        current threshold) at a location mode-2 energy overlaps ⇒ recommend C0g = gauge-fixed pseudo-arclength (cheapest)
        or shoot-from-sonic + L'Hôpital (most faithful).
      - **NON_SONIC_FOLD** if Step-4 `M_max < 0.7` ⇒ recommend C0g = gauge-fixed pseudo-arclength.
  - **CONDITIONING** — Step-2 `cosθ < 1e-2` OR Step-3 `σ_min²` quadratic/flat, AND the Step-5 LM probe gets ≥10×
    original-residual drop / same branch (overlap > 0.99). ⇒ recommend C0g = shifted-Newton LM (`J+λI`, λ∈[1,100] anneal)
    and/or PTC+SER + inverse-Hamiltonian (Sobolev) preconditioner on the matter block.
  - **MIXED / INCONCLUSIVE** — the primary tests DISAGREE (e.g. `cosθ>0.1` but `σ_min²` flat), OR any primary item is in
    its gray band, OR gauge/stencil modes dominate while mode 2 stays transverse with Mach inconclusive. ⇒ recommend
    gauge-fix + LM/PTC first, THEN rerun the spectral/Mach battery before committing to pseudo-arclength. (Prefer THIS
    over a forced clean verdict — the user distrusts clean calls.)
  - **DIAGNOSTIC_INCOMPLETE** — required evidence NOT_MEASURED (scipy probe hit the cap; `F_τ` h-unstable; σ_min
    non-monotone so the σ_min² fit is unreliable; too few converged states for a tangent/fit).

`verdict_support` must include: the Step-0 decomposition table; the per-state SVD table (σ_min, σ_min/σ_max, mode lanes);
the wᵀF_τ `cosθ` table; the σ_min² fit (slope, R², τ_fold, monotonicity flag); the Mach table (incl. `j_at_max`); the
capped-LM trajectory + overlap; the bordered cond(J)/cond(Jb)/ratio; and the commutator + per-mode `1−P_G`.

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. Step 0 obtains the τ=0.0290625 stalled iterate by GENUINE warm-start (RE-VERIFY the C0f2 provenance from JSON, or
   recompute; NO cold-load), decomposes the state-step `δ` onto gauge/near-null/complement and the residual `F` via LEFT
   vectors `u` (NOT right-state subspaces), both on the ORIGINAL residual, and applies the executable PREMISE GATE
   (`f_nn < 1e-2` ⇒ HALT with ONLY `SKIPPED_PREMISE_FAILED` in the JSON). Step 0b reports the best-α trend across τ.
2. The spectral tests (Steps 1, 2, 6) use the CORRECT gauge-complement construction — `G` from the named ANALYTIC
   generators (`_c0c_generators_for_state`, `_c0d_scalar_gradient_matrix`/`_c0d_build_gauge_subspace`,
   `_c0e_coupled_gauge_matrix`), QR→`Q_perp`, dense-SVD `J·Q_perp` with declared scaling — NOT `SVD((I−P_G)J)`; `w` is the
   left singular vector of `J·Q_perp` at σ_min.
3. wᵀF_τ is the SCALE-FREE `cosθ = |wᵀF_τ|/(‖w‖‖F_τ‖)` with `F_τ` a real centered τ-FD of the ORIGINAL residual PLUS a
   stepsize-stability check across ≥2 `h`; σ_min²(τ) fit reports slope + R² + τ_fold + the monotonicity flag.
4. The Mach map uses the GAUGE-COVARIANT current (`coupled_branch.py:339`/`:342`, the −(q/m)Aρ terms present),
   `c_s=√(5K/m)ρ²`, reports `j_at_max`, applies the executable sonic-context thresholds, and is labeled CONTEXT-ONLY
   (not a verdict gate).
5. The scipy probe runs in its OWN ≤600 s script with the CONCRETE caps (`maxfev=4(n+1)`, ≤300 s/method wall timer,
   `'lm'` then `'hybr'`, per-method `NOT_MEASURED` on cap), judged ONLY by the ORIGINAL `patha_closed_branch_residual`,
   with the SCALED lane-wise branch-overlap gate; a non-progressing result is "no evidence," not "fold proven."
6. Exactly one Step-8 verdict obeying the GRAY-ZONE rule + AGREEMENT of the primary tests (MIXED/INCONCLUSIVE on
   disagreement or any gray-band item), with the full falsifiable numeric `verdict_support`.
7. NO C0g FIX implemented (diff): no gauge-fix/pin in the solver, no pseudo-arclength, no trust-region/dogleg/LM build, no
   PTC, no Sobolev preconditioner, no analytic/sparse-assembly rewrite. Faithful operators + frozen physics + export guard
   + SOLVER LOGIC UNTOUCHED (diff); depth continuation τ-only; chunk gates pass; report + machine JSON emitted.

**Fail conditions:** cold-loading the stalled state (the C0-v2 sin) or skipping the C0f2 provenance re-verify; skipping or
soft-pedaling the Step-0 premise gate (or emitting Step 1–7 measurements on PREMISE_FAILED); **dense-SVD'ing `(I−P_G)J` in
full coordinates instead of `J·Q_perp` (injects fake gauge-zero singulars); building the gauge basis from the SVD modes
themselves (tautology); projecting the residual `F` onto right-state subspaces (use LEFT vectors); a gauge-basis
rank/control failure left unreported; a `cosθ` declared without the `F_τ` stepsize-stability check; the scipy probe run
uncapped or sharing the spectral script**; running the spectral tests on the RAW (un-projected) J (w would be a gauge mode
→ false "not a fold"); routing on raw det / an absolute O(1) bordered threshold / the k-biased curl ratio; judging the
scipy probe by its least-squares
merit instead of the original residual; declaring a clean FOLD or CONDITIONING verdict when the independent tests disagree
(must be MIXED); implementing ANY C0g fix; changing a frozen physics parameter or the operators/export/solver-logic;
masking NOT_MEASURED; raising the timeout cap.

## Out of scope (the gated C0g BUILD — shaped by THIS verdict)
- **Unconditional (the C0g build, not here):** gauge-fix PATH-ONLY (far-field ψ phase pin + A-sector gauge handling in
  solver/preconditioner coords — NEVER the frozen residual); analytic/sparse Jacobian assembly AFTER the color audit
  (handle the dense μ/mass/wall rows; 253 colors = deterministic radius-3 coloring, not a bug).
- **Conditioning branch:** shifted-Newton LM / PTC+SER + inverse-Hamiltonian (Sobolev) preconditioner.
- **Fold branch:** gauge-fixed pseudo-arclength continuation; optional shoot-from-sonic + L'Hôpital.
- Then (only after the crawl unblocks toward the physics-target depth): constitutive family → calibrated branch →
  multi-knob calibrate-predict (R0/J/W → anchor → SURPLUS) → `pathA_22`.

## Review (orchestrator, after Codex)
Fidelity agent (term-by-term, code-vs-spec): is the Step-0 stalled state a GENUINE warm-start (not cold-loaded — show the
re-verified C0f2 provenance + `prefer_existing=False`); is `δ` decomposed in state space and `F` via LEFT vectors, both on
the ORIGINAL residual; are Steps 1/2/6 the CORRECT `J·Q_perp` construction (G from the named analytic generators
`_c0c_generators_for_state`/`_c0d_*`/`_c0e_coupled_gauge_matrix`, QR→Q_perp, NOT `(I−P_G)J`, NOT G from SVD modes); is
`cosθ` the scale-free transversality with a real τ-FD `F_τ` + the stepsize-stability check; is the Mach current
gauge-COVARIANT with the −(q/m)Aρ terms and `c_s=√(5K/m)ρ²` (transliteration-fidelity: match `coupled_branch.py:339`/`:342`
+ `physics.py:46`); is the scipy probe in its own script, capped, judged on the original residual; are
operators/frozen/export/SOLVER-LOGIC untouched (diff) and NO C0g fix implemented? Adversarial agent: is the verdict
EARNED — is PREMISE_FAILED honestly applied (not skipped to keep the battery alive); is a FOLD/CONDITIONING call backed by
AGREEMENT of the independent tests (not one test cherry-picked, not a clean call forced over disagreement); is the
one-sided σ_min² fit trustworthy (monotone) or over-read; is the LM "no evidence" not spun into "fold proven"; any
can't-fail-gate / tautological projector / hardcoded threshold? Then gate the C0g build branch (fold → pseudo-arclength;
conditioning → LM/PTC+Sobolev; mixed → gauge-fix+LM then re-run).
