# Path-A throat-solver — consolidated findings + research brief (C0 → C0f2)

**Date:** 2026-06-20. **Purpose:** (1) the single self-contained record of everything the C0–C0f2 conditioning
investigation established; (2) a precise problem statement to drive an ONLINE LITERATURE SEARCH for how this class of
solve is done + what tricks apply; (3) the post-/compact resume anchor. Read this first on resume.

---

## 0. TL;DR
The Path-A "throat profile solve" (the bottleneck for the whole physics program) is a **parameter-continuation of a
near-singular nonlinear elliptic BVP**. After 6 diagnostic rounds we have a VERIFIED diagnosis:
- The original "wall at τ≈0.029" was **partly a crippled-solver config** — with code defaults the Newton continuation
  now produces GENUINE converged states down to **τ=0.029125** (Linf~1e-12), which the crippled run never had.
- **But a real stall remains just past 0.029125.** Plain Newton + Armijo line-search is **effectively stuck** (best
  step reduces the residual ~0.05%). Cause: the Jacobian is **near-singular** (σ_min≈1e-11 vs σ_max≈826) and the Newton
  step amplifies a near-null direction that has a **transverse, non-gauge component (mode 2)** → the full step
  overshoots by ~13 orders → line-search collapses.
- **A gauge-only fix is ruled out** (a pure gauge direction can't quadruple the gauge-invariant residual; confirmed mode
  2 is genuinely transverse). **A production linear solver is NOT indicated** (the linear solve is already clean to
  1e-12). This is a **Newton-globalization + possible-fold** problem.
- **Open question (one cheap measurement settles it):** is τ≈0.029 a genuine **fold / turning point** (~30% likely) or
  just **severe ill-conditioning** (~55%)? → the det(J)-sign / σ_min(τ) / branch-tangent probe.
- **Cost:** ~96–97% of each step is **Jacobian assembly** (colored finite-difference/JVP, ~1k–2.5k probes/build) — a
  ~100× algorithmic target (sparse/analytic assembly), but irrelevant until the stall is cured.

---

## 1. The numerical problem (precise — this is what to search the literature for)
- **PDE:** a stationary, axisymmetric **gauged Gross–Pitaevskii / nonlinear-Schrödinger** order parameter `ψ=ψ_r+iψ_i`
  coupled to a **Maxwell-like vector sector** `(a0, ar, aw)` (scalar potential + spatial A) plus a chemical potential
  `μ` and a throat-radius profile `r0(w)`. Quintic EOS: enthalpy `(5K/4)ρ⁴`, `ρ=|ψ|²` (so the matter nonlinearity is a
  C^∞ polynomial in (ψ_r,ψ_i) — **no √ρ / |ψ| / density division**; the residual is smooth on the positive-R0 domain).
- **Discretization:** finite-volume on a logically-rectangular **(r,w) grid, 16×16** for this work; state =
  `[ψ_r,ψ_i,a0,ar,aw]` per cell (5·256) + `r0[16]` + `μ[1]` = **1297 unknowns**. Axisymmetric vector divergence/curl
  use centered + one-sided-closure stencils (the discrete de Rham complex is NOT exact).
- **Maxwell block** = physical curl-curl `+ (1/ξ) grad(Z·div A)` soft gauge penalty (ξ=1, no boundary gauge term).
- **Solver:** Newton with Armijo backtracking line-search (`max_newton_iters=20`, `max_line_search_iters=20`),
  parameter **continuation in a depth knob τ** (default 8-step fine-halving sequence 0.03→0.028, `max_tau_backtracks=5`),
  warm-started each τ from the previous converged state. Jacobian assembled by **graph-colored JVP** (≈253 colors),
  solved sparse/dense. **Single-Arbiter principle:** the ORIGINAL physical residual is the sole convergence/merit
  arbiter; all conditioning aids (ε preconditioner schedule, Jacobi scaling) are path-only, final solve at zero ε.
- **Conditioning:** Jacobian cond ≈ 1e13. Near-null subspace (~5 modes): the **U(1) global-phase gauge mode** (exact
  null, σ≈7e-14), **Maxwell A-sector gradient/gauge modes** (modes 1/3/4, near-null), and **one transverse mode**
  (mode 2, ~17% non-gradient, σ≈2.6e-10).
- **Symptom:** continuation converges cleanly to τ=0.029125 then stalls; solution sensitivity `‖dψ/dτ‖` grows ~3.2×
  (R0-lane dominated) approaching the stall; reachable Linf asymptotes toward but does not cross the 1e-6 gate as the
  τ-step shrinks (turning-point-like creep).

**Search this as:** "pseudo-arclength continuation past folds in nonlinear elliptic BVP / Gross-Pitaevskii"; "Newton
globalization for near-singular Jacobian (Levenberg–Marquardt / trust-region / pseudo-transient continuation)";
"null-space deflation / gauge fixing in coupled EM–matter Newton solves"; "deflated continuation (Farrell)"; "computing
stationary GPE / vortex / dark-soliton ground states (Bao & Cai; Sobolev/preconditioned-gradient)"; "sparse Jacobian via
graph coloring (Gebremedhin, 'What color is your Jacobian?')"; "analog-gravity / acoustic-metric throat profile numerics".

---

## 2. The investigation chain (C0 → C0f2), compressed
- **C0** (2026-06-19): verdict `PRODUCTION_SOLVER_REQUIRED` — **REJECTED** by adversarial review (short-circuited:
  cold-loaded pre-existing solutions, broke the ε-loop on first fail, hardcoded admissibility, misdiagnosed).
- **C0b** (corrected crawl): `DIAGNOSTIC_INCOMPLETE` — established the wall is a **persistent near-null subspace** in the
  physics field block (NOT the mass/μ border, NOT a fold by its limited test). **BUT this run used the crippled config
  that caused the false "wall"** (2 Newton iters, 0 τ-backtracks, single τ-step) — see C0f.
- **C0c**: the worst near-null mode **IS the U(1) phase gauge mode** (confirmed: annihilation, overlap 1.0, equivariance).
- **C0d**: projected the 4 Maxwell-lane modes onto an A-only ∇λ basis → `MIXED` ("3/4 stiffness") — later judged
  OVERSTATED (hand-picked threshold).
- **C0e** (gauge-invariant curl gate): verdict `GAUGE_FRAMING_REFUTED` — **REJECTED** as a **k-biased-metric false
  negative** (the curl fraction `‖∂A‖/‖A‖` is dimensionally a wavenumber; it amplified the high-k remnant of
  99.9%-gradient modes). C0e-0 (Newton-step framing) showed the near-null subspace did NOT block the step **at τ=0.02899**.
- **C0f** (globalization probe, defaults): `DIAGNOSTIC_INCOMPLETE` — **the reframe**: the C0b "wall" was a crippled
  config; with defaults the crawl converges cleanly to τ=0.029125, then ONE step hit the 600s per-attempt timeout.
  Unbiased `1−P_G` re-confirm: modes 1/3/4 GAUGE, mode 2 the lone transverse candidate (settles the C0d/C0e dispute).
- **C0f2** (chunked, timed, user-authorized >600s): the **decisive timing + stall characterization** — §3.

---

## 3. VERIFIED findings from C0f2 (fidelity-audited faithful; numbers trustworthy)
Per-τ crawl (genuine warm-started continuation, NO cold-load — verified):

| τ | wall_s | iters | backtracks | Linf | status |
|---|---|---|---|---|---|
| 0.030000 | 170 | 0 | 0 | 3.4e-10 | converged |
| 0.029500 | 123 | 0 | 0 | 1.1e-12 | converged |
| 0.029250 | 126 | 0 | 0 | 6.6e-12 | converged |
| 0.029125 | 167 | 0 | 0 | 3.4e-11 | converged |
| 0.0290625 | 297 | 6 | 29 | 2.4e-5 | line-search exhausted |
| 0.02909375 | 425 | 9 | 56 | 1.0e-5 | line-search exhausted |
| 0.0291125 | 371 | 8 | 63 | 1.6e-6 | line-search exhausted |

- **Cost = 96–97% Jacobian assembly** (e.g. 286 s of 297 s); linear solve only 3–16 s. Assembly is a colored per-column
  JVP build (~1k–2.5k JVP probes/assembly). ⇒ sparse/analytic assembly is a ~100× target.
- **The stall is real and Newton is effectively stuck**, not "slow": at τ=0.0290625 the **merit sweep** (solve `Jδ=−F`,
  clean to rel-resid 1.05e-12; then evaluate ‖F(x+αδ)‖ over α):
  - α=1: predicted L2 = 9.2e-17 (linear model "predicts convergence") but ACTUAL L2 = 3.7e-4 (**4× worse** than the
    start) — a **13-order predicted-vs-actual gap** = the full step massively overshoots.
  - best α (1/128): L2 = 8.715e-5 vs initial 8.719e-5 = **0.05% reduction**. To reach the gate needs ~thousands of such
    steps. Honest label: **STALLED / GLOBALIZATION_INSUFFICIENT**, not "slow."
- **The amplified direction is NON-gauge.** A pure gauge step can't change the gauge-invariant ‖F‖; the full step
  *quadruples* it ⇒ the J⁻¹-amplified near-null direction has a **transverse component = mode 2**. ⇒ **gauge-only
  fixes (the C0d/C0e dream) are DEAD.**
- **Fold-onset is LIVE (~30%), not ruled out.** The fold detector's **backtrack-exhaustion condition fired**; only the
  growth-factor gate (3.2 < 10×, window cut off by timeout) held it at `FOLD_RISK`. The "smaller steps improve Linf"
  pattern is actually **turning-point creep** (the τ-backtracks move *toward* the converged 0.029125; reachable Linf
  cleans up only as |Δτ|→0, never crossing the gate). Confidence: ~30% genuine fold, ~55% curable ill-conditioning,
  ~15% mix.
- **Timing is MOOT until the stall is cured.** The ~1.67h-to-τ=0.02 extrapolation was built from the trivial
  (0-iteration) steps and assumes all remaining steps are trivial; the stall breaks that. Assembly speed is irrelevant
  to a step that won't converge.

---

## 4. THE OPEN QUESTION + the cheap test that settles it
**Is τ≈0.029 a genuine fold (turning point) or severe ill-conditioning?** Decisive, cheap (no new solver, no long run):
on the EXISTING converged states (0.03, 0.0295, 0.02925, 0.029125), compute **sign of det(J)** (or smallest-real-
eigenvalue sign), the **σ_min(τ) trend**, and the **branch tangent dx/dτ**:
- **Fold** ⇒ det(J) flips sign / `dτ/ds → 0` across the limit point ⇒ τ-continuation cannot pass ⇒ **pseudo-arclength**.
- **Ill-conditioning** ⇒ no sign flip, σ_min ~constant ⇒ curable by **better globalization** (LM / trust-region).

---

## 5. Candidate fixes (C0g) and status
| candidate | status | note |
|---|---|---|
| **Sparse / analytic Jacobian assembly** | **needed regardless** | ~100× cost win; assembly is 96–97% of each step |
| **Levenberg–Marquardt / trust-region** | **likely primary (if ill-conditioning)** | `+λI` lifts σ_min — regularizes the exact near-null/mode-2 direction that overshoots; bounds the step to the basin where the linear model holds. Same "medicine" as gauge insurance, done robustly. |
| **Pseudo-arclength continuation** | **pre-scoped, branch-if-fold** | the right tool IF the probe shows a fold; does NOT yet exist in the codebase |
| **Gauge / near-null deflation** | optional | subsumed by LM's λI; targeted rank-k only if uniform λ over-damps |
| **Production linear solver** | **OFF** | linear solve already clean (1e-12); this is a globalization, not linear-solve, problem |

---

## 6. Post-/compact plan (resume here)
1. **LITERATURE / ONLINE RESEARCH FIRST (user-requested).** Using §1's precise problem statement + search terms, find
   how this class of solve is done and what tricks apply BEFORE building C0g. Target method families: pseudo-arclength /
   deflated continuation (AUTO, pde2path, BifurcationKit.jl, Farrell deflation); pseudo-transient continuation (Kelley–
   Keyes); Levenberg–Marquardt / trust-region globalization; null-space deflation & gauge-fixing in coupled EM–matter
   BVPs; stationary-GPE/vortex/soliton solvers (Bao & Cai; Sobolev/preconditioned gradient; Newton-for-GPE); sparse
   Jacobian via graph coloring; analog-gravity throat-profile numerics; near-empty-core GPE conditioning.
2. **The fold probe** (§4) — cheap, decisive: fold vs ill-conditioning.
3. **C0g**, shaped by (1)+(2): sparse assembly + the chosen globalization (LM/trust-region) and/or pseudo-arclength;
   instrument det(J)-sign/tangent at each converged τ. Keep the Single-Arbiter principle + genuine continuation.
4. **If the crawl unblocks** to the physics-target depth ⇒ promote the constitutive family → calibrated branch →
   multi-knob calibrate-predict (R0/J/W → anchor → SURPLUS) → `pathA_22`.

## 7. Standing constraints (carry into any continuation)
Single-Arbiter (original `patha_closed_branch_residual` is the sole convergence/merit arbiter; aids path-only; final
solve zero-ε). GENUINE warm-started continuation — `prefer_existing_b2c_background_predictor=False`, NO cold-loading of
pre-existing solutions (the C0-v2 sin). Faithful operators / frozen physics / `physical_export_permitted` untouched;
depth continuation τ-only. The throat profile is a **calibration knob**, not derive-ab-initio (calibrate-predict).
**Hardware:** local GPU off; RunPod A100/H100 (FP64) is the confirmed high-res scale-up lever — but only AFTER the stall
is cured and only for production resolution (useless at 16×16). Reviews: clean fidelity + adversarial agents per run;
the user distrusts clean/surprising verdicts — verify, especially any "can't-do" (production-solver) or "all-fixed" call.

## 8. Artifacts (committed unless noted)
Directives: `directives/pathA_C0{b,c,d,e,f}_*.md` (+ C0f2 prompt in `_scratch/`). Code:
`src/stage1_solver/patha_c0_conditioning_spike.py` (+C0c/C0d/C0e), `patha_c0f_globalization_probe.py`,
`patha_c0f2_timing_rerun.py`. Reports: `reports/pathA_C0{b,c,d,e,f}_*.md`, `reports/pathA_C0f2_timing_rerun.md`, this
file. Runs/JSON under `runs/pathA_C0*`. Resume pointer: `decisions/13` §0 (points here).
