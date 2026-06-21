# Path-A throat-solver — literature synthesis (5-lane online research, 2026-06-20)

**Purpose:** answer the user's post-/compact request — *"look online, find research on how this problem has been
computed, what tricks apply."* Five parallel research agents each owned a distinct literature; this is the synthesis +
the convergent recommendation. Companion to `pathA_throat_solver_findings_and_research_brief.md` (the problem statement).

The problem, in one line: a **natural-parameter continuation of a near-singular nonlinear elliptic BVP** (gauged
quintic Gross–Pitaevskii ψ coupled to a Maxwell-like A-sector + μ + throat profile r0(w)), which converges cleanly to
τ=0.029125 then stalls — full Newton step overshoots ~13 orders, cond(J)≈1e13, near-null subspace = gauge modes + one
transverse "mode 2." Open question: genuine fold/turning-point vs curable ill-conditioning.

---

## 0. The convergent picture (where ≥3 independent lanes agree)

1. **Gauge-fix the EXACT null modes first — it is necessary, cheap, and a precondition for everything else.**
   Three lanes reach this independently:
   - *Continuation (Lane 1):* a single pseudo-arclength border **cannot** desingularize a ≥2-dimensional null space.
     Govaerts' GMBE theorem: the bordered matrix is well-conditioned *only if the borders span the entire null space.*
     We have a confirmed U(1) phase null vector + suspected Maxwell gradient modes ⇒ one arclength row leaves Jb
     singular ⇒ **PALC fails the same way natural continuation did unless the gauge modes are removed first.**
   - *Maxwell (Lane 4):* the standard, robust, cheap cure is **discrete gauge fixing by pinning** — tree-cotree pin on
     A (one A-DOF per spanning-tree edge removes exactly `range(grad)`) + a **single reference-node phase pin on ψ**
     (pin Im ψ at a far-field node where |ψ|≫0, *not* the throat core). This removes ~4 of our ~5 near-null modes.
   - *GPE (Lane 3):* every standard GPE solver *assumes the gauge away* (works modulo phase, carries no vector
     potential). The gauge fix is **Abelian-Higgs / Ginzburg–Landau territory, not GPE** — we must supply it ourselves.

2. **The remaining transverse "mode 2" needs a REGULARIZED-STEP globalization, not a gauge fix.** Lanes 2 + 4 agree:
   - *Globalization (Lane 2):* **pseudo-transient continuation (Ψtc)** `(V/δ + J)δx = −F` with SER step control, or
     **Levenberg–Marquardt in shifted-Newton form** `(J + λI)δx = −F`, λ=μ‖F‖^δ. Both are a **one-diagonal change** to
     our existing assembled-J solve; both lift σ_min off 1e-11, cap the overshoot, and recover quadratic Newton as the
     iterate enters the basin. *Use the shifted form, NEVER `JᵀJ+λI`* (that squares cond 1e13→1e26, below float64).
   - *Maxwell (Lane 4):* a **trust-region / line-search cap** bounds the step regardless of the amplified direction;
     for a genuine 1-D non-gauge kernel, two-step singular-Newton (Jacobian + projected Hessian) restores quadratic rate.

3. **A NEW, sharper hypothesis for the open question: the throat may be a PHYSICAL sonic/horizon critical point
   (Lane 5).** The deep-throat limit of a draining superfluid is *physically* where the flow speed crosses the local
   sound speed (v→c_s). At that sonic point the **stationary** flow equations are genuinely singular (0/0 indeterminate,
   hyperbolic↔elliptic type change, saddle/turning-point structure) — a textbook obstacle in transonic-flow /
   Bondi-accretion / de-Laval-nozzle solvers, and the analog-gravity literature identifies *acoustic horizon = sonic
   point = critical point of the steady flow.* If true, τ≈0.029 is a **genuine fold at horizon formation**, not a solver
   bug, and naive Newton continuation *cannot* pass it by construction. **This reframes "fold vs ill-conditioning" into
   a testable physical prediction** (the v/c_s diagnostic, §2).

> **Genuine tension to resolve by measurement, NOT by preference:** Lane 5 says "expect a real sonic fold"
> (→ pseudo-arclength / shoot-from-critical-point). Lane 3 says "the literature strongly predicts your wall is
> stiffness + gauge, not a true fold" (→ inverse-Hamiltonian preconditioner + gauge-fix). These make *opposite*
> predictions and are separated by the cheap diagnostics in §2. Do not adopt either until the measurement is in.

---

## 1. Per-lane top recommendations + key references

### Lane 1 — Continuation past folds (pseudo-arclength, deflated continuation, software)
- **Verdict:** pseudo-arclength continuation (PALC, Keller) is the correct fix **iff** τ≈0.029 is a genuine *simple*
  fold — and **only after** the gauge null space is removed (else the single border can't desingularize a multi-dim
  null space). Hand-rolling Keller PALC on our colored-JVP Jacobian is ~50–100 lines and needs **no new linear-algebra
  primitives** beyond residual evals + forward JVP (only fold *tracking* needs transpose-JVP). Deflated continuation
  (Farrell) is *not* a fold-passer — it finds disconnected branches; useful only to check for a partner branch.
- **Refs:** Keller 1977 (PALC origin); Dickson–Kelley–Ipsen–Kevrekidis, *Condition Estimates for Pseudo-Arclength
  Continuation*, SIAM JNA 2007 (arXiv:math/0603716) — proves Jb well-conditioned where F_x is singular;
  Govaerts 1991, *Stable Solvers and Block Elimination for Bordered Systems*, SIMAX (GMBE theorem); Moore–Spence
  turning-point system; **BifurcationKit.jl** (matrix-free/JVP, fold continuation + deflation — closest off-the-shelf
  fit), LOCA/Trilinos (Householder bordering, large-scale), AUTO-07p, pde2path. Farrell–Birkisson–Funke, *Deflation
  Techniques…*, SIAM JSC 2015.
- **Decisive measurement it contributes:** the **bordered σ_min(Jb)** test — border with the fold tangent; if
  σ_min(Jb)=O(1) while σ_min(F_x)≈1e-11 ⇒ fold ⇒ PALC passes. If Jb is *still* near-singular ⇒ it's gauge/conditioning,
  not a simple fold.

### Lane 2 — Newton globalization for near-singular Jacobians
- **Verdict (primary candidate overall):** **pseudo-transient continuation (Ψtc)** — near-drop-in over our existing
  τ-continuation+warm-start; `(V/δ+J)δx=−F` + SER (`δ_{n+1}=δ_n·‖F_n‖/‖F_{n+1}‖`, reject+halve δ if residual grows).
  Doubles as a **fold detector**: δ→∞ & ‖F‖→0 ⇒ it was conditioning (wall broken); δ stalls & ‖F‖ floors >0 ⇒ genuine
  fold ⇒ switch to PALC. **Levenberg–Marquardt** (shifted-Newton form, λ=μ‖F‖^δ) is the same machine with a
  self-annealing λ and a *proven* singular-J convergence theory (Yamashita–Fukushima; Fan–Yuan). Pick the regularizer
  in **[1,100]** (≫σ_min=1e-11, ≤σ_max=826) ⇒ regularized cond ≈ σ_max/λ ≈ 10²–10³.
- **Cheap off-the-shelf sanity check (do BEFORE coding):** run `scipy.optimize.root(F, x_stall, method='lm')` and
  `method='hybr'` from the stalled state. n=1297 runs fine as a diagnostic. If either makes real progress where
  Newton+Armijo stalls ⇒ regularized-step methods are the answer AND it's not a hard fold.
- **Refs:** Kelley & Keyes, *Convergence Analysis of Pseudo-Transient Continuation*, SIAM JNA 1998; Coffey–Kelley–Keyes
  2003 (PTC for DAEs — relevant: r0[16]+μ make our system DAE-like); Yamashita–Fukushima 2001; Fan–Yuan 2005; Fischer
  et al. 2024 (modern LM survey); Nocedal & Wright Ch.4 (LM↔TR↔regularized-Newton equivalence); Steihaug 1983
  (truncated-CG TR that survives singular J); PETSc SNES (NEWTONTR), KINSOL, scipy root.
- **Honest caveat:** if it's a *true saddle-node fold*, **none of LM/PTC/TR passes it** — PTC stagnates (δ stops
  growing), LM/TR "converge" to a nonzero-residual non-root (the residual flooring above zero is the only signal).
  Then the tool is PALC. Also: regularization papers over the *exact* gauge nulls with arbitrary gauge drift — cleaner
  to gauge-fix first, then let λ handle only the lone transverse mode.

### Lane 3 — Stationary GPE / vortex / deep-core solvers
- **Verdict:** **Sobolev / inverse-Hamiltonian preconditioning** is the single highest-leverage, lowest-risk fix for
  the *matter block* — it targets exactly "Newton overshoots as the core empties." Key result (Henning–Peterseim 2020):
  the right metric is `(linearized Hamiltonian)⁻¹`; in that metric the gradient = identity ⇒ the preconditioned step is
  order-1 conditioned, mesh-independent, and **robust to strong nonlinearity** (our deepening quintic core). Apply as a
  **left preconditioner to the colored-JVP linear solve** — one solve with the matter-block elliptic+potential part we
  already assemble. Secondary: **backward-Euler imaginary-time** (unconditionally energy-diminishing, no δt cap from the
  stiff core) as a robust per-τ *warm-start*; the **damped J-method** (Altmann–Henning–Peterseim; Jarlebring et al.) as
  a μ-aware globalized replacement for Newton+Armijo with a spectral shift to *select the depleted branch*.
- **Refs:** Bao & Cai, *Mathematical theory and numerical methods for BEC*, KRM 2013 (canonical review); Henning &
  Peterseim, *Sobolev Gradient Flow for the GPE Eigenvalue Problem*, SIAM JNA 2020 (arXiv:1812.00835); Danaila &
  Protas, *Computation of Ground States… via Riemannian Optimization*, SIAM JSC 2017; Altmann–Henning–Peterseim,
  *The J-method for the GPE eigenvalue problem*, Numer. Math. 2021; Danaila & Hecht 2010 (FE toolbox, vortices);
  Antoine & Duboscq 2017 (accelerated imaginary time).
- **Decisive measurement it contributes:** re-run the near-null SVD **in the preconditioned (a_u) metric** + a
  **Bogoliubov–de Gennes (BdG) zero-mode / Krein-signature** check. Pure ill-conditioning ⇒ the small singular-value
  *cluster collapses* in the energy metric. A genuine fold/bifurcation ⇒ a *single isolated* eigenvalue crosses zero
  smoothly and survives the metric change. **This is the GPE community's own fold-vs-conditioning discriminator.**
- **Honest caveat:** GPE methods do **not** touch the Maxwell/gauge sector (must be gauge-fixed separately); the
  "metric = identity" optimality is proven for the cubic nonlinearity, qualitative for our quintic; plain imaginary
  time relaxes to the *ground* state and will fight a depleted/excited branch (need the shifted J-method).

### Lane 4 — Maxwell curl-curl gauge / null-space handling
- **Verdict:** **(a)** pin ψ global phase at a far-field reference node (one line — strike the DOF); **(b)** tree-cotree
  pin on A (spanning tree of the grid graph, pin one A-DOF per tree edge — removes exactly the discrete-gradient
  kernel); **(c)** null-space deflation of the Newton update `δ ← δ − V(VᵀV)⁻¹(Vᵀδ)` as a cheap safety net for the
  residual gradient modes the non-exact stencil leaves slightly off `range(G)`. **Our soft `(1/ξ)grad(div A)` penalty
  is part of the problem, not the cure** — penalties provably leave residual spurious modes (too small) or damage
  eigenvalues/conditioning (too large); a *hard* pin is the literature-endorsed replacement.
- **The non-exact-complex finding:** our centered + one-sided-closure stencils break `curl_h∘grad_h = 0`, so
  gradient directions become *near*-null (rotated slightly out of `range(G)`) — this manufactures spurious near-null
  modes and *plausibly leaks the "17% non-gradient" into mode 2.* Cheap tests: **(1)** commutator norm `‖curl_h·G‖`
  (≈1e-15 if exact; O(h) localized to boundary rows if polluted); **(2)** overlap `‖P_G v‖/‖v‖` per mode (re-confirm
  1/3/4≈gauge, mode 2≈0.83); **(3)** refinement scaling of the near-null gap (O(h)→pollution; O(1)→real).
- **Refs:** Hiptmair–Xu 2007 (auxiliary-space/AMS — but *secondary* for us since our linear solve is already clean);
  Arnold–Falk–Winther, *Finite Element Exterior Calculus*, Acta Numerica 2006; Hyman–Shashkov (mimetic FD on logically
  rectangular grids); tree-cotree gauge (arXiv:2502.13588); Costabel–Dauge 2002 (weighted regularization — why the
  penalty is incomplete); Du–Gunzburger–Peterson 1992 (Ginzburg–Landau gauge choices); lattice Abelian-Higgs gauge
  fixing (arXiv:2305.15236); Pestana–Rees 2016 (null-space preconditioners); two-step Newton for deflation-one
  singular zeros (arXiv:2305.10803 — for a *genuine* 1-D kernel = mode 2).
- **Honest caveat:** gauge fixing removes **4 of 5** modes; **mode 2 survives** and still drives the overshoot ⇒ gauge
  fix is *necessary but not sufficient* — pair with a trust-region/LM/PTC cap (and possibly a bordered solve if mode 2
  is a near-bifurcation).

### Lane 5 — Sparse Jacobian assembly + analog-gravity throat numerics
- **Verdict A (assembly, 96–97% of each step):** **(1)** our **253 colors is likely 5–10× wasteful** — colors scale
  with the *stencil footprint*, not grid size (a 5-field compact stencil should color to ~30–60). **BUT check first**:
  a high count is *correct* if μ is a global Lagrange-multiplier row coupling to all unknowns, or r0(w) adds long-range
  coupling — a single dense row forces many colors. **(2)** the real win is **analytic/symbolic sparse stencil
  assembly** (residual is a smooth local polynomial ⇒ exact per-block partials, one vectorized pass, no probing ⇒ the
  ~100× target). **(3)** **JFNK is the WRONG lever here** — our linear solve is already clean+cheap, so matrix-free
  buys nothing and still needs an assembled preconditioner near the throat.
- **Verdict B (the big one):** **test whether the throat is a sonic/horizon critical point** — if `v/c_s → 1` at the
  stall, the wall is a *physical* fold (sonic point = acoustic horizon = saddle critical point of the steady flow),
  curable only by **pseudo-arclength** (cheapest) or **shoot-from-the-critical-point + L'Hôpital** (most faithful), NOT
  by more Newton iters or a faster Jacobian. C0f ("wall = crippled-defaults artifact") and "a real sonic fold deeper
  in" can BOTH be true: defaults fixed the τ=0.029 number; pushing further should hit the real turning point.
- **Refs:** Gebremedhin–Manne–Pothen, *What Color Is Your Jacobian?*, SIAM Review 2005; Coleman–Moré 1983; Knoll &
  Keyes, *JFNK survey*, JCP 2004; SparseDiffTools.jl; torch.func jacfwd+vmap. Analog gravity: Barceló–Liberati–Visser,
  *Analogue Gravity*, Living Rev. Rel. 2005/2011; *Hydraulic jump as a white hole* (arXiv:physics/0508215 — steady eqns
  degenerate / turning-point at the sonic transition); Das/Bilić, accretion-as-analogue-gravity (arXiv:0704.3618,
  1205.5506 — critical point = sonic point = acoustic horizon, saddle, L'Hôpital); transonic critical-point analysis
  (SIADS).
- **Honest caveat:** assembly speedup **does NOT cure the stall** (it makes a non-converging step cheaper, nothing
  more); the sonic-point picture is an *emergent leading-order* GPE description — quantum pressure can *regularize*
  (smear) the singularity, and the gauge sector can shift/split the critical point ⇒ **a testable prediction, not an
  identity.**

---

## 2. The cheap decisive diagnostic battery (run on EXISTING converged states — no new solver, no long run)

The research converges on a small set of measurements that *together* settle fold vs ill-conditioning vs sonic-point,
and tell us which C0g machinery to build. Cheapest/most-decisive first:

1. **v/c_s at the stall (Lane 5).** Compute radial flow speed vs local sound speed from the GPE/Madelung fields at the
   converged states approaching τ=0.029. `v/c_s → 1` ⇒ sonic horizon ⇒ genuine physical fold. *Single cheapest test;
   it can pre-empt everything.*
2. **Near-null vector localization (Lanes 3+5).** Where does the worst surviving (non-gauge) mode live? On the
   flow/density fields at the sonic radius (physics fold) or on the gauge DOFs (gauge)? We already know 1/3/4 are gauge;
   pin down where **mode 2** localizes.
3. **scipy LM/hybr from the stalled state (Lane 2).** Off-the-shelf, no torch code. Progress ⇒ regularized-step is the
   fix and it's not a hard fold. No progress ⇒ supports fold.
4. **det(J)-sign / σ_min(τ) / tangent τ̇-sign across the existing converged τ (Lane 1).** The defining fold test:
   τ̇→0 with a sign change ⇒ fold. (This is the original brief §4 probe, now sharpened.)
5. **Bordered σ_min(Jb) with the fold tangent (Lane 1).** σ_min(Jb)=O(1) while σ_min(F_x)≈1e-11 ⇒ fold PALC can pass.
6. **Commutator norm ‖curl_h·G‖ + per-mode overlap ‖P_G v‖/‖v‖ (Lane 4).** Quantifies discrete-complex non-exactness
   and whether part of mode 2 is a stencil artifact.

(Optional confirmatory: re-run the near-null SVD in the inverse-Hamiltonian metric + a BdG zero-mode/Krein check — the
GPE community's own discriminator; a cluster that collapses in the energy metric = conditioning, an isolated surviving
crossing = fold.)

---

## 3. The C0g build, shaped by the research (decision tree, gated on §2)

**Unconditional, do regardless (necessary precondition + cost):**
- **Gauge fix**: pin ψ global phase (far-field node) + tree-cotree pin on A; demote the soft grad-div penalty to a
  complement; deflate residual gradient modes from the Newton update. (Removes 4 of 5 near-null modes.)
- **Assembly**: first *diagnose* why 253 colors (check for the dense μ/r0(w) rows), then move to analytic sparse stencil
  assembly. NOT JFNK. (Cheaper steps; does not cure the stall.)

**Then branch on the §2 verdict:**
- **If conditioning / stiffness (no fold, v/c_s<1, σ_min(τ) flat, scipy-LM progresses):** add **pseudo-transient
  continuation** (primary) and/or **Levenberg–Marquardt** shifted-Newton (`J+λI`, λ∈[1,100] annealing) + the
  **inverse-Hamiltonian (Sobolev) preconditioner** on the matter block. This is the most likely outcome per Lanes 2+3.
- **If genuine fold / sonic-horizon (v/c_s→1, τ̇ sign change, det(J) flip, bordered σ_min(Jb)=O(1)):** build
  **pseudo-arclength continuation** (Keller bordered system on the colored-JVP Jacobian, ~50–100 lines) — *only after*
  the gauge fix, so the single border desingularizes the lone fold null vector. Most-faithful alternative if the throat
  is truly transonic: **shoot from the sonic point with L'Hôpital regularization.**

**Keep throughout:** Single-Arbiter principle (original `patha_closed_branch_residual` is the sole convergence/merit
arbiter; aids path-only; final solve at zero ε); GENUINE warm-started continuation (`prefer_existing…=False`, no
cold-load); faithful operators / frozen physics / physical-export guard untouched. Hardware: local CPU sparse-direct
for low-res; RunPod A100/H100 (FP64) only for the post-cure production resolution.

---

## 4. Recommended next action

1. **This synthesis → Codex consult** (Claude+Codex agree the diagnostic-battery directive + the C0g decision tree),
   then **→ GLM 5.2 review** (the established third-eyes loop), fold in.
2. **Run the §2 diagnostic battery** (cheap, decisive — settles fold vs conditioning vs sonic-point) as the next
   executable directive (C0g-diag), design-reviewed, then the user execution gate.
3. **C0g build** per §3, branched on the battery's verdict.

Resume pointer: `decisions/13` §0 → this file + the §2 battery as the next concrete step.

---

## 5. Codex consult — agreed corrections + final design (Claude + Codex, 2026-06-20)

Codex read the actual solver code and corrected several overclaims in §0–§4. **These corrections supersede the
prose above where they conflict.** This is the design we take to GLM 5.2.

### 5.1 Code-grounded corrections (what changed and why)
1. **Gauge handling is PATH-ONLY — never a physics change.** §0/§3 said "demote the soft grad-div penalty." That edits
   the *frozen residual/operators* → forbidden. The gauge pins (ψ far-field phase pin, A-sector gauge handling) live in
   the **solver coordinates / preconditioner**, not in `patha_closed_branch_residual` or `operators.py`. The original
   residual stays the sole arbiter; final solve at zero ε. (State packing confirmed:
   `[psi_real,psi_imag,a0,ar,aw,r0,mu]`, `coupled_branch.py:137`; arbiter `:512`; ε aids preconditioner-only
   `patha_c0_conditioning_spike.py:318`.)
2. **scipy `root(lm/hybr)` is a CAPPED SOLVER PROBE, not a cheap readout.** It builds a dense FD Jacobian at n=1297
   (minutes) and uses a *different* globalization + least-squares merit. Cap `maxfev`/wall-time; **judge only by the
   original residual.** Gate: progress = `Linf≤1e-6` OR original-residual `L2` drops ≥10×; anything less = "no
   evidence," NOT "fold proven."
3. **No det(J) sign-flip is observable from our states** — they are all on ONE side of the putative fold, so a sign
   change literally cannot appear. **Do not route on raw det.** Route on the *scaled / gauge-projected* `σ_min/σ_max`
   trend, the tangent `|τ̇|=|dτ/ds|→0`, and `‖x_τ‖` blow-up. (slogdet of a scaled, gauge-pinned J is confirmatory only.)
4. **Bordered `σ_min(Jb)` is meaningful only AFTER gauge projection** — raw Jb stays small while the U(1)/Maxwell gauge
   modes remain, giving a false "not a fold." Compare *scaled/projected* `cond(J)` vs `cond(Jb)`: **fold-support if
   `cond(J)>1e10` and `cond(Jb)<1e8`, or `σ_min(Jb)/σ_min(J) > 1e4`.** Absolute O(1) thresholds are scale-dependent —
   don't use them.
5. **Mode-2 SVD must be re-run at the CURRENT C0f2 converged states** — the C0d/C0e modes are stale (from C0b
   τ=0.02899/0.03). Track modes by vector overlap (indices swap). Avoid a hand-picked "throat window" (can't-fail-gate).
6. **253 colors is NOT a bug** — it's the deterministic radius-3, per-field coloring: `5·7² + 7 + 1 = 253`
   (5 cell fields × 49 + r0 wall 7 + μ 1; mass row inserted analytically). The lit-agent's "~30–60" assumed a
   distance-2/radius-2 coloring; the code uses `color_separation > 2·stencil_radius` (radius 3). The assembly win is
   **analytic/sparse assembly** (still applicable) + possibly a smaller coloring radius — and it **must** include the
   dense/low-rank rows for μ, the mass constraint, and the wall radial reductions (it is NOT a pure local 2D stencil).
   Verify in `_closed_coupled_color_groups`, `_closed_cell_column_rows`, `_closed_wall_column_rows`,
   `_constraint_row_entries`. **Reminder (§5 honesty): assembly speed does not cure the stall.**

### 5.2 The gauge-covariant Mach (sonic) diagnostic — agreed formula
Use the code's hydrodynamic variables (gauge-covariant current, NOT phase-gradient alone — `coupled_branch.py:330`):
```yaml
rho      : psi_real^2 + psi_imag^2
j_r      : (hbar/m)*(psi_real*d_r psi_imag - psi_imag*d_r psi_real) - (q/m)*A_r*rho
j_w      : (hbar/m)*(psi_real*d_w psi_imag - psi_imag*d_w psi_real) - (q/m)*A_w*rho
v_r, v_w : j_r/max(rho,floor), j_w/max(rho,floor)
speed    : sqrt(v_r^2 + v_w^2)
c_s      : sqrt(5*K/m) * rho^2          # from c_s^2 = rho*dh/drho/m, h=(5K/4)rho^4
Mach     : speed / c_s
```
Operational "sonic at a radius": full 2D Mach map on cell centers; report argmax under a density mask + `(r*,w*,r*/R0(w*))`
(R0 enters the confinement potential, not a hard boundary — do NOT assume the sonic point sits on R0). **Sonic-support
only if** `0.8 ≤ M_max ≤ 1.2`, **stable under a density-floor sweep**, and **M_max monotonically → 1 as τ decreases**.
Caveats in the directive: quantum pressure smears the hydrodynamic singularity; the Maxwell sector can shift/split the
critical mode; near-empty ρ makes `c_s→0` and `v=j/ρ` unstable; **if the current is ~0, the sonic hypothesis is not
represented in these stationary states** (report `j_at_max`). **`M≪1` rules out *sonic* but NOT a generic fold.**

### 5.3 Jacobian-spectral implementability (confirmed against code)
- `J·v` (colored JVP) AND an **explicit sparse `J`** (`assemble_closed_coupled_colored_sparse_jacobian`) are both
  available ⇒ **transpose-JVP NOT needed** (we have explicit J for any left-null/SVD work). `F_τ` = centered FD in τ at
  fixed x using the same residual + s_σ provider.
- Tangent: solve `J·x_τ = −F_τ` reusing the assembled/factored J; normalize `(x_τ,1)`; for decreasing-τ orientation
  `τ̇ = −1/√(1+‖x_τ‖_W²)`. `τ̇→0` is observable; a sign change is not (one-sided states).

### 5.4 The decision tree (refined, machine-gated)
```yaml
unconditional:                      # do regardless; preserve fidelity
  preserve: [original residual = sole arbiter, prefer_existing=False,
             no operator/physics/export-guard change during diagnostics]
  build:    [far-field phase pin (large rho), A-sector gauge handling in solver/preconditioner coords,
             analytic/sparse assembly AFTER the color audit]
  refrain:  [do NOT remove/demote grad-div from the residual here; any demotion is path-only or a separate decision]
conditioning_branch:                # most likely per Lanes 2+3
  if:   [max_Mach < 0.7 OR no stable approach to 1; projected sigma_min does NOT collapse >10x across tau;
         scipy bounded LM/hybr gets >=10x residual drop or hits original Linf gate; bordered Jb not better than J]
  build:[shifted-Newton/LM (J+lambda*I)dx=-F, lambda in [1,100] annealing; or PTC+SER;
         inverse-Hamiltonian/Sobolev preconditioner on the matter block]
fold_or_sonic_branch:
  if:   [Mach -> 1 at a stable non-vacuum location; mode-2 energy overlaps that location;
         projected sigma_min collapses AND |tau_s| -> 0; bordered projected cond improves >= 1e4]
  build:[gauge-fixed pseudo-arclength continuation; optional later: shoot-from-sonic + L'Hopital]
mixed_branch:
  if:   [gauge/stencil modes dominate but mode 2 stays transverse; Mach inconclusive]
  build:[gauge-fix + LM/PTC first, THEN rerun the spectral/Mach battery before committing to PALC]
```

### 5.5 Codex's standing disagreements (kept as live caveats)
- "v/c_s pre-empts everything" — **overstated**; confirms sonic if positive, cannot rule out a *non-sonic* fold.
- scipy root — capped probe, original-residual judged (5.1#2).
- raw det sign — too fragile pre-gauge-fix (5.1#3).
- bordered σ_min — needs scaling/projection (5.1#4).
- "253 wasteful" — false; deterministic radius-3 coloring (5.1#6).
- "demote grad-div" — conflicts with frozen operators; path-side only (5.1#1).
- "PALC ≈ 50–100 lines" — optimistic once τ-FD, gauge fixing, scaling, tangent normalization, and output gates are in.

### 5.6 Open questions for GLM 5.2 (ANSWERED in §6 — kept for the record)
1. **The one-sided-states limitation:** we can see `σ_min(τ)` trend, `τ̇→0`, `‖x_τ‖` blow-up — but NOT a det/τ̇ sign
   flip (no states past the fold). Is there a cheaper way to get a *both-sides* signal short of actually stepping
   through (e.g. a small deflated/perturbed solve, or a local quadratic-fold fit to σ_min(τ))? Or is the σ_min-collapse
   + bordered-Jb-improvement combination already sufficient to call a fold without a sign flip?
2. **Sonic indicator for a *gauged* GPE:** is the gauge-covariant Mach (5.2) the right critical-point indicator, or
   should the discriminant be the *characteristic/type-change* of the linearized stationary operator (eigenvalues of
   the principal symbol coalescing) rather than a hydrodynamic Mach ratio — given quantum pressure regularizes the
   hydrodynamic limit?
3. **Any battery member that still won't discriminate**, or a cheaper/more-decisive test we both missed?

---

## 6. LOCKED diagnostic-battery spec (Claude + Codex + GLM 5.2, 2026-06-20)

GLM's review (accepted in full — the math was checked: the wᵀF_τ transversality condition, the √|τ−τ_fold| ⇒
σ_min² linear scaling, and the always-elliptic principal-symbol argument are all standard and correct) added a
premise-validating **Step 0**, a cheaper decisive fold test (**wᵀF_τ**), a free fold-vs-bifurcation fit
(**σ_min²(τ)**), and closed the sonic-indicator question. **This §6 is the authoritative spec for the C0g-diag
directive; where it conflicts with §2/§5, §6 wins.**

### 6.1 What GLM changed (the deltas)
- **NEW Step 0 — framing check at the ACTUAL stall τ=0.0290625 (highest-value catch).** The "mode 2 drives the stall"
  premise rests on C0e-0's near-null fraction = 1.8e-21 — but that was measured at the *old crippled-solver* stall
  τ=0.02899, a DIFFERENT state. As τ↓ and J becomes more singular, the near-null direction can rotate to align with the
  residual (fraction grows — the fold signature). **If the near-null fraction is still ~0 at τ=0.0290625, the entire
  mode-2 / gauge-fix / fold analysis must be rethought.** Cost: 1 solve + a few residual evals on the saved stalled
  state. From the C0f2 merit sweep GLM estimated the broken quadratic model: ‖H(δ,δ)‖≈7.47e-4 ≈ 8.6×‖F‖ — the step ‖δ‖
  is large enough to break the model; Step 0 decides WHY (mode 2 in the step → LM is the fix; mode 2 absent → revise).
- **wᵀF_τ is the PRIMARY fold test; bordered cond(Jb) is demoted to confirmatory.** A simple fold ⟺ the bordered J_b is
  nonsingular ⟺ (i) gauge-projected J has a 1-D null space, (ii) wᵀv≠0 (automatic), (iii) **wᵀF_τ ≠ 0** (the fold
  regularity / transversality condition; a bifurcation has wᵀF_τ=0). w = left-null vector = u_min from the
  gauge-projected SVD (already computed). F_τ = centered FD (2 residual evals). Decisive yes/no from one-sided data;
  cost = 1 SVD + 2 residual evals + 1 dot product. *Must be gauge-projected first* (else w is a gauge mode).
- **σ_min²(τ) linear fit — free fold-vs-bifurcation + τ_fold estimate.** Simple fold: σ_min ∝ √|τ−τ_fold| ⇒ σ_min²
  linear in τ (zero-crossing ⇒ τ_fold). Bifurcation/conditioning: σ_min ∝ |τ−τ| ⇒ σ_min² quadratic/flat. Fit the 4
  gauge-projected σ_min(τ). Requires a clearly MONOTONE σ_min trend to be trustworthy (only 4 points, σ_min~1e-10–1e-11).
- **Sonic question CLOSED: the principal-symbol type-change test is a DEAD END.** In ψ-form the stationary operator is
  *always elliptic* — principal symbol = the Schrödinger Laplacian ħ²|ξ|²/2m (>0); gauge terms (a·∇ψ, |a|²ψ) are
  lower-order; Maxwell curl-curl+grad-div is positive-definite. No eigenvalue coalescence / type change to detect. The
  hydrodynamic elliptic→hyperbolic change at M=1 lives only in the Madelung-reduced Euler frame and is *smeared* (not
  eliminated) by quantum pressure. ⇒ **Mach = physical CONTEXT (why a fold might exist); the Jacobian test (wᵀF_τ /
  σ_min²) = the DISCRIMINANT (whether one exists).** Complementary, not alternatives.
- **scipy LM MUST carry a branch-overlap gate** (the cheapest both-sides signal). LM `(J+λI)δ=−F` can step THROUGH a
  fold. After LM converges at τ=0.0290625, overlap the result with the pre-stall solution: **>0.99 = same branch**
  (stall was conditioning, or the fold is deeper); **<0.5 = post-fold branch found** (fold confirmed); **LM fails** =
  supports fold (not definitive). Overlap = 1 dot product; it is a required gate, not optional.
- **Best-α descent trend across τ** (GLM addition B; nearly free, code exists): run the merit sweep at the converged
  states too. Best-α reduction degrading monotonically toward the stall ⇒ fold (direction → tangent); roughly constant
  ⇒ fixed conditioning. (Bonus already-diagnostic: the C0f2 predicted-vs-actual gap scales ~α^1.4, between O(α) and
  O(α²) — confirms near-null 1/σ_min amplification, not a stale Jacobian or nonsmoothness; residual is C^∞:
  (5K/4)|ψ|⁸ψ.)

### 6.2 The locked battery — ordering + per-step machine gates
| # | Test | Cost | Output / gate |
|---|---|---|---|
| **0** | **Framing check at the ACTUAL stall τ=0.0290625** (full-Newton step decomposed onto gauge vs near-null vs residual) + best-α merit sweep here | 1 solve + few residual evals | `near_null_component_fraction`, `gauge_component_fraction`, `transverse(mode2)_fraction`, best-α & %reduction. **PREMISE GATE:** if near-null fraction ~0 ⇒ HALT, rethink (stall ≠ mode-2 amplification). |
| **1** | **Gauge-projected SVD** at the 4 converged states (τ=0.03, 0.0295, 0.02925, 0.029125) | 4 SVDs (explicit sparse J available) | `σ_min(τ)`, `σ_min/σ_max`, left-null `w(τ)`, right-null `v(τ)`, mode tracking by overlap (indices swap). |
| **2** | **wᵀF_τ** (primary fold test) | 2 residual evals + 1 dot product (reuses w from #1) | `wᵀF_τ` per state. **\|wᵀF_τ\| clearly ≠ 0 ⇒ fold-support; ≈ 0 ⇒ bifurcation.** Gauge-projected. |
| **3** | **σ_min²(τ) linear fit** | free (from #1) | linear (nonzero slope) ⇒ fold + `τ_fold` estimate; quadratic/flat ⇒ bifurcation/conditioning. Trust only if σ_min monotone. |
| **4** | **Gauge-covariant Mach map** (§5.2) | cheap | `M_full_max`, `M_w_max`, `rho_at_max`, `j_at_max`, `(r*,w*,r*/R0)`, density-floor-sweep stability, monotone-→1 trend. Context-only; report `j_at_max` (≈0 ⇒ sonic not represented). |
| **5** | **scipy LM/hybr capped probe + branch-overlap** | minutes (capped `maxfev`/wall) | judged by ORIGINAL residual: progress = `Linf≤1e-6` or `L2↓≥10×`. **overlap gate:** >0.99 same branch; <0.5 post-fold branch (fold). |
| **6** | **Bordered gauge-projected cond(J) vs cond(Jb)** (confirmatory) | 1 factorization | fold-support if `cond(J)>1e10 & cond(Jb)<1e8` or `σ_min(Jb)/σ_min(J)>1e4`. Confirms #2. |
| **7** | **Commutator ‖curl_h·G‖ + per-mode P_G overlap** | cheap (op) + SVD modes | `‖C·G‖/‖G‖`, boundary-row fraction, per-mode `P_G`, residual curl. Mode-2 characterization (stencil-artifact vs real), NOT fold discrimination. |

### 6.3 Routing (the 2×2 + the fold tests)
- **Premise (Step 0) fails** (near-null fraction ~0 at the real stall) ⇒ **HALT and rethink** — do not proceed to the branch build.
- **Fold confirmed** = wᵀF_τ≠0 (gauge-projected) **AND** σ_min² linear **AND** (LM overlap<0.5 **OR** cond(Jb)≪cond(J)):
  - **+ Mach≈1 at a stable non-vacuum location, mode-2 energy overlaps it** ⇒ *sonic fold* ⇒ gauge-fixed
    pseudo-arclength (cheapest) or shoot-from-sonic + L'Hôpital (most faithful).
  - **+ Mach≪1** ⇒ *non-sonic fold* ⇒ gauge-fixed pseudo-arclength.
- **Conditioning** = wᵀF_τ≈0 **OR** σ_min² quadratic/flat, LM gets ≥10× drop / same branch ⇒ shifted-Newton LM
  (`J+λI`, λ∈[1,100] anneal) and/or PTC+SER + inverse-Hamiltonian (Sobolev) preconditioner on the matter block.
- **Mixed** (gauge/stencil modes dominate but mode 2 stays transverse; Mach inconclusive) ⇒ gauge-fix + LM/PTC first,
  THEN rerun the spectral/Mach battery before committing to PALC.
- **Unconditional regardless:** gauge-fix **path-only** (far-field ψ phase pin + A-sector gauge handling in
  solver/preconditioner coords — NEVER the frozen residual/operators); analytic/sparse assembly AFTER the color audit
  (handle the dense μ/mass/wall rows). Preserve Single-Arbiter + genuine warm-started continuation.

### 6.4 Next action
Cut the formal **C0g-diag directive** from §6 → Codex design-review (the directive, per the design-review standing
step) → user execution gate. Steps 0–3 are nearly free and premise-deciding; run them first and REPORT before the
costlier 4–7, so a premise failure (Step 0) or an early decisive verdict (Step 2/3) can short-circuit the rest.

