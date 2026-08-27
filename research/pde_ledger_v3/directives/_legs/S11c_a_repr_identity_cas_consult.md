# S11c-a T7 — from-spec CAS construction consult: are two surviving cross-engine residuals identities or findings?

You are one of two INDEPENDENT engines (the other builds the same objects separately; you will not see its
work). Two independently-built CAS engines computed the S11c-a **T7 interface shape-derivatives** — a SymPy
"PY" engine and a blind Wolfram "WL" engine. A reviewed comparator diffed all 39 tag families. After the
already-settled fixes, **two residual classes remain nonzero**, both a *difference in how the two engines
NAME/represent a quantity*. Your job: determine, **from the spec alone and in runnable CAS**, whether each
residual is **identically zero under the spec's own definitions** (a benign representational identity) or a
**genuine difference** (one engine has dropped a real dependence). Do not decide it in prose — build the
objects and compute.

## The spec (your ONLY source of truth — read it first, build from it)
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`. Relevant: §1a (coordinates, `ρ_br⁰≡ρ_4D⁰W₀`,
`e_W≡δW/W₀`), §1b (bulk current `j=ρ_4D v_bulk`, projection law `∫Ω(∂_tδρ+∇₄·δj)dw`), §1c (constitutive
derivatives are variational; `μ_θ ≡ (δU/δθ)|_{u,e_W,…fixed}`; `μ_s=μ_θ/ρ_br⁰`), §2a (background profiles
`W_bg(y)≡W̄₀[1+η w₁(ξ)]`, `μ_R,bg(y)≡μ̄_R[1+η m₁(ξ)]`, the **supplied first-derivative maps** `∂_{yᵢ}W_bg`,
`∂_{yᵢ}μ_R,bg`, and the truncation rule "first order in wave amplitude and first shape order in each
background bookkeeper"), §2b (the two density representatives RHO4-CONSTANT / RHOBR-CONSTANT), §2c (the two
anchorings LAB_HELD `Q_bg^L(x,t)≡Q_bg(x)` and MATERIAL_ADVECTED `Q_bg^M(x,t)≡Q_bg(χ(x,t))`), §3a (face maps;
`ρ_4D^α(x,t)≡ρ_4D,bg^{0,α}(x,t)(1+θ)`; `μ_s^α≡μ_θ^α/ρ_br,bg^{0,α}`), §3c (traced fields; no background
derivative may be introduced as a free premise).

## The two residuals (observable facts about the committed operands — NOT a verdict)
The committed run is `~/.s11_build/comparator_run.out` (per-case `operand_A` = PY, `operand_B` = WL,
`A_minus_B`; srepr). Extract the literal operands yourself if useful. The two representational differences:

**Q1 — the closure coefficient μ_θ (families TRACTION, CLOSURE_SHAPE_DERIV, VIRTUAL_WORK_SHAPE_DERIV).**
WL writes it as an APPLIED function `mu_theta_L(x1,x2,x3,time)` (and `mu_theta_M(…)`); PY writes it as a BARE
symbol `mu_theta_L` / `mu_theta_M`. Across those families the entire cross-engine residual is a sum of terms
each carrying the factor `( mu_theta_L − mu_theta_L(x1,x2,x3,time) )`.
- **Question:** does the spec make `μ_θ^α` a quantity with genuine in-plane (`x1,x2,x3`) and/or `time`
  dependence that the T7 shape derivative must act on — or a per-anchoring constant whose evaluation-point
  arguments are inert (so bare ≡ applied)? Ground your answer in what the spec **supplies** about `μ_θ`'s
  in-plane/time dependence (§1c, §2a, §3a) and §3c's free-premise rule. Then, in CAS, state the condition on
  `μ_θ` under which the residual is identically zero, and whether the spec supplies that condition.

**Q2 — the projected density perturbation (families PROJECTION_RESIDUAL / _SHAPE_DERIV / _STATIC_OPERAND /
_DYNAMIC_OPERAND / _TERM_ORIGINS).** In these, the in-plane current divergence `∂_iδj_i` appears in both
engines. What differs: PY carries the density-time term as a single symbol `delta_rho_4D_bulk_t`, while WL
carries it as an explicit expansion in `ρ_4D,bg^{0,α}` and `θ` (you will see `rho_br`, `W_0`/`W_bg`, `θ_t`,
and — in the MATERIAL_ADVECTED branch — material-velocity `u_i_t` terms).
- **Question:** per §3a `ρ_4D^α=ρ_4D,bg^{0,α}(1+θ)` (so `δρ_4D=ρ_4D,bg^{0,α}·θ` to first wave order), with the
  §2b representative and §2c anchoring, is PY's single symbol `delta_rho_4D_bulk` the SAME object as WL's
  explicit expansion? Build `δρ_4D` from the spec for **each of the four (branch × density) combinations**
  — LAB_HELD/RHO4_CONSTANT, LAB_HELD/RHOBR_CONSTANT, MATERIAL_ADVECTED/RHO4_CONSTANT,
  MATERIAL_ADVECTED/RHOBR_CONSTANT — including, for MATERIAL_ADVECTED, the time-derivative of the advected
  background `∂_t[ρ_4D,bg⁰(χ(x,t))]`, and using the §2a ansatz + first-shape-order truncation. Then determine
  whether the density-time residual vanishes under that spec-built `δρ_4D`, per combination. If it vanishes
  for all four → representational identity. If it does NOT vanish for some combination → a candidate finding
  (name which term/combination and why).

## Method — MANDATORY (a prose derivation is discarded)
- Read the spec, then **write a runnable SymPy script** that builds `μ_θ^α` (Q1) and `δρ_4D` per combination
  (Q2) FROM THE SPEC, and computes each residual condition. **Do not assume the identity** — compute whether
  it holds and report the **IFF condition** (e.g. "residual = 0 iff `μ_θ` is independent of `x,time`", "the
  density residual = 0 iff `δρ_4D = …`"). An honest conditional is the deliverable.
- You may read the committed operands from `~/.s11_build/comparator_run.out` to confirm the residual form,
  but **build `μ_θ`/`δρ` yourself from the spec** — do not import any name map, reconciler, or the comparator.
- Save BOTH your script AND its literal stdout to named absolute paths under `~/.s11_build/` (prefix
  `repr_identity_`), and report those paths. Without script + literal stdout your derivation is discarded.
- **Physics filter:** report a finding only if it catches a way the physics could be wrong (a real dropped
  dependence), not "PY would be wrong on a different input."
- Report per question: **IDENTITY** (residual ≡ 0 under a spec-supplied condition — name it) or **FINDING**
  (a genuine difference — name the surviving term and the combination), plus your IFF condition and the paths.
