# S11c-a T7 — from-spec ADJUDICATION: μ_θ representation and the projected density perturbation

You are one of two INDEPENDENT adjudicators (the other works separately; you will not see its work, and there
is a separate CAS consult you also will not see). Two independently-built CAS engines computed the S11c-a
**T7 interface shape-derivatives** — a SymPy "PY" engine and a blind Wolfram "WL" engine. Their cross-engine
comparison leaves **two nonzero residual classes**, each a *difference in how a quantity is represented*. Your
job is to **read the spec first, form your own view, then rule** — for each class — whether the spec makes the
two representations the **same physical object** (a benign notational identity, no finding) or whether one
engine has **dropped a real dependence the spec requires** (a finding, and which engine). This is a physics
decision grounded in the spec, not a code-style judgement.

## Read FIRST (form your own derivation before looking at the engine forms below)
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — especially:
- §1c: `μ_θ ≡ (δU/δθ)|_{u,e_W,and all other fields fixed}` (a variational derivative; the held-fixed qualifier
  is binding), and `μ_s=μ_θ/ρ_br⁰`.
- §2a: the background varies in-plane — `W_bg(y)`, `μ_R,bg(y)` — and the spec **supplies first-derivative maps
  ONLY for those** (`∂_{yᵢ}W_bg=σ_W∂_{ξᵢ}w₁`, `∂_{yᵢ}μ_R,bg=(μ̄_R/W̄₀)σ_W∂_{ξᵢ}m₁`). Truncation: "first order in
  wave amplitude and first shape order in each background bookkeeper."
- §2b/§2c: the two density representatives (RHO4-CONSTANT / RHOBR-CONSTANT) and the two anchorings
  (LAB_HELD `Q_bg^L(x,t)≡Q_bg(x)`, MATERIAL_ADVECTED `Q_bg^M(x,t)≡Q_bg(χ(x,t))`).
- §3a: `ρ_4D^α(x,t)≡ρ_4D,bg^{0,α}(x,t)(1+θ)`, `ρ_4D,bg^{0,L}≡ρ_4D,bg⁰(x)`,
  `ρ_4D,bg^{0,M}≡ρ_4D,bg⁰(χ(x,t))`, `μ_s^α≡μ_θ^α/ρ_br,bg^{0,α}`.
- §3c: every background face value/derivative in a trace "is obtained by differentiating a member of the
  supplied background state 𝔅⁰; none may be introduced as a free premise."

## Class 1 — the closure coefficient μ_θ
In TRACTION / CLOSURE_SHAPE_DERIV / VIRTUAL_WORK_SHAPE_DERIV, WL writes `μ_θ^α` as an **applied function**
`mu_theta_L(x1,x2,x3,time)` while PY writes a **bare symbol** `mu_theta_L`; the entire cross-engine residual
in those families is `Σ coeff·( mu_theta_L − mu_theta_L(x1,x2,x3,time) )`.
- **Rule on this:** does the spec endow `μ_θ^α` with in-plane (`x`) and/or `time` dependence that the T7 shape
  derivative must differentiate — in which case PY's bare symbol has **dropped** it (finding, PY-side) — or is
  `μ_θ^α` a per-anchoring constant for T7 (its evaluation-point arguments inert, so bare ≡ applied, benign)?
  Decide which from what §1c/§2a/§3a supply about `μ_θ`'s dependence and §3c's free-premise rule, and cite the
  spec sentences that force your answer.

## Class 2 — the projected density perturbation
In the PROJECTION families, the density-time term appears in PY as a single symbol `delta_rho_4D_bulk_t` and
in WL as an explicit expansion (`rho_br`, `W_0`/`W_bg`, `θ_t`, and — MATERIAL_ADVECTED — material velocity
`u_i_t`). The in-plane current divergence `∂_iδj_i` is present in **both** engines.
- **Rule on this:** per §3a `ρ_4D^α=ρ_4D,bg^{0,α}(1+θ)`, is PY's `delta_rho_4D_bulk` the spec's density
  perturbation `δρ_4D=ρ_4D,bg^{0,α}·θ`, i.e. the SAME object as WL's explicit expansion — for **each** of the
  four (anchoring × representative) combinations — including the MATERIAL_ADVECTED case, where the background
  `ρ_4D,bg⁰(χ(x,t))` is itself advected: does that background's time-dependence contribute to `δρ_4D`, and does
  each engine include that contribution? Or does one engine's density perturbation omit a term the spec
  requires (finding, and which engine)? Determine whether the spec's own `delta_rho_4D_bulk` (PY) is defined to
  carry the full advected time-dependence, or whether it must appear as a separate term.

## Method — MANDATORY (a prose derivation is discarded)
- Derive from the spec FIRST. Then **write a runnable SymPy script** that builds `μ_θ^α` and `δρ_4D` (all four
  combinations) from the spec and tests the two residual conditions; save the **script AND its literal stdout**
  to named absolute paths under `~/.s11_build/` (prefix `repr_adj_`) and report them. Prose-only derivations
  are discarded.
- Do NOT import or read any name map, reconciler, or the comparator source; build the objects yourself from the
  spec. You may inspect the committed operands in `~/.s11_build/comparator_run.out` only to confirm the residual
  form.
- **Physics filter:** a finding must catch a way the physics is wrong (a spec-required term genuinely dropped),
  not a code-style or different-input concern.
- **Deliverable per class:** verdict ∈ {IDENTITY (benign — same object, name the spec relation that makes bare
  ≡ applied / PY-symbol ≡ WL-expansion), FINDING (name which engine dropped which spec-required term, and the
  combination)}, the spec sentences you relied on, and the script + stdout paths. If your verdict differs
  between the four density combinations, say so per combination.
