# S11c-a — SHARED PHYSICS (background & geometry: tilted-face shape derivatives)

**Codex-authored replacement, 2026-08-23, under `CLAUDE.md` rule 15.** This is the physics authority read by
the two blind S11c-a engines: the SymPy engine imports the closed S11b `LEDGER`; the Wolfram engine imports
nothing and reconstructs the same objects from the supplied setup below. This document is an
**obligation-to-compute specification**, not a script and not a record of results.

## 0 · Scope

S11c-a supplies a non-uniform background and its two physical anchoring branches, then asks for the
first-order shape derivatives of all S11b face geometry and interface laws. Its exports are background
maps, face kinematics, geometric boundary operators, and interface-law shape derivatives.

The following remain outside this build:

- S11c-b: the variable-coefficient brane energy/operator, its off-diagonal kernel, new
  gradient-of-background invariants (`N15`), and the background admissibility residual reserved in §2d;
- S11c-c: solution of the curved two-face outgoing bulk problem and its nonlocal DtN/impedance/self-energy;
- S11c-d: any profile-conditioned scattering, Bloch, WKB, resonance, or spectral object; no global
  dispersion relation is requested (`N5`);
- S11c-e: the flux-normalized conversion form, leakage observable, confinement question, and bench-optics
  comparison. The magnitude requiring the throat interior remains assigned to `R1` by
  `V3_STEP_PLAN.md:1179`;
- all nonlinear DC, harmonic, sideband, intensity, and soliton questions (`G14`).

There is no acceptance value to withhold in S11c-a. There is no terminal `VERDICT`, `PASS`, or `FAIL`.

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

Everything in this section is supplied from the closed uniform S11b physics. S11c-a does not test it.

### 1a · Coordinates, fields, face degrees of freedom, and conventions

There are three in-plane coordinates `x = (x¹,x²,x³)`, one normal coordinate `w`, and time `t`;
`D_brane = 3`. The uniform reference thickness is `W₀`. The slab has bulk on both sides.

- `u(x,t)` is a three-component in-plane displacement and has no `w` component. The material map is
  `x(X,t) = X + u(X,t)` with `𝒥_x = det(∂x/∂X)`.
- `ζ₊` and `ζ₋` are independent upper- and lower-face displacements, each measured along global `+w`.
  Their independent parity coordinates are
  `δW ≡ ζ₊ − ζ₋`, `ζ_c ≡ (ζ₊ + ζ₋)/2`, and hence `ζ_s = ζ_c + s δW/2` for `s ∈ {+1,−1}`.
  No S11c-a computation may set `ζ_c = 0` or replace the two face variables by a thickness-only ansatz.
- `θ` is the Eulerian fractional densification,
  `ρ_4D = ρ_4D⁰(1+θ)`. The uniform integrated inertia is
  `ρ_br⁰ ≡ ρ_4D⁰W₀`; the fractional thickness field is `e_W ≡ δW/W₀`.
- `μ_R` is the inherited twist modulus. The other slab inputs are
  `B_ρ`, `B_ρ⁽³⁾ ≡ B_ρW₀`, `μ_W`, `k_W`, `κ_W`, and `C`.
- `ρ_m` and `c_s0` are the bulk mass density and sound speed. `v_bulk_normal_0` is the steady bulk-normal
  drain and remains only the inherited rest-frame scope limit; the convective bulk problem is not reopened.

Only the perturbations may use harmonic time dependence. The uniform S11b plane-wave table is not an
ansatz for any background profile in §2.

For the flat reference faces `w = sW₀/2`, the outward normals are `n̂_s = s ŵ`. The outward flat-face
velocity and relative flux are

```text
V_s = s∂_tζ_s ,               J_s = ρ_m(v_bulk,w−∂_tζ_s)s .
```

Relative flux is positive when mass leaves the slab, is measured along the outward normal, and is per unit
true face area. The outward sign appears once.

### 1b · Bulk acoustics and the projection law

The supplied rest-frame bulk fields obey

```text
v_bulk = ∇₄φ ,              δp = −ρ_m ∂_tφ ,
∂_t²φ = c_s0² ∇₄²φ .
```

The supplied four-dimensional mass current and conservation law are

```text
j = ρ_4D v_bulk ,           ∂_tρ_4D + ∇₄·j = 0 .
```

`Ω` is a smooth slab window, approximately one inside and zero outside. S11b's projection object is the
result of integrating this conservation law against `Ω` and integrating by parts in `w`. S11c-a uses the
dynamic, anchored window supplied in §3; it may not replace it by a static window before differentiation.

The outgoing/decaying branch prescription and the dynamic-window projection problem are exactly S11b's.
S11c-a exports geometric boundary operators only; it does not compute a curved-bulk response.

### 1c · Stored energy, constitutive definitions, balance-law method

Per unit projected in-plane volume, the supplied carried energy and kinetic energy are

```text
U = ½ μ_R |∇×u|² + ½ B_ρ⁽³⁾ θ² + C W₀ θ e_W
    + ½ k_W W₀² e_W² + ½ κ_W W₀² |∇(δW)|² ,
T = ½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_tδW)² .
```

The uniform-basis construction uses the fields and first gradients
`{u,∇u,θ,∇θ,e_W,∇e_W}` with in-plane translation invariance, the full in-plane `O(3)` including parity,
reflection `w→−w`, and equivalence modulo total in-plane divergences. Time reversal, positivity, and
boundedness are not supplied assumptions. The blind engine reconstructs any uniform invariant it needs
from these inputs. S11c-a does not declare the displayed carried list closed and does not construct the
variable-coefficient basis or its new invariants; S11c-b owns that computation.

Constitutive derivatives are **variational**, not ordinary partial derivatives. In particular,

```text
μ_θ ≡ (δU/δθ)|_{u, e_W, and all other fields fixed} .
```

The held-fixed qualifier is binding: `θ` may not be eliminated through a constraint before this derivative
is taken.

The supplied flat-face closure, affinity, response kernels, traction, and kinematic balance are

```text
J_s = Λ_A(ω) 𝒜_s + Λ_V(ω) V_s ,
Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,                         I ∈ {A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,
μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,
n̂_s·v_bulk,s = V_s + J_s/ρ_m .
```

The three `τ_I` remain independent. `J_s`, `V_s`, and `t_s` are face-oriented quantities; the outward sign
is not applied a second time.

The slab material velocity used by the physical mass balance is

```text
v ≡ ∂_t u .
```

It is distinct from `v_bulk`. With `Σ_E ≡ Σ ≡ ρ_4D W`, the supplied physical evolution law is

```text
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋)                         (flat faces).
```

This sourced evolution law is separate from the supplied virtual-displacement law. At one instant a
virtual displacement transfers no mass through a face:

```text
Σ_mat(X,t) ≡ Σ_E(x(X,t),t) 𝒥_x(X,t) ,
δ_vΣ_mat(X,t) = 0 ,
δ_vx(X,t) = δ_vu(X,t) .
```

On the uniform background its dimensionless linearisation is
`δ_vθ + δ_ve_W + ∇_x·δ_vu = 0`. The binding object is `δ_vΣ_mat = 0`, never `δ_vΣ_E = 0`.

Equations of motion are obtained with balance laws, the virtual-displacement rule, and external virtual
work—not by putting an irreversible response kernel in an ordinary action. The supplied flat breathing
virtual work is

```text
δ_v𝒲_bulk = Σ_s t_s·δ_vx_s ,       δ_vx_s = n̂_s δ_v(δW)/2 .
```

S11c-a replaces the last face-displacement specialization by the branch-specific face maps in §3 and then
computes the resulting shape derivative. No extra direct `J_s` generalized force is supplied.

## 2 · Supplied non-uniform background ansatz

### 2a · Constant bindings, contrast, and independent profile scale

The constant references in the ansatz are the imported S11b constants:

```text
W̄₀ ≡ W_0 ,                    μ̄_R ≡ mu_R .
```

`W_0` and `mu_R` remain constant ledger keys. Define the fresh varying profiles on an anchor coordinate
`y` by

```text
ξ ≡ y/L_W ,
W_bg(y)   ≡ W̄₀ [1 + η w₁(ξ)] ,
μ_R,bg(y) ≡ μ̄_R [1 + η m₁(ξ)] ,
σ_W       ≡ η W̄₀/L_W .
```

`w₁` and `m₁` are dimensionless `O(1)` profile functions; `L_W` is an independent length. The supplied
first-derivative map is

```text
∂_{yᵢ}W_bg = σ_W ∂_{ξᵢ}w₁ ,
∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W ∂_{ξᵢ}m₁ .
```

Thus the zero-jet contrast bookkeeping uses `η`, while the supplied first-jet bookkeeping uses `σ_W`.
`η` and `σ_W` are varied independently by varying `η` and `L_W`; no engine may replace `σ_W` by `η` or
assign a common order to them. These are grades of the supplied ansatz, not preassigned grades of normals,
tractions, traces, window derivatives, or any other output.

Wave perturbations `u`, `ζ_c`, `δW`, `θ`, and the first variations of the face and bulk fields carry the
independent amplitude bookkeeper `ε`. Every computed object must be multigraded by `(ε,η,σ_W)` from its
actual data dependency. No term is removed merely because it contains both a wave and a background
bookkeeper; the requested truncation is first order in wave amplitude and first shape order in each
background bookkeeper.

The local thickness variable and affinity normalization are supplied through the exact maps

```text
e_W,bg(y) ≡ δW/W_bg(y) = [W_0/W_bg(y)] e_W ,
μ_s(y)    ≡ μ_θ/ρ_br,bg⁰(y) .
```

These maps remain present in every branch and density representative until after differentiation. Their
bookkeeping is computed from the ansatz; this document assigns no resulting term or order.

### 2b · The two density representatives (`N11b`)

Let `ρ_4D,ref⁰ ≡ rho_br/W_0`. Carry both supplied representatives, using fresh names for varying fields:

```text
RHO4-CONSTANT:
    ρ_4D,bg⁰(y) ≡ ρ_4D,ref⁰ ,
    ρ_br,bg⁰(y) ≡ ρ_4D,bg⁰(y) W_bg(y) ;

RHOBR-CONSTANT:
    ρ_br,bg⁰(y) ≡ rho_br ,
    ρ_4D,bg⁰(y) ≡ rho_br/W_bg(y) .
```

For each representative, construct `Σ_E⁰(y) ≡ ρ_4D,bg⁰(y)W_bg(y)` and compute its full in-plane gradient.
Do not simplify either representative into the other before T-g, T-h, or T-i. Emit the two computed maps
under `S11CA_BACKGROUND_DENSITY_MAP`.

### 2c · Two physical anchoring branches

Let `χ(x,t)` be the inverse material map, defined only by `χ(x(X,t),t)=X`. For every background profile
`Q_bg ∈ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰}`, the two supplied physical branches are

```text
LAB_HELD:          Q_bg^L(x,t) ≡ Q_bg(x) ,
MATERIAL_ADVECTED: Q_bg^M(x,t) ≡ Q_bg(χ(x,t)) .
```

These are distinct physical anchorings, not alternate Eulerian/material representations of one branch.
Both branches must be computed and emitted separately. Their single permitted face parameterizations are
given in §3; an engine may not substitute another parameterization.

### 2d · Background state, support, and the reserved admissibility computation

For each anchoring and density representative, name the supplied background state

```text
𝔅⁰ ≡ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰, θ⁰, p_s⁰, J_s⁰, boundary loads} .
```

Because `θ` is defined relative to the selected local background density, set `θ⁰ ≡ 0`; the spatial
background density is already carried by `ρ_4D,bg⁰` and is not counted again through `θ⁰`.
The profile is declared **externally held** by the named support bundle
`𝒮_hold⁰ ≡ {f_hold⁰(x), t_hold,s⁰(x)}`. The face backgrounds `p_s⁰(x)` and `J_s⁰(x)` are not set to zero,
and their in-plane gradients are not excluded. Emit the state and premise, unchanged, as supplied objects:

```text
S11CA_BACKGROUND_STATE
S11CA_ADMISSIBILITY_PREMISE
```

S11c-a does not test stationarity. Reserve the can-fail S11c-b comparison, with these exact names:

```text
S11CB_ADMISSIBILITY_OPERATOR_OPERAND  := the S11c-b variable-coefficient stationary operator on 𝔅⁰,
S11CB_ADMISSIBILITY_SUPPORT_OPERAND   := the declared support 𝒮_hold⁰,
S11CB_ADMISSIBILITY_RESIDUAL          := operator operand minus support operand.
```

Those three objects are not emitted or guarded in S11c-a. In particular, do not insert
`W_bg−W_0` into the uniform S11b perturbation equations as an admissibility test.

## 3 · Supplied face geometry and laws for the shape derivative

All definitions in this section are inputs. Their expansions, component pairings, coefficients, orders,
and cancellations are outputs of §4.

### 3a · One parametric face map per anchoring, with both face degrees of freedom

At fixed material label `X`, use exactly

```text
R_s^L(X,t) = ( x(X,t),  s W_bg(x(X,t))/2 + ζ_s(X,t) ) ,
R_s^M(X,t) = ( x(X,t),  s W_bg(X)/2        + ζ_s(X,t) ) ,
ζ_s(X,t)   = ζ_c(X,t) + s δW(X,t)/2 .
```

Emit these two supplied maps without alteration under `S11CA_FACE_MAP_LAB_HELD` and
`S11CA_FACE_MAP_MATERIAL_ADVECTED`, respectively.

The corresponding Eulerian graph heights and level sets are fixed by the same maps:

```text
h_s^L(x,t) = s W_bg(x)/2       + ζ_s(χ(x,t),t) ,
h_s^M(x,t) = s W_bg(χ(x,t))/2  + ζ_s(χ(x,t),t) ,
F_s^α(x,w,t) ≡ w − h_s^α(x,t) ,                 α ∈ {L,M} .
```

The slab interior is selected by `sF_s^α < 0`. The outward unit normal `n̂_s^α` is the unit vector normal
to `F_s^α=0` satisfying the binding orientation law

```text
s(n̂_s^α·ŵ) > 0 .
```

The bare sign of `∇₄F_s^α` is not an orientation prescription. Define the graph measure per projected
in-plane area by

```text
a_s^α(x,t) ≡ sqrt(1 + |∇_x h_s^α(x,t)|²) .
```

### 3b · Face displacement, velocity, flux, traction, and balance

Write each total face pressure as `p_s^α = p_s^{0,α} + δp_s^α`; the first variation `δp_s^α` is the
quantity denoted `δp_s` in the uniform S11b law. Decompose every other face field about its declared §2d
background in the same way, but obtain the variation of `J_s^α` from the single relative-flux definition
below rather than introducing an independent flat flux.

For each branch, the only permitted virtual face displacement and face velocity vector are obtained from
the supplied parameterization:

```text
δ_vx_s^α ≡ δ_vR_s^α|_X ,              v_face,s^α ≡ ∂_tR_s^α|_X ,
V_s^α ≡ V_{n,s}^α ≡ n̂_s^α·v_face,s^α .
```

Use the same `V_s^α` in every object below. With all bulk quantities traced at `R_s^α`, define once

```text
J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α ,
n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m ,
𝒜_s^α ≡ μ_s^α − p_s^α/ρ_m ,
t_s^α ≡ −(p_s^α + Λ_X𝒜_s^α)n̂_s^α ,
J_s^α − Λ_A𝒜_s^α − Λ_VV_s^α = 0 .
```

`J_s^α` is per unit true face area. T-c, T-c′, T-h, and T-i must use this one object; no flat Cartesian
`J_s` or separately reconstructed `V_s` may appear in any of them.

Per unit projected in-plane area, the supplied face contribution to virtual work and the supplied physical
mass balance are

```text
δ_v𝒲_bulk^α ≡ Σ_s a_s^α t_s^α·δ_vx_s^α ,
∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α ,       v = ∂_t u .
```

Neither `a_s^α` may be replaced by one before the shape derivative is taken. The background values
`p_s⁰` and `J_s⁰` remain in these laws during differentiation.

### 3c · Shifted traces and the dynamic window

For any bulk field `f(x,w,t)` with a nonzero background face value or derivative, use the trace
linearisation law

```text
δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰) .
```

Apply it to every traced bulk field used by T-c through T-i; do not evaluate first at `w=sW_0/2` and then
discard the face shift.

Let `G_s^α ≡ sF_s^α`. Supply one fixed smooth two-argument window function `𝒪` and define

```text
Ω^α(x,w,t) ≡ 𝒪(G_+^α(x,w,t), G_-^α(x,w,t)) ,
```

with `𝒪` approximately one when both arguments are negative and tending to zero outside. This fixes the
window anchoring from the same face maps. Its time dependence is retained until after the projection
identity and its shape derivative have been computed.

## 4 · Objects to compute and emit

Every item is computed for both anchoring branches, both faces, both independent face DOFs
`{δW,ζ_c}`, and both density representatives wherever density enters. Every payload carries the object's
computed `(ε,η,σ_W)` multigrade and restored physical dimension. Do not state a component, term, order,
coefficient, parity, survival, or cancellation in prose.

- **T-0 · Background density map:** construct `Σ_E⁰` and its in-plane gradient from §2b.
  ⇒ `S11CA_BACKGROUND_DENSITY_MAP`.
- **T-a · Outward face normal:** the oriented unit-normal object from §3a.
  ⇒ `S11CA_FACE_NORMAL`.
- **T-a′ · Conormal:** the face operator `n̂_s·∇₄`, including evaluation on the supplied graph.
  ⇒ `S11CA_CONORMAL_DERIV`.
- **T-a″ · True-area measure:** the face measure and its full shape derivative from §3a.
  ⇒ `S11CA_FACE_MEASURE_SHAPE_DERIV`.
- **T-b · Outward face velocity:** the normal-velocity object obtained from `R_s^α` in §3b.
  ⇒ `S11CA_FACE_VELOCITY`.
- **T-c · Relative flux:** the full shape derivative of the one relative-flux law in §3b, without a
  component filter.
  ⇒ `S11CA_RELATIVE_FLUX`.
- **T-c′ · Interfacial kinematic balance:** the shape derivative of the bound kinematic law in §3b.
  ⇒ `S11CA_KINEMATIC_BALANCE`.
- **T-d · Traction and virtual work:** compute the traction object and the shape derivative of
  `δ_v𝒲_bulk^α`, taking `δ_vx_s^α` only from the applicable `R_s^α` and differentiating its face measure.
  Which virtual-displacement pairings occur is part of the computation.
  ⇒ `S11CA_TRACTION`, `S11CA_VIRTUAL_WORK_SHAPE_DERIV`.
- **T-e · Shifted-face trace:** the face-evaluation operator obtained from §3c for every bulk field consumed
  downstream.
  ⇒ `S11CA_FACE_SHIFT`.
- **T-f · Dynamic-window projection:** the shape derivative of S11b's projection identity using `Ω^α`,
  together with the independently computed static-flat projection operand, their residual, and a computed
  provenance inventory for every term in the dynamic and static operands.
  ⇒ `S11CA_PROJECTION_SHAPE_DERIV`, `S11CA_PROJECTION_STATIC_OPERAND`,
  `S11CA_PROJECTION_DYNAMIC_OPERAND`, `S11CA_PROJECTION_RESIDUAL`,
  `S11CA_PROJECTION_TERM_ORIGINS`.
- **T-g · Virtual material constraint:** linearise the dimensionless object
  `δ_vΣ_mat^α/Σ_mat^{0,α}` about the non-uniform background, using the virtual `δ_vu`, the selected density
  representative, and the exact `e_W,bg` map. Do not substitute the physical `u` for `δ_vu`; do not add a
  dimensionful term to the dimensionless object; do not remove `e_W,bg` in either density representative.
  ⇒ `S11CA_VIRTUAL_CONSTRAINT`.
- **T-h · Physical sourced mass balance:** linearise the true-area balance in §3b about the declared
  background using the slab velocity `v=∂_tu`, the same T-c fluxes, and the full spatial dependence of
  `Σ_E⁰`. Preserve the provenance of every computed term from the exact balance.
  ⇒ `S11CA_EVOLUTION_MASS_BALANCE`, `S11CA_EVOLUTION_TERM_ORIGINS`.
- **T-i · Face-closure shape derivative:** shape-differentiate only the supplied face closure in §3b,
  using T-c's `J_s`, T-b's `V_s`, the shifted pressure trace, and the selected
  `μ_s=μ_θ/ρ_br,bg⁰` map. Keep the variable-coefficient `μ_θ` operator as S11c-b's named operand rather
  than constructing it by substitution into the uniform energy. This object is **not** B0c's bulk-response assembly
  `δp=Z·v_bulk`; no bulk DtN, impedance, or pressure-response solve belongs in T-i.
  ⇒ `S11CA_CLOSURE_SHAPE_DERIV`.

No task may be replaced by a hand-typed expression carrying its anticipated content.

## 5 · Independent routes and controls

Every comparison emits operand A, operand B, and the computed residual. No residual value is supplied.

### 5a · Representation-invariance routes (`N4`/`N6`)

For each physical anchoring, compute T-g, T-c, T-d, and T-i by two independent routes:

1. direct Eulerian graph/level-set differentiation from §§2c and 3;
2. material-coordinate differentiation using `x=x(X,t)` and the exact face-flattening coordinate

   ```text
   w′ = [w−ζ_c(X,t)]/[W_bg(x(X,t))+δW(X,t)]       (LAB_HELD),
   w′ = [w−ζ_c(X,t)]/[W_bg(X)+δW(X,t)]            (MATERIAL_ADVECTED).
   ```

Map route 2 back to the Eulerian variables before comparison. The anchoring branch is held fixed across
the two routes. For every compared T-object emit the two operators and their difference under
`S11CA_REP_INVARIANCE_EULERIAN_OPERAND`, `S11CA_REP_INVARIANCE_MATERIAL_OPERAND`, and
`S11CA_REP_INVARIANCE_RESIDUAL`, keyed by object and anchoring.

For the one-sided independence control, mutate only one route at its source and recompute downstream:

- for T-g, replace `Σ_E(x(X,t),t)` by `Σ_E(X,t)` only in the direct-route definition of `Σ_mat`;
- for T-c, T-d, and T-i, reverse only the `x¹` first jet of `W_bg` in the upper-face direct-route
  `F_+^α`/`R_+^α` source, leaving the material-coordinate route and every other source unchanged.

Emit the uncorrupted route operand, corrupted route operand, and their residual under
`S11CA_CONTROL_INDEPENDENCE_BASE_OPERAND`, `S11CA_CONTROL_INDEPENDENCE_CORRUPTED_OPERAND`, and
`S11CA_CONTROL_INDEPENDENCE_RESIDUAL`. The control does not authorize editing an already-emitted normal,
traction, work term, closure term, or residual.

### 5b · C-1 source-level form ablations, one direction at a time

For each `i ∈ {1,2,3}` separately, create a formal first-jet ablation of the supplied face source in which
only `∂_{xᶦ}W_bg` is set to zero in `F_s^α` and the corresponding `R_s^α`; retain the other two components,
both faces, both face DOFs, all background face fields, and every constitutive law. Recompute T-a through
T-i from that source, including T-d and T-i. Do not ablate an already-computed `n̂_s`, do not drop all three
directions simultaneously, and do not use an `η` rescaling as this control.

For every T-object and each direction emit the baseline operand, the independently recomputed ablated
operand, and their residual under `S11CA_CONTROL_FORM_BASE_OPERAND`,
`S11CA_CONTROL_FORM_ABLATED_OPERAND`, and `S11CA_CONTROL_FORM_RESIDUAL`, keyed by
`{object,anchoring,face,direction}`.

### 5c · Uniform-limit regression smoke test

For each S11c-a object, independently obtain the `(η,σ_W)→(0,0)` operand and the corresponding S11b object.
Emit both and their residual under `S11CA_UNIFORM_LIMIT_S11CA_OPERAND`,
`S11CA_UNIFORM_LIMIT_S11B_OPERAND`, and `S11CA_UNIFORM_LIMIT_RESIDUAL`. This is a regression smoke test,
not the form control of §5b.

## 6 · Method, dimensions, and script obligations

The derivation method is S11b's: balance laws, the binding virtual-displacement rule, variational
derivatives with held-fixed fields named, the supplied sign conventions, and prescribed external virtual
work. It is not an ordinary action-principle derivation of retarded response.

The structural rule is binding: physical symbols may be combined by hand only in the supplied setup,
background ansatz, face maps, and supplied laws in §§1–3. Every §4–§5 expression is reached by computation.
Every control re-enters at an ansatz, map, level set, or supplied law, never at a result.

Restore units and compute the `[L,T,M]` dimension of every emitted object. `η`, `σ_W`, `w₁`, and `m₁` are
dimensionless; `W_0` and `L_W` have length dimension. For every homogeneity comparison, emit both operands
and the residual, and include a source-level dimension corruption demonstrating that the check can fail.
⇒ `S11CA_DIMENSIONS`, `S11CA_HOMOGENEITY_BASE_OPERAND`, `S11CA_HOMOGENEITY_CONTROL_OPERAND`,
`S11CA_HOMOGENEITY_RESIDUAL`.

The three script clauses are exact obligations:

1. A script prints computed CAS objects and never states conclusions.
2. For every comparison it emits operand A, operand B, and `A−B` before any guard. A physics disagreement
   emits and continues with exit status 0; nonzero exit is reserved for operational failure only.
3. Interpretation belongs to the step record, not either engine.

Emission is never conditional on a payload's value. A boolean remains a typed CAS object with its operands;
it is not serialized as a host-language native boolean. No script emits a terminal judgement.

## 7 · Names, F9 reservations, and parallel tag grammar

Every spatially varying object has a fresh name. Reserved varying names include `W_bg`, `w1_profile`,
`L_W`, `sigma_W`, `mu_R_bg`, `m1_profile`, `rho_4D_bg_rho4_constant`,
`rho_br_bg_rho4_constant`, `rho_4D_bg_rhobr_constant`, `e_W_bg`, `eta_bg`, and the later S11c kernels and
observables. They must not reuse the imported keys `W_0`, `mu_R`, `rho_br`, `e_W`, or `v_0`.
`v_bulk_normal_0` remains reserved for the drain and is never aliased to `v_0`. The constant bindings in
§2a do not authorize either varying profile to use a constant key. The importing engine applies the
inherited F9 collision check; a disagreement emits and continues.

Both engines use the exact grammar

```text
<ENGINE>_S11CA_<QUANTITY>
```

where `<ENGINE>` is `PY` or `WL`. Each base name written as `S11CA_<QUANTITY>` in §§2–6 is emitted by
replacing its leading `S11CA` with `PY_S11CA` or `WL_S11CA`; do not duplicate the `S11CA` component. Emit
one tag per named object. A single object's branch/face/DOF/density/direction cases are a keyed CAS map in
that object's payload, not separately invented tag names. Both engines emit parallel tag sets. Any
unavoidable engine-local tag has `_LOCAL_` immediately after `S11CA`, and each engine emits one local-tag
inventory. The SymPy engine imports and carries the S11b `LEDGER`; the Wolfram engine re-derives the
supplied §§1–3 inputs without an import.

## 8 · Supplied versus computed; builder report

**Supplied and unfalsifiable here:** all of §1; the constant bindings, profile ansätze, derivative-order
map, density representatives, anchoring definitions, background state, and support premise of §2; and all
face maps, level sets, orientation, measures, trace law, dynamic window, and interface laws of §3. The
admissibility state/premise are supplied, not tested; their can-fail comparison is the explicitly reserved
S11c-b three-object package in §2d.

**Computed here:** every T-object in §4, including its multigrade and dimension; both operands and residual
for the representation, projection, homogeneity, per-direction form-ablation, one-sided independence, and
uniform-limit packages in §§4–6.

The builder's report is under 25 lines and gives the files written, line and tag counts, tasks run or
skipped, runtime, all emitted tag names, and any ambiguity or non-computable requested object. It states
that §§1–3 and the admissibility premise were supplied and unfalsifiable in this build. It reports no
computed value and no conclusion about the physics.
