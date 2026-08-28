# S11c-b — SHARED PHYSICS (variable-coefficient brane operator: the off-diagonal coupling kernel)

**Orchestrator-authored, 2026-08-27, under the S11c decision list (`directives/S11c_decisions.md`, N-series)
and `CLAUDE.md`.** This is the physics authority read by the two blind S11c-b engines: the SymPy engine
imports the closed S11b `LEDGER` (`scripts/S11b_exports.py`) and the S11c-a exports
(`scripts/S11c_a_exports.py`); the Wolfram engine imports nothing and reconstructs every object from the
supplied setup below and its two cited sibling specs. This document is an **obligation-to-compute
specification**, not a script and not a record of results. There is no acceptance value to withhold beyond
the value-free packages named below; there is no terminal `VERDICT`, `PASS`, or `FAIL`.

## 0 · Scope

S11c-b takes the non-uniform background, its two anchoring branches, and the first-order tilted-face
shape-derivative substrate produced by S11c-a, and constructs the **variable-coefficient brane (slab)
operator** and its **off-diagonal transverse→{θ,e_W,u_L} coupling kernel** — the first S11c physics object.
Concretely it asks for: the variable-coefficient stored-energy basis and any **new gradient-of-background
invariants** it admits (`N15`); the **divergence-form** slab operator obtained from that energy by S11b's
balance-law method with x-dependent coefficients; the off-diagonal block of that operator coupling the
transverse (photon) displacement sector to the thickness/longitudinal sector; and the **background
admissibility residual** reserved by S11c-a §2d (`N12`). Every object is multigraded by `(ε,η,σ_W)` and
carries restored physical dimension.

The following remain outside this build:

- S11c-c: solution of the curved two-face outgoing bulk problem and its nonlocal DtN/impedance/self-energy;
  no bulk-response solve (`δp = Z·v_bulk`), no permeability/memory kernel, no face-parity impedance object
  belongs here. S11c-b is slab-side geometry and the slab operator only (`N11a`);
- S11c-d: any profile-conditioned scattering, Born, Bloch, WKB, resonance, or spectral object; ⛔ **no global
  dispersion relation** `ω(k)` for a generic `W₀(x)` is requested (`N5`). S11c-b emits the coupling
  **operator/kernel**, not its spectrum;
- S11c-e: the flux-normalized dimensionless conversion form, the leakage observable, the confinement
  interpretation (`N13`), and the bench-optics comparison (`N7`);
- all nonlinear DC, harmonic, sideband, intensity, and soliton questions (`N10`/`G14`).

⛔ The background-flow correction `O(v_bulk_normal_0·|q_n|/ω)` is **not** an S11c-b task (`N11a`): S11c-b
inherits it as a standing rest-frame limit, and every emitted object is implicitly conditional on the
S11b-derived smallness domain `|q·v_bulk_normal_0/ω| ≪ 1`. `v_bulk_normal_0` is the inert rest-frame scope
limit of §1 and is never aliased to `v_0`.

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

### 1a · Inheritance, degrees of freedom, and the transverse/longitudinal split

The supplied setup, background ansatz, anchoring branches, face maps, level sets, orientation, measures,
trace law, dynamic window, and interface laws are **S11c-a §§1–3 imported unchanged**
(`directives/S11c_a_SHARED_PHYSICS.md`). The tilted-face first-order shape derivatives T-a..T-i are the
S11c-a §4 exports. The SymPy engine imports them from `scripts/S11c_a_exports.py`; the Wolfram engine
re-derives them from the S11c-a spec, importing nothing. This document amends and extends that setup with
the variable-coefficient constructions of §§3–5; where it does not amend, the S11c-a text governs.

The slab degrees of freedom are S11b's (`S11b_SHARED_PHYSICS.md:69–80`), with the internal slab fields
`{u,δW,θ}` distinguished from the two independent face variables `{ζ_+,ζ_-}` (S11c-a §3a):

```text
u(x,t)     in-plane displacement, three in-plane components, no w-component ;
θ(x,t)     Eulerian densification, ρ_4D = ρ_4D⁰(1+θ) ;
ζ_+, ζ_-   the two independent face variables, combined as
           δW ≡ ζ_+ − ζ_- (thickness) ,   ζ_c ≡ (ζ_+ + ζ_-)/2 (centre shift) ,
           ζ_s = ζ_c + s δW/2 ,   e_W ≡ δW/W₀ .
```

⛔ `ζ_c` is an independent face DOF; no S11c-b computation may set `ζ_c=0` or replace the two face variables
by a thickness-only ansatz (S11c-a §3a), except in an explicitly named centre-fixed uniform regression.

For the sector language used below, supply the in-plane displacement split

```text
u = u_L + u_T ,        ∇×u_L = 0 ,        ∇·u_T = 0 .
```

The split is an operator **label by the local differential structure**, not a global spectral projection:
`u_T` is the part acted on by `∇×`, `u_L` the part acted on by `∇·`. ⛔ No global Helmholtz projector,
plane-wave projector, or inverse-Laplacian is supplied or requested (`N5` removes the plane-wave setting).
On the **uniform** background `u` enters the stored energy only through its gradients (`∇×u`, `∇·u`) — the
translation-invariance consequence of §1c. ⚠ This is a uniform statement only: the non-uniform background
**breaks** in-plane translation invariance, so `u` (either part) may enter **undifferentiated** when
contracted with a background gradient — these are the `N15` spurion couplings §3a constructs, and are exactly
the channels this step exists to emit. Which parts a given operator block connects, and whether `u` enters
undifferentiated, is computed in §3a/§3c, not stated here. The label **transverse sector** means `u_T`; the
**{θ,e_W,u_L} sector** groups the two scalars with the longitudinal displacement.

`v_bulk_normal_0` (the bulk normal drain, `S11b_SHARED_PHYSICS.md:104`) is a scope-limit parameter, not an
active DOF, and appears in no derived operator (§0).

### 1b · Bulk acoustics and the projection law

Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
`∂_t²φ=c_s0²∇₄²φ`; the current and conservation law `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; and the
dynamic, anchored slab window `Ω` supplied in S11c-a §3. S11c-b performs no curved-bulk response solve (§0).

### 1c · Stored energy, constitutive definitions, symmetry inputs, and the balance-law method

Per unit projected in-plane volume, the supplied carried and kinetic energies are S11b's
(`S11b_SHARED_PHYSICS.md:257–263`), stated with the **constant** uniform coefficients:

```text
U = ½ μ_R |∇×u|² + ½ B_ρ⁽³⁾ θ² + C W₀ θ e_W
    + ½ k_W W₀² e_W² + ½ κ_W W₀² |∇(δW)|² ,
T = ½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_tδW)² ,
    with  B_ρ⁽³⁾ ≡ B_ρ W₀ ,   ρ_br⁰ ≡ ρ_4D⁰ W₀ .
```

⚠ The displayed five carried terms are **the set S11b carried, not a closed basis** — S11b's own basis
construction returned a larger independent set (a non-unique quotient; `S11b_interface_coupling_law.md:96–121`).
S11c-b constructs its own variable-coefficient basis in §3a and does not treat the displayed list as closed.

The uniform-basis inputs are the fields and their first gradients `{u,∇u,θ,∇θ,e_W,∇e_W}`, with the S11b
symmetry group in full (`S11b_SHARED_PHYSICS.md:280–288`): **in-plane translation invariance** (so `u` enters
only through its gradients, never undifferentiated), **in-plane `O(3)` isotropy and parity**, **reflection
`w→−w`**, **equivalence modulo total in-plane divergences**, and **no time-reversal, positivity, or
boundedness** assumption. Independence of invariants is judged as field bilinears with B1's constraint **not**
applied; every independent invariant is carried with its own free symbolic coefficient. These are the same
construction rules S11b task B0-energy used; §3a inherits them with variable coefficients.

⚠ **In-plane translation invariance and in-plane isotropy are properties of the UNIFORM background.** The
non-uniform profile breaks both. §3a inherits the construction METHOD but carries the background first jets
`{∂W_bg, ∂μ_R,bg}` as symmetry-breaking spurions, so undifferentiated-`u` couplings and anisotropic couplings
become admissible at first background-jet order (the `N15` new invariants). The reflection `w→−w`, the
total-in-plane-divergence quotient, and the no-time-reversal / no-positivity assumptions are unbroken —
but ⚠ the divergence quotient does **not** lift trivially to variable coefficients (§1d/§3a).

Constitutive derivatives are **variational**, not ordinary partial derivatives
(`S11b_SHARED_PHYSICS.md:349`):

```text
μ_θ ≡ (δU/δθ)|_{u, e_W, and all other fields fixed} ,   p_W ≡ (δU/δe_W)|_{…} .
```

The held-fixed qualifier is binding: `θ` may not be eliminated through a constraint before this derivative
is taken.

The supplied flat-face closure, affinity, response kernels, traction, kinematic balance, virtual-constraint,
sourced mass balance, and virtual work are exactly S11c-a §1c / §3b (`S11b_SHARED_PHYSICS.md:199–224,367–372,
642–645`):

```text
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
δ_vΣ_mat = 0 ,   (uniform linearisation)  δ_vθ + δ_ve_W + ∇_x·δ_vu = 0 .
```

The three `τ_I` are independent. Equations of motion are obtained by S11b's method — balance laws, the
binding virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed
external virtual work — **not** by putting an irreversible response kernel in an ordinary action.

### 1d · The uniform decoupling S11c-b generalizes (supplied context, not a target)

On the uniform background S11b found the transverse sector **decoupled** from the thickness sector: the mixed
second variation of the stored energy between `u_T` and `e_W` is inherited as the uniform baseline against
which S11c-b's variable-coefficient operator is compared in the §5c regression. Its uniform value and the
uniform transverse dispersion are S11b objects (`S11B_TRANSVERSE_COUPLING`, `S11B_TRANSVERSE_DISPERSION`,
`S11B_TRANSVERSE_DISSIPATION`); S11c-b imports/re-derives them only as the uniform-limit operand of §5c and
states nothing here about what the non-uniform operator returns.

⚠ The uniform transverse stiffness is **basis-representative-dependent** — an artifact of the non-unique
S11b energy-basis quotient (modulo total in-plane divergences). This is a property of the UNIFORM limit and
is reconciled ONLY in the §5c uniform-limit regression, under the pinned schema
(`steps/S11c_a_interface_shape_derivatives.md`, "reconciliation schema"). ⛔ It is **not** pre-registered as a
variable-coefficient comparator fold, and ⛔ the uniform quotient does **not** lift trivially to variable
coefficients: integrating by parts a variable coefficient generates first-background-jet terms
(`c∇·F ≡ −(∇c)·F` modulo a boundary term), so representatives that were equivalent uniformly differ by a
first-jet invariant that is **physics in the operator/kernel**, not a representational identity. §3a and the
comparator therefore treat any variable-coefficient representative difference as a computed object adjudicated
after the run, never as an inherited fold. No withheld uniform coefficient identity is stated here.

---

## 2 · Supplied non-uniform background ansatz and the activated admissibility computation

### 2a · The background ansatz — inherited

The background ansatz is S11c-a §2 imported unchanged: the constant bindings `W̄₀≡W_0`, `μ̄_R≡mu_R`; the
fresh varying profiles on the anchor coordinate `y`,

```text
ξ ≡ y/L_W ,   W_bg(y) ≡ W̄₀[1+η w₁(ξ)] ,   μ_R,bg(y) ≡ μ̄_R[1+η m₁(ξ)] ,   σ_W ≡ η W̄₀/L_W ,
∂_{yᵢ}W_bg = σ_W ∂_{ξᵢ}w₁ ,   ∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W ∂_{ξᵢ}m₁ ;
```

the two density representatives (`RHO4-CONSTANT`, `RHOBR-CONSTANT`, S11c-a §2b); the two physical anchoring
branches `LAB_HELD`/`MATERIAL_ADVECTED` (S11c-a §2c); and the background state and support bundle (S11c-a
§2d):

```text
𝔅⁰ ≡ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰, θ⁰, V_s⁰, J_s⁰, 𝒜_s⁰, boundary loads} ,   θ⁰ ≡ 0 ,
𝒮_hold⁰ ≡ {f_hold⁰(x), t_hold,s⁰(x)} ,   V_s⁰ = J_s⁰ = 𝒜_s⁰ = 0 .
```

**Which quantities vary (`N12`).** The varying background fields are exactly `W_bg`, `μ_R,bg`, and — through
§2b — the density representatives `ρ_4D,bg⁰`/`ρ_br,bg⁰`. Every explicit `W₀` factor in `U` and in the
bindings `B_ρ⁽³⁾≡B_ρ W₀`, `ρ_br⁰≡ρ_4D⁰ W₀` is the physical background thickness and becomes `W_bg(y)`. The
constitutive moduli `B_ρ`, `C`, `k_W`, `κ_W`, `μ_W` and the bulk constants `ρ_m`, `c_s0` do **not** vary
in-plane. ⛔ This geometric `W₀→W_bg(y)` promotion is legitimate but does **not** by itself define the
variable-coefficient energy: §3a constructs the full admissible basis, and the newly admitted
gradient-of-background invariants (`N15`) are **additional** to the promoted uniform terms — the substitution
is not the whole answer (`N15`).

The zero-jet contrast bookkeeper is `η`; the supplied first-jet bookkeeper is `σ_W`. They are varied
independently (by `η` and `L_W`); ⛔ no engine may replace `σ_W` by `η` or assign them a common order. Wave
perturbations carry the amplitude bookkeeper `ε`. **Every computed object is multigraded by `(ε,η,σ_W)` from
its actual data dependency**; the requested truncation is **first order in wave amplitude `ε` and first
shape order in each background bookkeeper `η` and `σ_W`**. ⛔ No term is removed merely because it carries
both a wave and a background bookkeeper, and no object's grade is stated here — it is computed and emitted.

### 2b · The activated admissibility computation (`N12`; reserved by S11c-a §2d)

S11c-a supplied `S11CA_BACKGROUND_STATE` and `S11CA_ADMISSIBILITY_PREMISE` and declared the profile
**support-stabilised** without testing stationarity. S11c-b performs the reserved test with the exact names
S11c-a fixed:

```text
S11CB_ADMISSIBILITY_OPERATOR_OPERAND  := the BACKGROUND-order (ε⁰, zero wave amplitude) balance residual of
                                          the variable-coefficient energy-and-geometry on 𝔅⁰ — the
                                          generalized body force and per-face traction the profile itself
                                          sources — expressed in the SAME ordered generalized-force pairing
                                          (bulk-DOF body force + per-face traction) as 𝒮_hold⁰ (§3d) ,
S11CB_ADMISSIBILITY_SUPPORT_OPERAND   := the declared support bundle 𝒮_hold⁰ = {f_hold⁰, t_hold,s⁰} in that
                                          same pairing ,
S11CB_ADMISSIBILITY_RESIDUAL          := operator operand − support operand (dimension- and pairing-matched) .
```

Emit operator operand, support operand, and their residual, per anchoring and density representative, with
multigrade and dimension. ⛔ This residual is **not** asserted zero and **not** guarded: it is the computed
statement of whether the declared support balances the variable-coefficient background, and a nonzero
residual is an admissible emitted outcome (a background the support does not hold would source spurious
coupling — the reason the test exists). ⛔ It is a **background-order** object (§3d): do **not** define it as
the `ε→0` limit of the §3b first-order wave operator — that operator is bilinear in the perturbations, so its
`ε→0` limit is identically zero and the test would be vacuous — and do **not** insert `W_bg−W_0` into the
uniform S11b perturbation equations (S11c-a §2d).

---

## 3 · The variable-coefficient operator construction

All definitions in §§1–2 are inputs. The expansions, coefficients, invariants, orders, blocks, and
cancellations below are **outputs** computed from them; ⛔ none is stated here.

### 3a · The variable-coefficient energy basis and its new invariants (`N15`)

Construct the stored-energy density admissible on the non-uniform background. The construction rule is
S11b's B0-energy, extended: enumerate the DOF fields and their first gradients `{u,∇u,θ,∇θ,e_W,∇e_W}` —
**one thickness coordinate only**, with the exact map `e_W,bg=(W_0/W_bg)e_W` (S11c-a §2a) imposed **before**
any independence/rank test so the thickness sector is not double-counted — **together with the supplied
background first jets** `{∂W_bg, ∂μ_R,bg}` treated as symmetry-breaking spurion data (the divergence-form
operator needs only first coefficient jets; no higher jet is introduced), and construct **every** scalar
bilinear in the DOF data allowed by the S11b symmetry group with those spurions carrying their
transformation. Independence is judged as field bilinears with B1's constraint not applied;
carry every independent invariant with its own free symbolic coefficient.

⭐ **Emit as RESULTS** the constructed basis, its count, and — separately — every invariant whose coefficient
carries a background first jet (a **new gradient-of-background invariant**), naming any new constant it
introduces. In-plane translation invariance and in-plane isotropy that a uniform background enjoyed are
broken by the profile; report which uniform invariants acquire a variable coefficient and which invariants
are newly admitted. ⛔ Do **not** obtain the variable-coefficient energy by the substitution `W₀→W_bg(y)`,
`μ_R→μ_R,bg(y)` into the uniform `U` — that smuggles the answer and omits the newly admitted terms; ⛔ and
do **not** forbid a symmetry-allowed new invariant. Whether any new invariant appears, and with what
coefficient, is the computed result.

```text
⇒ S11CB_ENERGY_BASIS_VARIABLE , S11CB_ENERGY_BASIS_COUNT , S11CB_ENERGY_BASIS_NEW_INVARIANTS ,
  S11CB_ENERGY_BASIS_OMISSIONS .
```

### 3b · The divergence-form variable-coefficient slab operator

From the §3a energy and S11b's balance-law method (§1c), compute the first-order equations of motion for the
slab DOFs `{u,θ,e_W}` on the non-uniform background, in **divergence form** — the variable coefficients sit
inside the in-plane divergences, so their first jets appear explicitly. Retain the full spatial dependence of
every background coefficient (`μ_R,bg`, `W_bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, and the `Σ_E⁰` map) and its first jet;
do not freeze a coefficient at its constant binding before differentiation. Use the S11c-a tilted-face
shape-derivative substrate (T-a..T-i) for every boundary/face contribution to the operator; the face and
flux terms are not recomputed here but consumed from that substrate. Emit the operator and, separately, the
computed provenance of each term (bulk-energy vs face/flux vs advective), per anchoring and density
representative, multigraded and dimensioned.

```text
⇒ S11CB_SLAB_OPERATOR , S11CB_SLAB_OPERATOR_TERM_ORIGINS ,
  S11CB_MU_THETA_OPERATOR  (the variable-coefficient μ_θ = δU/δθ operand, kept as a named operand,
                            not constructed by substitution into the uniform energy) .
```

### 3c · The off-diagonal coupling kernel

Extract the block of the §3b operator that couples the transverse structure (the `∇×u` part of the `u`
equations and its reciprocal) to the `{θ,e_W,∇·u}` structure — the object whose uniform limit is S11b's
decoupled zero (§1d). The sectors are the local differential-operator labels of §1a (curl vs divergence),
not a global spectral projection. Emit both off-diagonal blocks (transverse→thickness and thickness→
transverse) and their adjointness relation **with respect to the supplied variational pairing** (the
stored-energy/kinetic bilinear form of §1c, on the stated in-plane domain with the inherited face boundary
conditions), together with the kernel's `(ε,η,σ_W)` multigrade and dimension, per anchoring and density
representative. ⛔ Do not filter the kernel to a single channel and do not state its coefficient, sign,
parity, or grade — the operator's own block extraction is the computation.

```text
⇒ S11CB_COUPLING_KERNEL , S11CB_COUPLING_KERNEL_TERM_ORIGINS .
```

### 3d · The background-order balance (the admissibility operator operand)

Separately from the §3b first-order wave operator, compute the **background-order (ε⁰)** balance of the
variable-coefficient energy-and-geometry on `𝔅⁰`: the first variation of the total background functional with
respect to the brane configuration at `𝔅⁰`, yielding the generalized **body force and per-face traction** the
varying profile sources, in the **same ordered pairing** as `𝒮_hold⁰` (bulk-DOF body force + per-face
traction). This is the operator operand of §2b. Retain the full spatial dependence and first jets of `W_bg`,
`μ_R,bg`, and the `Σ_E⁰` map; keep both anchorings and both density representatives. ⛔ Do not reduce it to the
`ε→0` limit of the §3b wave operator, and do not insert `W_bg−W_0` into the uniform perturbation equations.
Whether the profile is stationary (residual against `𝒮_hold⁰` zero) or not is the computed result.

```text
⇒ S11CB_ADMISSIBILITY_OPERATOR_OPERAND  (the background-order balance operand of §2b) .
```

---

## 4 · Objects to compute and emit

Every item is computed for both anchoring branches and both density representatives wherever density enters,
carries the object's computed `(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no
component, term, order, coefficient, parity, survival, or cancellation in prose.

- **The variable-coefficient energy basis and new invariants** (§3a) ⇒ `S11CB_ENERGY_BASIS_VARIABLE`,
  `S11CB_ENERGY_BASIS_COUNT`, `S11CB_ENERGY_BASIS_NEW_INVARIANTS`, `S11CB_ENERGY_BASIS_OMISSIONS`.
- **The divergence-form slab operator** and its term provenance (§3b) ⇒ `S11CB_SLAB_OPERATOR`,
  `S11CB_SLAB_OPERATOR_TERM_ORIGINS`, `S11CB_MU_THETA_OPERATOR`.
- **The off-diagonal coupling kernel** (§3c) ⇒ `S11CB_COUPLING_KERNEL`,
  `S11CB_COUPLING_KERNEL_TERM_ORIGINS`.
- **The background admissibility package** (§2b) ⇒ `S11CB_ADMISSIBILITY_OPERATOR_OPERAND`,
  `S11CB_ADMISSIBILITY_SUPPORT_OPERAND`, `S11CB_ADMISSIBILITY_RESIDUAL`.
- **The representation-invariance package** (§5a) ⇒ `S11CB_REP_INVARIANCE_EULERIAN_OPERAND`,
  `S11CB_REP_INVARIANCE_MATERIAL_OPERAND`, `S11CB_REP_INVARIANCE_RESIDUAL`.
- **The one-sided independence control** (§5a) ⇒ `S11CB_CONTROL_INDEPENDENCE_BASE_OPERAND`,
  `S11CB_CONTROL_INDEPENDENCE_CORRUPTED_OPERAND`, `S11CB_CONTROL_INDEPENDENCE_RESIDUAL`.
- **The source-level form ablations** (§5b) ⇒ `S11CB_CONTROL_FORM_BASE_OPERAND`,
  `S11CB_CONTROL_FORM_ABLATED_OPERAND`, `S11CB_CONTROL_FORM_RESIDUAL`.
- **The uniform-limit regression** (§5c) ⇒ `S11CB_UNIFORM_LIMIT_S11CB_OPERAND`,
  `S11CB_UNIFORM_LIMIT_S11B_OPERAND`, `S11CB_UNIFORM_LIMIT_RESIDUAL`.
- **Dimensions and homogeneity** (§6) ⇒ `S11CB_DIMENSIONS`, `S11CB_HOMOGENEITY_BASE_OPERAND`,
  `S11CB_HOMOGENEITY_CONTROL_OPERAND`, `S11CB_HOMOGENEITY_RESIDUAL`.

No task may be replaced by a hand-typed expression carrying its anticipated content. Every comparison emits
operand A, operand B, and `A−B`; no residual value is supplied.

---

## 5 · Independent routes and controls

### 5a · Representation-invariance routes (`N4`/`N6`) — the genuine control

For each physical anchoring, compute the slab operator, the coupling kernel, and the admissibility operator
by two independent routes:

1. **direct Eulerian** graph/level-set linearization from the §§2–3 setup;
2. **material-coordinate** construction, using `x=x(X,t)` and the exact face-flattening coordinate of
   S11c-a §5a,

   ```text
   w′ = [w−ζ_c(X,t)]/[W_bg(x(X,t))+δW(X,t)]   (LAB_HELD) ,
   w′ = [w−ζ_c(X,t)]/[W_bg(X)+δW(X,t)]        (MATERIAL_ADVECTED) ,
   ```

   then mapped back to Eulerian variables (`N4`: `Δρ = δρ_E + u·∇ρ⁰`) before comparison.

The anchoring branch is held fixed across the two routes. Emit both operators and their difference under
`S11CB_REP_INVARIANCE_EULERIAN_OPERAND`, `S11CB_REP_INVARIANCE_MATERIAL_OPERAND`,
`S11CB_REP_INVARIANCE_RESIDUAL`, keyed by object and anchoring.

For the one-sided independence control (`N6`), mutate **only one route at its source** and recompute
downstream — the independence test between the two same-order channels (tilt, `N3`; advection, `N4`):

- replace `Σ_E(x(X,t),t)` by `Σ_E(X,t)` only in the direct-route material constraint feeding the operator; or
- reverse only the `x¹` first jet of `W_bg` in the upper-face direct-route source, leaving the
  material-coordinate route and every other source unchanged.

Emit the uncorrupted-route operand, the corrupted-route operand, and their residual under
`S11CB_CONTROL_INDEPENDENCE_BASE_OPERAND`, `S11CB_CONTROL_INDEPENDENCE_CORRUPTED_OPERAND`,
`S11CB_CONTROL_INDEPENDENCE_RESIDUAL`. The control does not authorize editing an already-emitted operator,
kernel, invariant, or residual.

### 5b · Source-level form ablations, one direction at a time

The two independent background profiles get **separate** one-source ablations — `w₁` and `m₁` are
independent with separate derivative maps (S11c-a §2a), so ablating both together cannot separate
thickness-slope from modulus-gradient coupling (the very channels `N6` distinguishes). For each direction
`i∈{1,2,3}` and each source `S∈{W_bg, μ_R,bg}` **separately**, create a formal first-jet ablation in which
only `∂_{xᶦ}S` is set to zero, **holding the other profile's jet fixed**, retaining the other two directions,
both anchorings, both density representatives, and every constitutive law and invariant; the density gradient
definitionally induced by `W_bg` (via §2b) co-varies only with the `W_bg` ablation. Recompute the §3a basis,
the §3b operator, and the §3c kernel from that source. ⛔ Do not ablate an already-computed operator, do not
drop all three directions simultaneously, and ⛔ do not use an `η` rescaling as this control (a coefficient
rescale tests arithmetic; only a form change tests physics). Emit the baseline operand, the independently
recomputed ablated operand, and their residual under `S11CB_CONTROL_FORM_BASE_OPERAND`,
`S11CB_CONTROL_FORM_ABLATED_OPERAND`, `S11CB_CONTROL_FORM_RESIDUAL`, keyed by
`{object,anchoring,density,source,direction}`.

### 5c · Uniform-limit regression smoke test

For the slab operator, the coupling kernel, and the transverse dispersion, independently obtain the
`(η,σ_W)→(0,0)` operand and the corresponding S11b uniform object (`S11B_TRANSVERSE_COUPLING`/
`S11B_TRANSVERSE_DISPERSION` and the uniform slab EOM). Emit both and their residual under
`S11CB_UNIFORM_LIMIT_S11CB_OPERAND`, `S11CB_UNIFORM_LIMIT_S11B_OPERAND`, `S11CB_UNIFORM_LIMIT_RESIDUAL`.

⛔ This is a **regression smoke test only** — it detects a forbidden gradient-**independent** term that
survives at `(η,σ_W)=0`. It does **not** validate the coupling's gradient coefficient, sign, or parity, and
`(η,σ_W)→0` is **not** an accepted corruption of the §5a/§5b controls (`N6`: a "coupling `∝∇W_bg` ⇒ vanishes
at `∇W_bg=0`" statement is the vacuous uniform limit renamed). The genuine controls are §5a and §5b.

---

## 6 · Method, dimensions, and script obligations

The derivation method is S11b's: balance laws, the binding virtual-displacement rule, variational
derivatives with held-fixed fields named, the supplied sign conventions, and prescribed external virtual
work. It is not an ordinary action-principle derivation of retarded response.

The structural rule is binding: physical symbols may be combined by hand only in the supplied setup,
background ansatz, face maps, energy, and supplied laws in §§1–3. Every §3–§5 expression is reached by
computation. Every control re-enters at an ansatz, map, level set, energy, or supplied law, **never at a
result**.

Restore units and compute the `[L,T,M]` dimension of every emitted object. `η`, `σ_W`, `w₁`, `m₁`, `θ`,
`e_W` are dimensionless; `W_0` and `L_W` have length dimension. For every homogeneity comparison emit both
operands and the residual, and include a source-level dimension corruption demonstrating that the check can
fail ⇒ `S11CB_DIMENSIONS`, `S11CB_HOMOGENEITY_BASE_OPERAND`, `S11CB_HOMOGENEITY_CONTROL_OPERAND`,
`S11CB_HOMOGENEITY_RESIDUAL`.

The three script clauses are exact obligations:

1. A script prints computed CAS objects and never states conclusions.
2. For every comparison it emits operand A, operand B, and `A−B` before any guard. A physics disagreement
   emits and continues with exit status 0; nonzero exit is reserved for operational failure only.
3. Interpretation belongs to the step record, not either engine.

Emission is never conditional on a payload's value; a boolean is a typed CAS object with its operands, never
a host-language native boolean. No script emits a terminal judgement, `VERDICT`, `PASS`, or `FAIL`.

---

## 7 · Names, F9 reservations, and parallel tag grammar

Every spatially varying object and every new kernel/invariant has a **fresh** standard name. The varying
names reserved by S11c-a are inherited unchanged (`W_bg`, `w1_profile`, `L_W`, `sigma_W`, `mu_R_bg`,
`m1_profile`, `rho_4D_bg_rho4_constant`, `rho_br_bg_rho4_constant`, `rho_4D_bg_rhobr_constant`, `e_W_bg`,
`eta_bg`). S11c-b's new objects — the variable-coefficient operator, the coupling kernel, the μ_θ operator,
and any new gradient-of-background constant emitted by §3a — take fresh names and ⛔ must **not** reuse the
imported constant keys `W_0`, `mu_R`, `rho_br`, `e_W`, or `v_0` (`N14`/`G3`: reusing a constant key for a
varying object proves a false `F9` equal and, for `rho_br`, silently freezes `∇Σ_E⁰=0` and drops the
advective channel). `v_bulk_normal_0` remains reserved for the drain and is never aliased to `v_0`. The
importing SymPy engine applies the inherited F9 collision check against the S11b and S11c-a keys; a
disagreement emits and continues.

Both engines use the exact grammar

```text
<ENGINE>_S11CB_<QUANTITY>
```

where `<ENGINE>` is `PY` or `WL`. Each base name written as `S11CB_<QUANTITY>` above is emitted by replacing
its leading `S11CB` with `PY_S11CB` or `WL_S11CB`; do not duplicate the `S11CB` component. Emit one tag per
named object; a single object's branch/face/DOF/density/direction cases are a keyed CAS map in that object's
payload, not separately invented tag names. Any unavoidable engine-local tag has `_LOCAL_` immediately after
`S11CB`, and each engine emits one local-tag inventory. The SymPy engine imports and carries the S11b
`LEDGER` and the S11c-a exports; the Wolfram engine re-derives the supplied §§1–3 inputs and the S11c-a
substrate without an import.

The S11c-b comparator is this sub-step's own frozen T7 join (`N8`, inherited verbatim from the S11b/S11c-a
contract — join by object name with the axis-typed keys of the S11c-a reconciliation schema, pair residual
operands, reject a native boolean as a residual operand, three-valued, repoint ablation). It computes and
prints, deciding nothing (rule 2). ⛔ No representational fold is pre-registered into the S11c-b comparator:
any cross-engine representative difference — including one inherited from the uniform basis quotient (§1d) —
is a computed residual **adjudicated after the run** by the hand-coded layer and the step record, as in
S11c-a, never comparator name-map surgery and never a pre-declared identity. ⚠ A pre-registered fold in a
variable-coefficient comparison could mask genuine first-background-jet coupling (§1d).

---

## 8 · Supplied versus computed; builder report

**Supplied and unfalsifiable here:** all of §1 (the inherited DOF, decomposition, energy, symmetry inputs,
and laws); the background ansatz, density representatives, anchoring definitions, background state, and
support premise of §2a (inherited from S11c-a §§1–3); and the S11c-a shape-derivative substrate T-a..T-i,
imported by SymPy and re-derived by Wolfram.

**Computed here:** the variable-coefficient energy basis and its new invariants (§3a); the divergence-form
slab operator and its term origins (§3b); the off-diagonal coupling kernel (§3c); the background
admissibility operator, support operand, and residual (§2b); and both operands and the residual for the
representation-invariance, one-sided independence, per-direction form-ablation, uniform-limit, and
homogeneity packages (§§5–6) — each with its computed multigrade and dimension.

The builder's report is under 25 lines and gives the files written, line and tag counts, tasks run or
skipped, runtime, all emitted tag names, and any ambiguity or non-computable requested object. It states
that §§1–2 and the S11c-a substrate were supplied and unfalsifiable in this build, and reports no computed
value and no conclusion about the physics.
