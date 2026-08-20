# S11b — SHARED PHYSICS  (the interface coupling law)

⛔⛔ **This file is read by BOTH engines. An error here defeats dual-engine agreement by construction: both
engines will compute the same wrong thing and agree.** It is reviewed as its own artifact — its own two
legs — before either build starts.

⭐ It supplies **the setup, the field content, the governing equations, the supplied constitutive laws and
derivation routes, the premises, and the list of quantities to compute.** ⛔ It supplies **no result** — no
face response, no impedance, no dispersion root, no sign, no stability boundary, no coefficient of any
coupling, no dimension tuple, no admissibility region, and **no reduction value any check is expected to
reproduce**. Where a check compares your computation against an independently-derived reference, the
reference is **held by the orchestrator and diffed on our side** — ⛔ it does not appear here.

⚠ **This is ONE unified step on the export chain**, subsuming the two historical execution stages S11b-A
(the bulk's response to moving faces) and S11b-B (the homogeneous interface assembly). ⭐ There is one
shared spec, two blind engines, one comparator, one step record, one card. ⛔ The **non-uniform**
variable-coefficient problem — whether light's confinement is unconditional — is a **separate later step
(S11b-C)** and is not built here.

---

## 0 · What this step is, and what it is not

Assemble the bulk's response to the slab's moving faces, the brane's in-plane sector, and the slab's
thickness degree of freedom into **one linear system on a UNIFORM background**, and determine the
**longitudinal mode's fate**: does it propagate freely, decay, **grow**, or fail to exist as a normal mode.
All background quantities are uniform constants.

⭐⭐ **DECAY AND GROWTH ARE ADMISSIBLE OUTCOMES AND MUST EACH BE REPORTED AS SUCH.** The moduli in §5 are
free symbols and ⛔ **no boundedness condition is imposed on them or on the cross term `C`**, so the
quadratic form is not assumed positive-definite and a root with `Im ω > 0` is a possible result of this
calculation. ⛔ **If you find a growing or decaying root, do not discard it, do not re-branch to remove it,
and do not add a stability assumption.** Report it, and report the condition on the moduli that separates
the two outcomes.

⚠⚠ **A GROWING OR DECAYING ROOT CAN ALSO BE MANUFACTURED BY A MISTAKE.** A derivation-route error (§6) can
turn a sink into a source; the branch re-selection error (§1b) can turn a damped resonance into an apparent
instability; each mirror error manufactures the opposite. Each mistake pair has a **named, mechanical
diagnostic** in its section. ⭐ **Run both diagnostics before reporting any root with `Im ω ≠ 0`, and report
their outcomes alongside every such root.** ⛔ This is **not** permission to discard either outcome —
distinguish each reported non-real root from an artifact.

⛔⛔ **IN SCOPE, and settled HERE** (uniform background): the bulk face response and its three regimes
including the grazing threshold; the permeable-face constitutive response under the supplied closure; the
assembled longitudinal mode's fate and the behaviour of the coupled system **at and near the sound-cone
grazing threshold** (`q_out → 0`); the transverse mode's coupling on this uniform background; the
passivity, reciprocity and admissibility regions; and the breathing-mode stability question.

⛔ **DEFERRED to S11b-C** (non-uniform only): the full variable-coefficient slab spectrum; actual leakage
rates in the non-uniform problem; and whether light's confinement is **unconditional**. ⭐ Task **B6**
computes the transverse coupling on this uniform background; whatever it returns, ⛔ that result does **not**
settle the unconditional question either way, and ⛔ no tag here may phrase it as if it does.

⛔ **DEFERRED to the nonlinear light programme** (⛔ **not** S11b-C): the DC / harmonic / sideband radiation
audit and any nonlinear intensity coupling. Nothing here is linear-in-amplitude evidence about those.

⛔ **Do not consult any file for the answers. Everything needed is below.** If something needed is missing,
emit the tag as `NOT_ESTABLISHED` and name what is missing — ⛔ do not go looking. A refusal is a valid and
valuable output.

---

## 1 · Geometry, fields, conventions

Four spatial dimensions `(x¹, x², x³, w)` and time `t`; `x = (x¹, x², x³)`, `D_brane = 3`. A slab of
thickness `W(x,t) = W₀ + δW` centred on `w = 0`, faces at `w = ±W₀/2`. `w` is normal to the slab. Bulk on
**both** sides, `|w| > W₀/2`.

| symbol | meaning |
|---|---|
| `u(x,t)` | in-plane displacement, **3 components**, ⛔ no `w`-component |
| `δW(x,t)` | thickness perturbation; faces displace to `ζ_± = ±δW/2` |
| `θ(x,t)` | **densification** — the **Eulerian** fractional perturbation of the brane material's local 4D density, `ρ_4D = ρ_4D⁰(1 + θ)`. ⚠ Eulerian, ⛔ not material |
| `e_W ≡ δW/W₀` | fractional thickness change |
| `ρ_br⁰ ≡ ρ_4D⁰ W₀` | background slab-integrated inertia density — ⭐ **see §11 for the mapping to the imported `rho_br`** |
| `μ_R` | curl-type (twist) modulus — ⭐ **an inherited input, see §11** |
| `B_ρ` | local 4D compression modulus; `B_ρ⁽³⁾ ≡ B_ρ W₀` |
| `μ_W`, `k_W`, `κ_W` | thickness inertia, restoring stiffness, gradient stiffness — **symbolic** |
| `C` | the symmetry-allowed cross term between densification and thickening (§5) |
| `ρ_m`, `c_s0` | bulk mass density and sound speed — ⭐ `c_s0` **inherited (§11)**; `ρ_m` **originates in this step** |

**Conventions — fixed here so the two engines produce comparable signs.** ⚠ Reported real and imaginary
parts are meaningless for comparison unless these are shared.

| | convention |
|---|---|
| harmonic dependence | every perturbation `∝ exp[i(k·x − ωt)]`, `k` real |
| face displacements | `ζ₊`, `ζ₋` = displacement of the **upper / lower** face, **both measured along global `+w`** |
| parity combinations | **thickness** `δW ≡ ζ₊ − ζ₋` · **centre shift** `ζ_c ≡ (ζ₊ + ζ₋)/2` |
| outward normals | upper face `n̂₊ = +ŵ` · lower face `n̂₋ = −ŵ` |
| outward face velocity | `V_± ≡ (∂_t ζ_±)·(±1)` — ⭐ **face-odd**; use this, ⛔ not the global `∂_t ζ_±`, wherever a quantity is along the outward normal |
| response ratio | `Z ≡ (bulk pressure perturbation at a face) / (bulk material's OUTWARD normal velocity at that face)` |
| ⭐⭐ radiation condition | in each half-space retain **only** waves carrying energy **away** from the slab (real normal wavenumber) or **decaying** away from it (imaginary). ⛔ **There are NO incoming waves from infinity.** State explicitly the sign of the normal wavenumber selected in each half-space under the harmonic convention above. ⚠ Boundedness alone does **not** select a branch when the normal wavenumber is real. |

⛔ `ζ_c` is **NOT** the `h`-branon: the ledger's `h` is dimensionless (`ξ_w = ℓh`); `ζ_c` here is a
**length**. They are different objects and must never be identified without the normalisation relating them.
⛔ `x` is a position, ⛔ never a displacement.

### ⭐⭐ The background drain — a NAMED quantity, distinct from any in-plane background velocity

There is a **steady background transfer** of material across the interface. Let **`v_dr`** be the resulting
steady **normal** flow speed in the bulk near a face (**symbolic**). ⭐ **`v_dr` is the bulk's NORMAL drain
across a face; it couples to the bulk NORMAL wavenumber `q_out`.** Its standard ledger name is
`v_bulk_normal_0`. ⛔⛔ **It is a DIFFERENT physical quantity from any in-plane / tangential background
brane velocity, and it must NOT share that quantity's name or key.** ⛔ This step does not set the brane
moving in-plane; `v_dr` here is neither that quantity nor an override of any premise about it.

⚠ **`v_dr`'s only role in this step is as a scope limit.** §2's bulk acoustics are linearised about **rest**;
a nonzero `v_dr` adds convective terms and a first-order flux. ⭐ Both engines carry the rest-frame operators
and **record the discarded convective correction as a stated scope limit** (Task **B9**);
⛔ `v_dr` does **not** appear as an active term in any derived operator, dispersion determinant or root.

## ⭐⭐ 1b · COMPLEX FREQUENCY — SUPPLIED. If a root is complex, the branch fixes its continuation

⚠⚠ **This section is supplied, and it is load-bearing.** Let `q_out(ω,k)` be the outward bulk normal
wavenumber — the root of §2's bulk dispersion selected by §1's radiation condition. Its branch points sit
at the **sound cone `ω = ±c_s0|k|`**, **on the real axis**. If B5 yields a complex root, ⛔ **nothing here
says on which side of the real axis it lies; that is B5's to determine, and §0 admits either.** Two
continuation paths that wind differently around a branch point reach **different sheets**, where `q_out`
differs by a factor of `−1`, exchanging a **normal mode** for a **leaky resonance** at the *same* `ω`. ⛔
**The physical requirements below do NOT by themselves fix the continuation.**

**On the REAL axis, `q_out` is fixed by:**
1. `q_out² > 0` (propagating): the bulk solution carries energy **away** from the slab on both sides. ⚠
   Read this as an **energy-flux** condition, valid for **both signs of `ω`**. ⛔ It is **not** a
   phase-velocity condition — that reading breaks for `ω < 0`.
2. `q_out² < 0` (evanescent): the solution **decays** as `|w| → ∞`.

**At COMPLEX `ω`, `q_out` is DEFINED as follows.** ⛔ Not derived, not re-selected — defined:

> `q_out(ω,k)` is the analytic continuation of its real-axis value reached from `ω + i0⁺` by moving
> **downward along the ray of fixed `Re ω`**. Equivalently: deform the inverse-Fourier-in-time contour
> downward from above the real axis while the branch points `ω = ±c_s0|k|` stay fixed and their cuts are
> dragged **vertically downward**, so `q_out` is single-valued on the lower half-plane cut along
> `Re ω = ±c_s0|k|`.

> For `Im ω > 0`, use the analytic continuation **upward from the same upper-rim values `ω + i0⁺`**. No
> dragged cut enters the upper half-plane. This supplies `q_out` there without re-imposing either
> real-axis requirement.

⛔⛔ **REQUIREMENTS 1–2 MUST NOT BE RE-IMPOSED AT COMPLEX `ω`.** Whatever `|w| → ∞` behaviour the
continuation produces there is a **RESULT to report**, ⛔ never a criterion for re-selecting the root.
⚠⚠ **THIS IS THE DIAGNOSTIC §0 REFERS TO:** an engine that re-applies "must decay" at a complex pole lands
on the opposite sheet and can turn a **damped resonance into an apparent instability**; the mirror wrong
re-selection can turn a growing object into apparent decay. ⇒ ⭐ **If you report a growing or decaying
root, state explicitly that you did not re-impose 1–2 to obtain it.**

⭐ **Verify, and report:** that this definition reproduces requirements 1–2 on the real axis. ⛔ If it does
not, report the disagreement rather than adjusting either side.
⭐ **Report the degenerate point** `ω = ±c_s0|k|`: there `q_out = 0`, the two bulk solutions **coalesce**
(the second going linear in `w`), and ⛔ neither requirement selects anything — continuity supplies it.
⭐ **If a root's trajectory crosses `Re ω = ±c_s0|k|` under parameter variation, report that it has LEFT
this sheet** — ⛔ do not re-select it onto one — and say whether the object is a normal mode or a resonance.
⭐⭐ **MAKE ANY DEPENDENCE MEASURABLE:** if B5 yields a complex root, report its imaginary part **also on the
opposite sheet**, and report the ratio where it is defined; if there is no complex root or the ratio is
undefined, report that instead. ⛔ **No expectation is supplied for the ratio.**
**Wrong derivations caught:** leaving the upper-half-plane sheet undefined; re-selecting it by spatial
decay; mishandling the coalescence point; hiding a complex root's dependence on the opposite sheet.
⇒ `S11B_BRANCH_REALAXIS_CHECK`, `S11B_BRANCH_DEGENERATE_POINT`, `S11B_BRANCH_SENSITIVITY`,
  `S11B_SHEET_OF_EACH_ROOT`

## 2 · The bulk sector — supplied acoustics

A scalar superfluid linearised to acoustics, with **no shear modulus**:

```
v = ∇₄ φ ,     δp = −ρ_m ∂_t φ ,     ∂_t² φ = c_s0² ∇₄² φ
```

`ρ_m` is the bulk mass density, `c_s0` the bulk sound speed. Both half-spaces `|w| > W₀/2` are bulk.
⛔ **The BULK's response depends on the slab only through the motion of its faces** — nothing about the
slab's interior enters §2. ⚠ The slab's interior — its stored and kinetic energy — **is** supplied for the
assembly in §5; ⛔ do not read this as "the interior is not needed" for the step as a whole.

⛔⛔ **DO NOT type the DERIVED bulk response — the impedance `Z`, its three regime forms, or the per-face
inertial loading — as a supplied result.** They are **computed objects** (Task **B0b**). ⭐ The bulk normal
wavenumber `q_out` and its branch structure (§1b) follow directly from §2's acoustics under §1's radiation
condition and are used by §1b and B0b; ⛔ what is withheld is the **response** `Z` and its regime
**classification**, ⛔ not the existence or branch points of `q_out`. ⚠ In the two historical stages `Z` and
its regimes were derived in stage A and *imported* into stage B as "established input"; in this unified step
there is no such import — ⭐ **they are derived here**.

## 3 · The window and the projection current

`Ω(w)` is a smooth **window function**, `≈ 1` inside the slab and `→ 0` outside, used to project
four-dimensional equations onto a three-dimensional description. Unless a task says otherwise take `Ω(w)`
**even** in `w`.

The four-dimensional **mass-conservation law** the projection integrates is `∂_t ρ_4D + ∇₄·(ρ_4D v) = 0`,
with `ρ_4D` the 4D mass density of §1 (`ρ_4D = ρ_4D⁰(1 + θ)` in the slab, the bulk density in the
half-spaces) and `v = ∇₄φ` the §2 velocity. `j^w` is the `w`-component of the four-dimensional **MASS**
current `j = ρ_4D v` (⛔ not a number current — the two differ by a conversion factor that would change
reported dimensions).

⭐ **Relative flux.** Where material crosses a **moving** face, the physically meaningful flux is measured
**relative to that face and along its outward normal**:

```
J_± ≡ ρ_m ( v_w − ∂_t ζ_± ) · (±1)        evaluated at the upper (+) / lower (−) face
```

⚠ **Note the plain minus.** The outward-normal sign appears **once**, in the trailing `(±1)`; applying it
again to the face velocity would give a face co-moving with the bulk a nonzero relative flux.
⚠ Report explicitly which signed combination of `J₊` and `J₋` is **net accretion by the slab** and which is
**through-flow**, under the single tag `S11B_FLUX_CHANNELS` (⭐ one named object — ⛔ do not let the two
engines choose two names). ⛔ Do not use a sum of global-`w` fluxes in place of relative ones.

## ⭐⭐ 4 · THE FACE CLOSURE — SUPPLIED, and the affinity is the canonical driver

⭐ **Signs fixed here, because two engines can otherwise obtain opposite damping:** `J_±` is the **mass flux
per unit face area leaving the slab through that face**, measured along that face's **outward** normal
(`J_± > 0` removes material from the slab). `δp_±` is the **bulk** pressure perturbation at that face.

```
J_± = Λ_A(ω) 𝒜_±  +  Λ_V(ω) V_± ,
Λ_I(ω) = Λ_I⁰ / (1 − iωτ_I) ,               I ∈ {A, V, X}
```

with `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰`, `τ_A`, `τ_V`, `τ_X` real and free, each `τ_I ≥ 0`. ⭐⭐ **Keep the three
relaxation times INDEPENDENT throughout.** The affinity and its normalization are **SUPPLIED**:

```
𝒜_± ≡ μ_s − δp_±/ρ_m ,
μ_s ≡ μ_θ/ρ_br⁰ ,       μ_θ ≡ (δU/δθ)|_{u, e_W and all other fields fixed} .
```

⛔ **Use this affinity exactly; do not derive, re-normalize, or re-select it. `μ_θ` and `μ_s` are not
interchangeable.** Derive the dimension of every term (Task **B7**) before combining them.

⭐⭐ **THE AFFINITY `𝒜` IS THE CANONICAL DRIVER OF THE FLUX**, and it carries a **slab-side** contribution
(`μ_s`) that a raw bulk-pressure closure does not. ⚠ **Setting the slab-side contribution to zero
(`μ_s = 0`) turns the affinity-driven part of the flux into a raw-pressure-driven part** — ⛔ the velocity
channel `Λ_V(ω) V_±` is **unchanged**, ⛔ **not** removed. ⚠ That `μ_s = 0` case (a raw-pressure-plus-velocity
closure) is a special case, ⛔ **not** the general law and ⛔ must not be taken as the default; its equal-time
reduction is an **acceptance obligation** (Task **B2c**), ⛔ not a premise.

⭐ **`Λ_X⁰`, `τ_X` are SUPPLIED constitutive inputs** — a **prescribed** reciprocal mechanical response
(§6), ⛔ not derived into existence and ⛔ not inherited from any pressure-only closure. ⭐ This step DERIVES
their **consequences** (the power form, the admissibility region, and **whatever relation between `Λ_X` and
`Λ_V` is forced — if any**), ⛔ never the channel's existence. ⛔ **Assert no such relation in advance.**

⭐ **`V_±` is the OUTWARD face velocity of §1**, ⛔ not the global `∂_t ζ_±`: `J_±` is face-odd by
construction, so pairing it with a face-even quantity through a single scalar coefficient would give
reflection-related faces opposite constitutive terms and spuriously mix `δW` into `ζ_c`.
⭐ **`Λ_I⁰` and `τ_I` are REAL** — the frequency dependence above is the **only** source of complexity in
the closure. ⛔ Do not carry them as generic complex quantities; a dissipation classification asks which
terms are real and in phase with velocity, and that question is vacuous if the constants may themselves be
complex.
⛔⛔ **DO NOT set any `τ_I = 0`, and do not treat an instantaneous law as the default.** ⚠ `τ_I → 0` is the
**memoryless limit** and asserts that conversion across the interface is *instantaneous* — close to the
quantity this programme exists to determine; imposing it would answer a rate question by assumption.
⭐ **Report the `τ_I → 0` limit as a special case**, ⛔ not as the premise, and report behaviour across the
whole range of each `ωτ_I` — small, order unity, large.

## 5 · The brane's stored energy — construct the basis, ⛔ do not assume it closed

Per unit `x`-3-volume, with `e_W ≡ δW/W₀`, the terms **carried** are

```
U = ½ μ_R |∇×u|²  +  ½ B_ρ⁽³⁾ θ²  +  C W₀ θ e_W  +  ½ k_W W₀² e_W²  +  ½ κ_W W₀² |∇(δW)|²
```

with kinetic energy `½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)²`.

⭐⭐ **`C` is the symmetry-allowed CROSS term between densification and thickening, included deliberately.**
⚠ A diagonal energy would already impose that the two channels are separable — which is part of what B2c
determines, ⛔ not an input. **Report how every result depends on `C`, and what changes when `C = 0`.** ⛔ Do
not set it to zero by default. ⚠ `κ_W` is included because "restoring stiffness" alone is ambiguous — a
thickness stiffness may act on `δW` or on `∇δW`, with different `k`-dependence; **report how each result
depends on `κ_W`, and what changes if it vanishes.** ⛔ Do not set it to zero by default.

⛔⛔ **THE LIST ABOVE IS THE SET OF TERMS CARRIED. IT IS NOT ASSERTED TO BE A CLOSED BASIS — CONSTRUCT THE
BASIS YOURSELF AND CHECK IT** (Task **B0-energy**). ⚠ An omitted allowed term can change the dispersion
while a control that removes a *listed* term still reports cleanly. ⭐ **Before deriving anything:**
enumerate the fields and their first gradients (`u`, `∇u`, `θ`, `∇θ`, `e_W`, `∇e_W`) and construct **every**
scalar quadratic in them allowed by the symmetry group below; compare that basis against the list above,
and **report the count you obtain (tag `S11B_ENERGY_BASIS_COUNT`) and any omission you find** — ⛔ the count
is a result, not given here.

⛔⛔ **THE SYMMETRY GROUP, STATED IN FULL.** Isotropy and `w → −w` alone are **not enough**: under those
alone a **pinning term `½K|u|²` is allowed**, and an engine that includes it **gaps the modes** while one
that does not finds them gapless — both obeying the text. The group is:
- **In-plane translation invariance** ⇒ `u` may enter **only through its gradients**, never undifferentiated.
- **In-plane isotropy _and_ parity** — the full `O(3)` acting on the three in-plane directions.
- **Reflection `w → −w`.**
- **Equivalence modulo total divergences** — two densities differing by a total in-plane divergence are the
  same term; ⛔ do not count both.
- ⛔ **Time-reversal is NOT assumed**, and ⛔ no positivity or boundedness is assumed (§0).

⭐⭐ **JUDGE INDEPENDENCE AS FIELD BILINEARS, with B1's constraint NOT applied.** ⚠ Redundancy *"modulo B1's
constraint"* has no single meaning — B1 is sourced, carries memory, and changes rank at `ω = 0`. ⇒ **Carry
EVERY independent invariant with a free symbolic coefficient** and report its effect on B2c/B5. ⭐ **Separately**,
report which basis elements become redundant **once the constraint is applied**, and ⚠ **whether that set
differs between the impermeable and the flux-on case** — ⛔ report the difference rather than choosing one.
⚠ In particular, if the basis construction yields **one or more independent invariants containing `∇·u`**,
⭐ carry **every one** with its own free coefficient — ⛔ do not drop any — and report how each changes B2c's
identification. Which reading your construction produces is itself a result to report.
⇒ `S11B_ENERGY_BASIS`, `S11B_ENERGY_BASIS_COUNT`, `S11B_ENERGY_BASIS_OMISSIONS`,
  `S11B_ENERGY_BASIS_INDEPENDENT_TERMS`, `S11B_BASIS_REDUNDANCY_UNDER_CONSTRAINT`

## ⭐⭐ 6 · HOW TO DERIVE THE EQUATIONS — **BALANCE LAWS**. ⛔ NOT an action principle

⚠ **Under ORDINARY single-copy variation a retarded kernel cannot be varied.** In a bilinear
`∫ λ 𝒴[∂_t δW] dt`, variation transposes the operator, whose symbol is `Y(−ω)` — the **advanced** kernel,
moving its response pole from the lower to the upper half-plane. ⇒ **an irreversible flux may not be placed
inside the varied functional.** ⚠⚠ This is a statement about **ordinary** variation only — ⛔ doubled-variable
(in-in / Galley) constructions **do** yield genuinely retarded kernels, and ⭐ using one as an
**independent cross-check** on the route below is allowed and encouraged.

⛔⛔ **THREE FORBIDDEN ROUTES**, each wrong under ordinary variation, ⛔ not merely disfavoured:
- **(i)** Substituting the mass balance **with its flux source** into `U` and then varying. With
  `(Λ_A⁰, Λ_V⁰) ≠ (0,0)` the relation carries a history functional and may not enter a stored energy.
  ⚠ Signature: an extra root raising the determinant's degree by one — an invented mode.
- **(ii)** Varying `J_±` or `δp_±` inside the action. Manufactures an anti-causal, **energy-generating**
  kernel. ⚠ Signature: a finite pole inherited from a response kernel in the upper half-plane.
- **(iii)** Any route in which a response kernel is differentiated with respect to a field.

### ⭐⭐ THE VIRTUAL-DISPLACEMENT RULE — binding

Use `δ_v` for a virtual variation, to distinguish it from the perturbation `δW`. Let `X` label an in-plane
material point; define the current in-plane map, its Jacobian, the Eulerian slab density, and the
material-area mass:

```
x(X,t) = X + u(X,t) ,                         𝒥_x ≡ det(∂x/∂X) ,
u(X,t) = u(x,t) + O(2) ,
Σ_E(x,t) ≡ Σ(x,t) ≡ ρ_4D(x,t) W(x,t) ,
Σ_mat(X,t) ≡ Σ_E(x(X,t),t) 𝒥_x(X,t) .
```

**A virtual variation is taken at one instant and transfers no mass through either face.** Its binding
constraint is the **material** (⛔ not Eulerian) equation `δ_v Σ_mat(X,t) = 0`. On the uniform background,

```
Σ_E = ρ_br⁰(1 + θ + e_W) + O(2) ,            𝒥_x = 1 + ∇_x·u + O(2) ,
δ_v e_W = δ_v(δW)/W₀ ,
δ_v θ + δ_v e_W + ∇_x·δ_v u = 0 .            (VIRTUAL CONSTRAINT)
```

⛔⛔ **The term `∇_x·δ_v u` MUST NOT BE ABSENT.** The equation `δ_vΣ_E = 0` is **not** the material
constraint and must not be used. Thus `δ_vθ`, `δ_v(δW)`, `δ_vu` are not independent.
⛔ **Do NOT vary `U` with `θ` held fixed.** Impose the VIRTUAL CONSTRAINT either by eliminating one virtual
variation or by a Lagrange multiplier you then eliminate. ⭐ **Report which you used and the multiplier's
physical identity.** ⚠ The **same** multiplier supplies the in-plane restoring force and the thickness
term; ⛔ you may not keep one and drop the other.

⭐ **USE THIS AND NO OTHER:**
1. **State variables `u`, `δW`, `θ`**, with `δ_vθ + δ_ve_W + ∇_x·δ_vu = 0` binding their **variations**.
2. **Constitutive quantities are VARIATIONAL derivatives of `U`:** `μ_θ ≡ δU/δθ` and `p_W ≡ δU/δe_W`.
   ⛔⛔ **Functional, not ordinary partial** — if §5's basis carries gradient terms such as `|∇θ|²`, an
   ordinary partial drops their contribution while the variational derivative keeps a `k²` term.
   ⭐ **State, for every derivative you write, what is held fixed.**
3. **In-plane momentum balance** for `u`.
4. **Thickness equation**, from the variation under the rule above — ⛔ **not** from a verbal description of
   which forces "push the faces apart".
5. **The mass balance WITH its `J_±` source is a separate EVOLUTION equation**, restored **after** the
   variation. ⛔ It is not the same object as `δ_vΣ_mat = 0`.
6. **`δp_±` and the affinity traction proportional to `Λ_X𝒜_±` are PRESCRIBED mechanical responses**
   entering through the external virtual work below. **`J_±` enters the mechanics through the separate mass
   evolution and through its contribution to the prescribed `δp_±`; it has no direct mechanical
   generalized-force term:** `Q_J^direct = 0`. ⛔ Substitute the prescribed responses only after 3–5 are
   complete; ⛔ vary or differentiate them with respect to nothing; ⛔ add no term proportional to
   `J_± δ_v(δW)` or `J_± ∇·δ_v u`, and no mechanical face term beyond the two supplied below.

### ⛔⛔ GEOMETRY AND SIGN CONVENTIONS — a sign error here is an ENERGY-SOURCE error

```
ζ_± = ±δW/2 ,               n̂_± = ±ŵ ,
V_± ≡ (∂_tζ_±)(±1) = ½ ∂_t(δW) ,
δ_vx_± = n̂_± δ_v(δW)/2 ,    t_± = −(δp_± + Λ_X(ω)𝒜_±) n̂_± ,
δ_v𝒲_bulk ≡ Σ_{s=±} t_s·δ_vx_s
            = −½ [δp_+ + δp_− + Λ_X(ω)(𝒜_+ + 𝒜_−)] δ_v(δW) .
```

⭐ Both faces move outward under positive `δ_v(δW)`. ⛔ **Use the displayed virtual-work equation; do not
guess a pressure force and do not append another flux-force term.** ⚠⚠ Taking this sign the other way
reverses the pressure-work term and can manufacture an instability; the mirror bookkeeping error — applying
the reversal to a source-valued traction contribution and treating it as a sink — can manufacture spurious
decay.

### ⭐⭐ THE CAUSALITY DIAGNOSTIC — two bounded checks, ⛔ neither classifies a mode root

In the **equations of motion, closure, face response and dispersion determinant**, every response kernel
must retain the **retarded analytic orientation** supplied here. Run both checks.

1. **Kernel-orientation identity — decisive for every active memory channel.** Before eliminating any face
   quantity, extract, from the equations actually used,

   ```
   K_A,s ≡ coeff_{𝒜_s}(J_s)|_{V_s=0} ,     K_V,s ≡ coeff_{V_s}(J_s)|_{𝒜_s=0} ,
   K_X,s ≡ coeff_{𝒜_s}(−n̂_s·t_s − δp_s) .
   ```

   This is coefficient extraction, ⛔ not variation of a response. For each `I ∈ {A,V,X}`, put
   `K_I,s(ω) − Λ_I(ω)` over a common denominator, cancel its polynomial gcd, and report whether its
   numerator is identically zero as a rational function of symbolic `ω`, **before** specializing any
   coefficient or time. Then carry independent inert placeholders `ℓ_A, ℓ_V, ℓ_X` from these checked
   coefficients through the face elimination and balance-law substitution, replace them only at the end by
   `Λ_A(ω), Λ_V(ω), Λ_X(ω)`, and report the symbolic differences from the face response and equations
   actually used; recompute the dispersion determinant from those checked equations and report its
   difference from the reported determinant under the same row normalization. ⛔ Do not conjugate a
   placeholder or send `ω → −ω` during this propagation. Judge **rational-function identity**, ⛔ not literal
   symbol appearance; a check at the single point `ω = 0` is invalid. **Coverage boundary:** for either sign
   of any nonzero `Λ_I⁰`, this detects a transposition whenever `τ_I > 0`; when `τ_I = 0` retarded and
   transposed kernels coincide, and when `Λ_I⁰ = 0` the channel is absent — report those as
   **indistinguishable**, ⛔ not as checked passes.
2. **Reduced-object pole inventory — partial coverage only.** Cancel removable factors in each downstream
   object, track every finite bare-kernel pole, and report whether it is retained, cancelled or
   feedback-displaced, with the location of every retained or displaced pole. A retained bare-kernel pole
   above the real axis is an advanced response. ⛔ A feedback-displaced pole is **not** an orientation
   verdict: with the declared free coefficients its half-plane is sign-indeterminate for either orientation.
   Conjugated kernels in the power forms of §4/§6 and all zeros of the dispersion determinant are outside
   this inventory.

⚠ **If either check fails, re-check rule step 6 and report both outcomes.** ⛔ Do not use either diagnostic
to delete, suppress or reclassify any growing or decaying dispersion root.
⇒ `S11B_CAUSALITY_CHECK`, `S11B_KERNEL_ORIENTATION_IDENTITIES`, `S11B_KERNEL_PROPAGATION_RESIDUALS`,
  `S11B_KERNEL_POLE_LOCATIONS`

### ⭐⭐ TWO MANDATORY CONVENTION CROSS-CHECKS — ⛔ run them, ⛔ they must be able to fail

⚠ The variational convention has **two candidate readings that give different dispersion relations**; one
wrong reading can manufacture growth, its mirror can conceal growth or manufacture apparent decay. ⛔ Report
both outcomes explicitly.

**(a)** The in-plane equation your variation produces must carry the restoring force `−∇(δU/δθ)`. ⭐ This
single check selects the convention uniquely. ⛔ If your in-plane equation lacks it, your variational rule
is wrong — ⛔ fix the rule, do not patch the equation. **Wrong derivation caught:** omitting `∇·δ_vu` from
the VIRTUAL CONSTRAINT, or varying at fixed `θ`, removes this contribution.

**(b)** With NO bulk, `Λ_A⁰ = Λ_V⁰ = Λ_X⁰ = 0`, `κ_W = 0`, `k = 0`: report `ω²` for the thickness mode and
**the inequality on `B_ρ⁽³⁾`, `C`, `k_W` under which it is positive.** ⭐ It must be positive for **every**
`U` that is positive-definite — ⛔ not merely on a convenient sub-region. ⭐ **State explicitly whether
`B_ρ⁽³⁾` appears in that `ω²` at all.** The kinematic reduction is

```
J_+ = J_- = 0 , k = 0  ⇒  ∂_t(θ + e_W) = 0  ⇒  θ + e_W = constant ,
K_check ≡ d² U(θ = constant − e_W, e_W)/de_W² .
```

⭐ Report whether the thickness stiffness from the equation equals `K_check`. ⛔ **Do not print `K_check`'s
value or the stability boundary here — it is a computed object** (⇒ §12, withheld). **Wrong derivation
caught:** a fixed-`θ` thickness variation, or omitting the common constraint multiplier, gives a stiffness
inconsistent with the independently reduced stored energy.

⛔⛔ **SCOPE OF CHECK (b):** confined to **no bulk, no permeation, no reciprocal traction, positive-definite
`U`**, where the system is conservative and growth or decay would be a **derivation error**. ⛔ It says
**nothing** about the full problem and ⛔ **must not** be used to reject a growing or decaying root of B5.
⇒ `S11B_CONVENTION_CHECK_INPLANE`, `S11B_CONVENTION_CHECK_CONSERVATIVE`,
  `S11B_CONSERVATIVE_POSITIVITY_INEQUALITY`

### ⭐ ENERGY ACCOUNTING — three discriminators that can fail

⚠ Enumerating loss channels **from the same equations used to derive them** and checking power balance
against that enumeration is an identity for any system built that way — **it cannot fail.** Use these
independent discriminators.

1. ⭐ **Compute `d/dt(T + U)` from YOUR equations.** Report every signed external exchange term — sink or
   source — and name the transport process each corresponds to. ⛔⛔ **If a term corresponds to no transport
   process, REPORT IT** as a defect. **Wrong derivation caught:** a response-kernel or flux-force term not
   among the supplied pressure, reciprocal-traction and transfer terms creates an exchange with no
   transport counterpart.
2. ⭐⭐ **INDEPENDENT PRESSURE-WORK SIGN CHECK.** Restrict **only this diagnostic** to real `ω`, propagating
   `q_out² > 0`, impermeable faces, and the form cut `Λ_X⁰ = 0`. In this sub-case B0b's impermeable
   kinematics give `J_± = 0`, `v_bulk,± = V_±`, `δp_± = Z V_±` (⭐ with `Z` **B0b's derived** object, ⛔ not
   supplied here), `P̄_bulk,± = ½ Re(δp_± V_±*)`, and `d(T+U)/dt|_pressure = −Σ_{s} P̄_bulk,s`. Compute the
   left side **off shell** by pairing the in-plane momentum equation with `(∂_t u)*`, the thickness equation
   with `(∂_tδW)*` (⛔ not with `V_s*`), and the unnormalized linearized mass balance with `μ_s` in the
   `½ Re(μ_s J_s*)` convention, **leaving the harmonic amplitude free**. ⛔ Do not impose the homogeneous
   thickness equation or dispersion relation, and ⛔ **do not first replace the left side by the literal
   period average of an exact total derivative**. Compute the **right** side independently **from the
   outgoing bulk flux**, and report the symbolic difference without changing either derivation. **Wrong derivation caught:** reversing the traction or
   face-displacement sign turns slab pressure work into a source while the independent bulk power keeps its
   sign. ⛔⛔ **SCOPE:** this checks **only** the pressure-work sign in that sub-case; it classifies no root.
3. ⭐⭐ **FULL TWO-PORT BALANCE CHECK.** At real `ω`, keep `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰` and all three `τ_I`
   symbolic. **Off-shell coefficient check:** pair each equation as in check 2, keep the harmonic amplitudes
   and the face quantities `δp_s`, `J_s`, `𝒜_s` algebraically free, ⛔ do not impose the homogeneous
   equations, determinant or **any on-shell amplitude relation**, and ⛔ **do not replace the result by the
   literal period average of an exact total derivative**; decompose the slab power by order in `Λ_X⁰`, and at
   every order compare it face-by-face and channel-by-channel with the independently supplied slab-side
   exchange
   `−½ Σ_{s} Re[ (δp_s + Λ_X(ω)𝒜_s) V_s* + μ_s J_s* ]`. Report every symbolic difference without changing
   either derivation. ⚠ **Verified blind spot:** the face response, closure and `δp/ρ_m` normalization
   cancel from this comparison — it tests the slab traction's **sign, factor and face count**, the
   constraint/multiplier bookkeeping and the slab-side `μ_s` conjugate normalization, ⛔ not the traction's
   analytic structure; a transposed `Λ_X` kernel used identically on both sides is invisible here. This
   check classifies no root.
⇒ `S11B_ENERGY_SINKS`, `S11B_ENERGY_SOURCES`, `S11B_UNATTRIBUTED_SINK_TERMS`,
  `S11B_UNATTRIBUTED_EXCHANGE_TERMS`, `S11B_PRESSURE_WORK_SIGN_CHECK`, `S11B_FULL_TWO_PORT_BALANCE_CHECK`

⛔⛔ **DO NOT USE THE ZERO-INTERFACE-COEFFICIENT LIMIT AS A CHECK ON THE ROUTE.** All candidate routes
**coincide** at `Λ_A⁰ = Λ_V⁰ = Λ_X⁰ = 0`, so that limit cannot discriminate a correct derivation from a
forbidden one. ⚠ It is also **not** the lossless limit — a **bulk radiation** loss (energy carried to
infinity by §2's bulk, present even with impermeable faces) survives it; ⛔ do not describe it as
"reversible".

---

## 7 · ⭐⭐⭐ THE STRUCTURAL RULE — verbatim, non-negotiable

> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION/ENERGY, the
> ANSATZ, and the SUPPLIED constitutive laws and routes above. Every other expression involving them must
> be REACHED BY COMPUTATION. Every control re-enters the chain at the ENERGY or a SUPPLIED law, ⛔ never at
> a result.**

⛔⛔ **A hand-typed CAS object is still hand-typed.** `emit[Simplify[(alpha + beta)^3]]` is a genuine CAS
object with **no data dependency on the derivation**: delete the solve and the output does not move. ⭐ **The
test: if you deleted the computation, would this tag change?** If not, it is dead. ⚠ Swapping the symbols of
a real answer is **not** enough — a list of the right length with a zero in the right slot states the result
whatever the other entries are called.

## 8 · ⭐⭐⭐ THE THREE CLAUSES, THE COROLLARIES, AND NO VERDICT

> **1. A script may PRINT computed objects. It may NOT state conclusions.** Every `emit`/`Print` payload is
> a CAS object — an expression, a solved root, a matrix, an integer from a rank computation, a boolean from
> a symbolic test. ⛔ Never author prose describing a result.
>
> **2. EMIT BOTH OPERANDS AND THE RESIDUAL, then guard.** ⛔ A residual **asserted** zero always prints `0`
> and carries no information — asserting `residual == 0` writes down the expected output. Emit operand A,
> operand B, **and** `A − B`; emit **first**, and only then — ⛔ **on an OPERATIONAL failure only** — may you
> hard-stop. ⛔ A physics disagreement (a nonzero residual) **EXITS 0** (see NO VERDICT).
>
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

**Corollaries:**
1. ⛔ **A hand-typed CAS object is still hand-typed** — see §7.
2. ⛔ **The tag NAME is output too.** A name may name **the object**; ⛔ never its value, sign, count, ratio,
   parity, region, or the shape of the answer. Schematically — **placeholder** content, ⛔ not this step's:
   `..._IS_POSITIVE`, `..._VANISHES`, `..._EQUALS_ZPERM`, `..._STABLE` are forbidden **name** patterns
   whatever the payload. ⭐ `..._RESPONSE`, `..._DISPERSION`, `..._ROOTS`, `..._RESIDUAL`, `..._COUPLING`
   are fine.
3. ⛔ **No tautological residual.** Before emitting a difference, ask whether the two operands came from
   **independent routes**. ⭐⭐ **The test:** mutate the **source** the residual is supposed to police,
   recompute downstream, and see whether the residual moves. A residual that moves only when you corrupt an
   **already-computed operand** tests your arithmetic, ⛔ not the object. ⭐ Where no second route exists,
   emit the objects and ⛔ emit no residual. ⭐ Verify independence by **one-sided corruption**: break route
   A only; if route B's object moves too, they were never independent.
4. ⛔⛔ **EMISSION MUST NEVER BE CONDITIONAL ON A PAYLOAD'S VALUE.** Whether a tag appears may depend only on
   **which named object and which control package** it belongs to. ⚠ A value **present and identical across
   controls is a RESULT**; a value **absent** is indistinguishable from *never computed*. ⛔ Never
   de-duplicate payloads across controls. Tag **names** are unique; **payloads may legitimately repeat.**

### ⛔⛔ NO VERDICT

⛔ **There is no `VERDICT` tag, no `PASS`, no `FAIL`, no summary judgement anywhere in either engine.** ⛔ **A
physics finding must EXIT 0.** Exit non-zero **only** on operational failure — an exception, an unsolvable
system. ⚠ Otherwise a builder iterating to exit 0 can make a genuine disagreement disappear, and the
disagreement is the most valuable output available. ⚠⚠ **This is a change from the two historical S11b
stages, which ended in a `VERDICT` line.** A boolean-valued test is emitted **as the CAS object the test
returned**, with its operands beside it, ⛔ never serialised as a host-language native boolean.

### ⛔⛔ THE LOCUS PROTOCOL — used wherever a task asks for a locus, and ⛔ never varied

⚠⚠ **"Solve this over the reals" is NOT a specification.** On the CAS both engines use, a multivariate solve
**ignores** real declarations and returns solutions in the algebraic closure; the single-variable real
solver **rejects** a tuple; and the ordinary solver returns the **same empty token** for a system that is
identically satisfied and one that is inconsistent. ⇒ ⛔ "solve over the reals" makes *checked-and-empty*
indistinguishable from *true-everywhere*.

⭐⭐ **Therefore, wherever this file asks for a genuine equation-system locus** — the **degenerate loci of
B0c** (coefficient values at which the face closure loses solvability) — **emit ALL of these and nothing
less:**

| suffix | the object |
|---|---|
| `_EQUATIONS` | the defining equation system, as CAS relations, **before** any solve |
| `_SOLUTION` | the solution set as your CAS returns it, solve variables named, domain left unrestricted |
| `_IDENTICALLY_SATISFIED` | the typed symbolic test of whether **every** equation simplifies to zero identically in its variables |
| `_INCONSISTENT` | the typed symbolic test of whether the system admits **no** solution at all |
| `_REAL_ADMISSIBLE` | for **every** branch in `_SOLUTION`: that branch and the typed symbolic test of whether it admits a point satisfying **every** §1/§4/§5 premise in force |

Each `_IDENTICALLY_SATISFIED` / `_INCONSISTENT` payload has the ordered fields
`STATUS_TOKEN` (one of `PROVED_TRUE · PROVED_FALSE · UNDECIDED`), `TEST_OBJECT` (the live CAS boolean, in
re-parseable form), `OPERANDS` (the live operands). Each `_REAL_ADMISSIBLE` entry, in `_SOLUTION` branch
order, has `BRANCH`, `STATUS_TOKEN` (`ADMISSIBLE · EXCLUDED · UNDECIDED`), `TEST_OBJECT`, `OPERANDS`.
⛔ **A BRANCH IS NEVER SILENTLY DROPPED.** ⭐ `UNDECIDED` is an explicit **coverage finding** — neither
agreement nor disagreement — and forbids the corresponding claim of real emptiness, non-emptiness, branch
admission or exclusion. ⛔ Emit every object whatever it comes out to.

⛔⛔ **THIS PROTOCOL IS FOR EQUATION-SYSTEM LOCI ONLY.** ⚠ B2d's admissibility is a **region / inequality** —
report it **in whatever form it takes** (eigenvalues, principal minors), ⛔ **not** through
`_EQUATIONS`/`_SOLUTION`, which would collapse it to the boundary equality. ⚠ B5's grazing is a
**classification** (bound / threshold / BIC / resonance), ⛔ not a solve. ⛔ Do not force either through this
protocol.

---

## 9 · What to compute — the task list

⛔ **This section says what to produce. It does not say what any of it equals.** Tasks are ordered by
dependency: the bulk face response (B0) is derived first, then the assembly consumes it.

### B0 · The bulk face response — DERIVE it (this is the historical A-material, recomputed here)

**B0a · Projection identity.** Integrate §2's four-dimensional continuity equation against `Ω(w)` over `w`,
integrating by parts to isolate the term carrying `j^w`. Report the source term for a **finite** interval
`[w₁, w₂]` and for an **infinite** one. With `Ω(w)` **even**, evaluate the source for `j^w(w)` (a) even and
(b) odd in `w`; ⭐ **report the interval used and whether it is symmetric about `w = 0`**, and for each
result whether it is **exact and on what interval** — ⛔ an oddness argument does not by itself fix an
integral over an asymmetric interval. Then repeat with a **dynamic** window `Ω = Ω(w; x, t)` and enumerate
**every** term present then and absent from the static case.
⇒ `S11B_PROJECTION_FINITE`, `S11B_PROJECTION_INFINITE`, `S11B_PARITY_EVEN_JW`, `S11B_PARITY_ODD_JW`,
  `S11B_PARITY_INTERVAL`, `S11B_DYNAMIC_WINDOW_EXTRA_TERMS`

**B0b · Bulk response to moving faces — impermeable.** Solve §2 with both faces displaced as `ζ₊`, `ζ₋`,
imposing that the bulk normal velocity at each face equals that face's normal velocity; treat both
half-spaces under §1's radiation condition and §1b's branch object. Report `Z` (as defined in §1) and its
real and imaginary parts, **separately for the `δW` and `ζ_c` combinations**. ⭐⭐ **There are THREE
regimes** — the bulk normal wavenumber squared positive, negative, **or zero**. Report all three, including
the behaviour of every reported quantity **at** the zero (grazing) case, where some may be singular. ⚠
Omitting the third case is a known prior defect. ⭐ **Inertial loading**, defined so both engines report the
same object: in a regime where `Z` is purely imaginary, `m_add` is the coefficient in
`δp|_face = m_add · ∂_t² ζ_±`, reported **PER FACE**. ⛔ Do not report a two-face sum, and ⛔ do not fold in
any half-amplitude convention. ⚠ **Trap (measured):** the per-face inertial loading acts with the **same
sign against each face's OUTWARD acceleration**; the global-`w` signed pair `(X, −X)` is a convention
artifact ⇒ ⛔ do not consume it as an inertia. ⛔ **Report the magnitude as the computed object** — ⛔ do not
type it here.
⇒ `S11B_Z_IMPERMEABLE`, `S11B_Z_BY_REGIME`, `S11B_Z_BY_PARITY`, `S11B_ADDED_MASS`, `S11B_GRAZING_BEHAVIOUR`

**B0c · Permeable face response — from the §4 closure.** Two relations close the face:
- the bulk, `δp = Z · v_bulk,±`, with `v_bulk,±` the **bulk material's** outward normal velocity;
- **interfacial mass balance** (kinematics, ⛔ not a result): `v_bulk,± = V_± + J_±/ρ_m`.

⭐ **Solve these together with the §4 closure**, keeping `τ_A`, `τ_V` independent, and report the face
pressure `δp` as a function of `V_±` and `μ_θ`; ⭐ **report it as what it is** — ⛔ do not force it into the
shape of a single velocity-impedance if the algebra does not give one, and report **whatever coefficients it
actually has**. ⭐ **Report any coefficient loci** where a face equation loses its dependence on the bulk
amplitude (a driven face then having **no** solution, an undriven one a **free** amplitude): solve for them
and report the set your CAS returns — ⭐ the **empty set** if there are none — via the §8 locus protocol.
⛔ Do not select a value for any coefficient, and ⛔ do not report a single "the" permeable response.
⭐ **Permeable-response dissipation audit.** For **each of B0b's three regimes and each parity combination**,
report whether a **dissipative part** (real, in phase with velocity) is present, **on which coefficient it
depends**, and **how it varies with each `ωτ_I`** — both limits included. ⛔ Report what the calculation
gives; ⛔ assert nothing in advance about which regime is dissipative.
⇒ `S11B_FACE_RESPONSE`, `S11B_FACE_RESPONSE_COEFFS`, `S11B_DEGENERATE_LOCI_*` (locus protocol),
  `S11B_PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY`, `S11B_PERMEABLE_DISSIPATION_VS_OMEGA_TAU`

### B1 · The mass-balance constraint

⛔⛔ **The exact balance is GIVEN — kinematics, not a result — and your job is to linearise it.** Let
`Σ(x,t) ≡ ρ_4D(x,t)·W(x,t)` be the slab-integrated mass per unit `x`-3-volume and let the slab material move
with `v = ∂_t u`:

```
∂_t Σ  +  ∇_x·( Σ v )  =  −( J₊ + J₋ )
```

with `∇_x·` the divergence in the **three** in-plane directions and `J_±` from §3. ⭐ **Linearise** to first
order in `θ`, `e_W`, `u` and the flux, and ⭐⭐ **report every term of the linearised constraint together
with the term of the exact balance it came from** — ⛔ nothing dropped without account, ⛔ no term added that
is not a consequence. Then state how many **independent internal degrees of freedom** survive and why; ⚠
**define what you count** (fields, amplitudes, or independent initial data) before and after the constraint,
at fixed `(k, ω)`.
⇒ `S11B_CONSTRAINT`, `S11B_CONSTRAINT_TERM_ORIGINS`, `S11B_INTERNAL_DOF_COUNT`,
  `S11B_DOF_COUNTING_CONVENTION`

### B2 · The equations of motion, the thickness response, the compressional response

**B2a · EOM.** Following §6, derive the in-plane equation and the thickness equation, including the force
the bulk exerts on **both** faces via the B0c face response and the reciprocal traction proportional to
`Λ_X`. Report both operators, isolate the `Λ_X` term in the thickness operator, and report the symbolic
difference from the `Λ_X⁰ = 0` operator.
⇒ `S11B_INPLANE_EOM`, `S11B_THICKNESS_EOM`, `S11B_BULK_FORCE_ON_THICKNESS`,
  `S11B_RECIPROCAL_TRACTION_THICKNESS_EFFECT`

**B2b · Thickness response.** Solve the thickness equation for its response. ⭐ **State explicitly what the
response is a ratio of** (which output field to which driving quantity) and give its dimensions. Report the
bulk's contribution to the thickness operator **in each regime of B0b**, decomposed into a part in phase
with `∂_t²δW` (inertia-like) and a part in phase with `∂_tδW` (velocity-like). ⛔⛔ **Do NOT fold a
velocity-like contribution into an "effective inertia":** a velocity-like bulk load is a resistance, and
calling it a mass smuggles damping into inertia. ⭐ Report which regimes, if any, admit a mass
interpretation at all — ⛔ do not assume in advance which do.
⇒ `S11B_THICKNESS_RESPONSE`, `S11B_RESPONSE_NORMALIZATION`, `S11B_BULK_OPERATOR_BY_REGIME`,
  `S11B_MASS_INTERPRETATION_VALID_WHERE`

**B2c · Compressional response, and the ACCEPTANCE OBLIGATION.** Eliminate the thickness degree of freedom
and report the in-plane compressional stress response. ⭐ **State explicitly what this response is a ratio
of** — which stress component to which measure of deformation — **before** reporting it. ⚠ **Check whether
B1's constraint changes rank at exactly `ω = 0`**; if an integration constant survives, say what fixes it,
and report whether dividing by `ω` before or after taking the limit changes the answer. ⭐ Report its
behaviour as `ω → 0` and `ω → ∞` **along a stated path in `(ω,k)`** — ⛔ the limits need not commute; name
the path and report whether another path differs. Then report **where a modulus measured with the thickness
held fixed would sit**, or that no consistent identification exists.

⭐⭐ **THE ACCEPTANCE OBLIGATION (a reduction the orchestrator diffs — ⛔ the target is NOT printed here).**
Form the **equal-time, slab-side-off slice** of your B0c face response: set the **slab-side** contribution
to `𝒜_±` to zero (`μ_s = 0`, so the affinity-driven flux becomes raw-pressure-driven with coefficient
`Λ_p⁰`; ⛔ the velocity channel `Λ_V(ω)V_±` is **unchanged**) and specialize `τ_A = τ_V = τ`. ⭐ **First
compute the coefficient mapping the supplied affinity forces** — `S11B_ZPERM_SLICE_MAP`, the relation
between `Λ_p⁰` and `Λ_A⁰` implied by the `μ_s = 0` reduction of §4 (⭐ a value **computed from the supplied
affinity**, ⛔ not independently supplied and ⛔ not withheld) — **before** solving. Then emit the specialized
face response as `S11B_ZPERM_SLICE`. ⛔⛔ **DO NOT print any expected form for `S11B_ZPERM_SLICE`, and ⛔ do
not adjust your derivation to reach one** — an independently-derived reference response is held by the
orchestrator and the diff happens on our side; a disagreement is a **finding**, ⛔ not a build failure.
⛔ Emit **no in-engine residual** against a reconstructed target (that would invite the "iterate until it
matches" loop); emit only `S11B_ZPERM_SLICE_MAP` and `S11B_ZPERM_SLICE` for the orchestrator-side diff. ⚠
This reduction checks the algebra combining the bulk relation, interfacial mass balance and closure; it
does **not** validate the affinity's sign/normalization, which is fixed by the supplied affinity. ⚠ It
applies **only** after the stated specialization; the **general independent-time face response** is emitted
under **B0c** and is **unchecked by this obligation** — state that plainly.
⇒ `S11B_COMPRESSIONAL_RESPONSE`, `S11B_LIMITS_AND_PATH`, `S11B_FROZEN_THICKNESS_IDENTIFICATION`,
  `S11B_ZPERM_SLICE_MAP`, `S11B_ZPERM_SLICE`

**B2d · Passivity, reciprocity, admissibility.** ⛔⛔ **THE TRANSFERRED-MASS PRESSURE WORK MUST NOT BE
ABSENT.** With `P̄ = ½ Re(force × velocity*)` and positive `P̄_out,±` energy leaving the slab, the complete
two-port identity is

```
P̄_out,±  = ½ Re[ (δp_± + Λ_X 𝒜_±) V_±* + μ_s J_±* ]
         = ½ Re[ δp_± v_bulk,±* + 𝒜_± J_±* + Λ_X 𝒜_± V_±* ] ,     v_bulk,± = V_± + J_±/ρ_m .
```

⛔ The mixed expression `½ Re(δp_±V_±* + 𝒜_±J_±*)` is incomplete by `½ Re(δp_±J_±*/ρ_m + Λ_X𝒜_±V_±*)`; the
`Λ_X⁰ = 0` comparison standard is obtained **solely** by making that substitution, every other coefficient
symbolic.

1. **Port dissipativity.** On real `ω`, with input vector `a_± ≡ (V_±, μ_s)ᵀ` and
   `P̄_out,± ≡ ½ a_±† H_±(ω,k) a_±`, `H_± = H_±†`: construct `H_±` including `Λ_X`, and compute the
   necessary-and-sufficient condition `H_± ⪰ 0` for every `a_±` (eigenvalues or principal minors). Report
   whether that condition depends on parameters alone or also on `(ω,k)` through `Z`. ⛔ Do not emit
   `NOT_ESTABLISHED`.
2. **Thermodynamic admissibility and reciprocity.** With `f_X,± ≡ Λ_X(ω)𝒜_±` and interfacial power
   `½ Re(𝒜_±J_±* + f_X,±V_±*)`, the supplied mixed representation is

   ```
   ( J_±   )   ( Λ_A  Λ_V ) ( 𝒜_± )
   ( f_X,± ) = ( Λ_X   0  ) ( V_±  ) .
   ```

   Construct its Hermitian power form for arbitrary complex `(𝒜_±, V_±)` and report the
   necessary-and-sufficient condition on `Λ_A⁰, Λ_V⁰, Λ_X⁰, τ_A, τ_V, τ_X` for nonnegative interfacial
   entropy production. Time-reversal invariance remains **not assumed**. For the conditional Onsager–Casimir
   test, `𝒜_±` and `f_X,±` are even and `J_±`, `V_±` odd under time reversal; ⛔ **do not apply
   Onsager–Casimir by transposing the displayed mixed matrix.** Use this **binding** conversion: for a
   general mixed law `(J₁; X₂)ᵀ = ((a,b),(c,d))(X₁; J₂)ᵀ`, with `d ≠ 0`, solve the second row for `J₂` to
   obtain the all-flux law `(J₁; J₂)ᵀ = ((a−bd⁻¹c, bd⁻¹),(−d⁻¹c, d⁻¹))(X₁; X₂)ᵀ`, apply Onsager–Casimir to
   **that** all-flux matrix at the same real `ω` using each state variable's parity (the opposite of its
   displayed rate's), ⛔ inserting no conjugation merely because the Hermitian form uses one. For the
   interfacial law use a formal scalar `d = ε ≠ 0` only for this conversion, clear all powers of `ε`, then
   take `ε → 0`. **Cross-check** by instead solving the first row where `a ≠ 0` (all-force resistance form),
   clearing `a`, extending to `a = 0`; the two conversions must give the **same** relation — report a
   disagreement rather than choosing.

   ⭐ Compute and report the relation between `Λ_X` and `Λ_V` **if one is forced**. Report the conditional
   reciprocity region and the unconditional admissibility region **separately**, then their set relation
   (equality, strict nesting either way, or incomparability). Separately for admissibility and for
   conditional reciprocity, report whether any relation among `τ_A, τ_V, τ_X` is forced. Report the
   necessary-and-sufficient coefficient condition **in whatever form it takes**, both without and (if
   applicable) with reciprocity. ⛔ The supplied laws determine the admissibility form — do not emit
   `NOT_ESTABLISHED` — and ⛔ **assert no outcome in advance**, ⛔ treat no coefficient as the expected one.

⛔⛔ **DO NOT IMPOSE EITHER CONDITION TO REMOVE A GROWING OR DECAYING ROOT.** ⭐ Instead, report whether each
growing or decaying root of B5 lies inside each computed region, with the set relation between the two
regions; admissibility is a **classification, ⛔ not an acceptance gate**.
⇒ `S11B_TWO_PORT_POWER_IDENTITY`, `S11B_PORT_DISSIPATIVITY`, `S11B_PORT_CONDITION_KIND`,
  `S11B_ONSAGER_CONDITION`, `S11B_ONSAGER_RECIPROCITY`, `S11B_ONSAGER_DETERMINABLE`,
  `S11B_RELAXATION_TIME_RELATIONS`, `S11B_COEFFICIENT_ADMISSIBILITY`, `S11B_GROWTH_INSIDE_ADMISSIBLE`,
  `S11B_DECAY_INSIDE_ADMISSIBLE`

⚠ **Two more traps, both measured in the historical stages:** (1) a **real** part of `Z` in the
**propagating** case is **bulk radiation** — energy carried to infinity by §2's bulk, present **with
impermeable faces** — so ⛔ **do not read any `Re Z` as evidence of transfer *through* the interface**;
which real part, if any, the **closure** creates is B0c/B2d's to determine, ⛔ not stated here. (2) The
rest-frame linearisation discards a relative correction that is **B9's to derive** from the convective
operator — ⛔ this file states neither its order nor its value.

### B5 · The longitudinal mode — its FATE, and the grazing threshold (⭐ IN SCOPE here)

Assemble and report the **dispersion relation** of the coupled longitudinal system. Report whether it admits
a closed-form `ω(k)`, whether roots are real, and for any complex root its imaginary part and **which
physical ingredient makes it nonzero** — ⭐ **attributing** the loss either to bulk radiation (present with
impermeable faces) or to a loss the closure creates, per the trap above. ⭐ **If a root fails to exist as a
normal mode, report that.** ⭐⭐ **Report the
SIGN of every imaginary part and classify each root as decaying or GROWING**, and report the condition on
the moduli and on `C` that separates the two. ⛔ A growing root and a decaying root are **both reportable
results** — ⛔ do not suppress either, re-branch to remove it, or assume the quadratic form is
positive-definite.

⭐⭐ **THE GRAZING THRESHOLD IS THIS STEP'S QUESTION** (uniform background). Report the coupled system's
behaviour **at and near `q_out → 0`** — §1b's degenerate point, the sound-cone branch locus of the supplied
bulk operator: whether the assembled longitudinal object there is a bound mode, a threshold state, a bound
state in the continuum, or a resonance, and the mechanism. ⛔ Do not defer this to S11b-C — only the **non-uniform** spectrum is
C's; the uniform grazing behaviour is settled here.

⭐ Carry `Λ_X` symbolically through the determinant and roots; compare with the form cut `Λ_X⁰ = 0` and
report which roots, multiplicities and classifications change. ⛔ Assert no outcome in advance.

⭐⭐ **FOR EVERY GROWING OR DECAYING ROOT, REPORT ALL THREE (the burden is identical for either sign):**
1. **Both §6 causality checks** — every kernel-orientation and propagation residual; separately inventory
   every finite response-related pole after removable factors are cancelled and classify it bare, cancelled
   or feedback-displaced. ⛔ Do not count the power-form conjugations or the dispersion root itself as a
   kernel pole.
2. **The sheet the root sits on** (§1b), and confirmation that requirements 1–2 of §1b were ⛔ **not**
   re-imposed at complex `ω`.
3. ⭐⭐ **Whether it lies INSIDE each admissible region of B2d**, with the set relation between the two
   regions. ⛔ Do not use either region as a gate that removes or re-branches any root.
4. ⭐ **If it lies OUTSIDE the passive region** (a non-passive growing/decaying root): via §6 energy-
   accounting discriminator 1, report the energy-exchange term that powers it and the transport process it
   corresponds to. ⛔ A non-passive root is a **reportable finding conditioned on a named power source**,
   ⛔ never labeled "forbidden" or "dead". ⚠ The background drain `v_dr` is the candidate reservoir, recorded
   as a scope limit (§1, B9) — ⛔ it is not an active operator term here.
⇒ `S11B_LONGITUDINAL_DISPERSION`, `S11B_ROOTS`, `S11B_IMAGINARY_PART`, `S11B_DISSIPATION_ORIGIN`,
  `S11B_ROOT_STABILITY_CLASS`, `S11B_STABILITY_CONDITION`, `S11B_GRAZING_MODE_CLASSIFICATION`,
  `S11B_ROOT_POWER_SOURCE`, `S11B_GROWTH_ARTIFACT_DIAGNOSTICS`, `S11B_DECAY_ARTIFACT_DIAGNOSTICS`,
  `S11B_RECIPROCAL_TRACTION_ROOT_EFFECT`

### B4 · The breathing-mode stability slice — the `k = 0`, impermeable reduction (⭐ IN SCOPE)

⭐ On the slice `k = 0`, **impermeable** faces (`Λ_A⁰ = Λ_V⁰ = 0`) and **no** reciprocal traction
(`Λ_X⁰ = 0`): assemble and report the **thickness-mode dispersion**, ⚠ **keeping the bulk load** — at
`k = 0` the bulk still reacts on the breathing faces; ⛔ do not drop it. Report whether the root is real or
complex, and the condition on the moduli and on `C` that **separates a stable from a growing** breathing
mode. ⛔ Do not assume the stored energy is positive-definite, ⛔ do not print the stability quantity or its
boundary (⇒ §12, withheld — a computed object), and ⛔ do not re-branch to remove a growing root.
⚠ **This slice is NOT the conservative check (b) of §6**, which removes the bulk **entirely**; here the bulk
radiation is retained, so the two are different objects — report both.
⇒ `S11B_BREATHING_SLICE_DISPERSION`, `S11B_BREATHING_STABILITY_CONDITION`

### B6 · The transverse mode — computed, ⛔ not asserted

On this uniform background, compute the coupling between the transverse in-plane mode and the thickness
degree of freedom **from B1's constraint and §5's energy**, ⛔ **not by asserting a divergence-free
argument.** Report the coefficient and any modification to the transverse dispersion. ⭐ **State explicitly
what the coefficient couples to what**, and its normalization, ⛔ **before** assigning it a value or a
dimension; ⚠ **if it vanishes identically, say so, and say that its normalization is then undetermined.** ⭐
Report whether the transverse mode acquires an imaginary part, and its dependence on
`Λ_A⁰, Λ_V⁰, Λ_X⁰, ωτ_A, ωτ_V, ωτ_X` and the slab-side part of `𝒜` across the full range.
⚠ ⛔ **Whatever this returns, it does NOT settle whether confinement is unconditional** — that is a
non-uniform question, S11b-C's. ⛔ Do not phrase it as if it does.
⇒ `S11B_TRANSVERSE_COUPLING`, `S11B_TRANSVERSE_DISPERSION`, `S11B_TRANSVERSE_DISSIPATION`

### B7 · Dimensions and homogeneity

Derive, from the equations above and ⛔ **not** from any table or registry, the `[L,T,M]` exponents of:
`Z`, `m_add`, `ρ_m`, `c_s0`, `v_dr`, `ρ_br⁰`, `B_ρ`, `B_ρ⁽³⁾`, `μ_W`, `k_W`, `κ_W`, `C`, `μ_R`, B2b's
response, B2c's response, B6's coefficient, **the coefficient of any additional independent invariant you
carried under §5**, `Λ_A⁰, Λ_V⁰, Λ_X⁰, τ_A, τ_V, τ_X`, `𝒜`, `μ_θ`, `μ_s`, the projection source term, and
the B0c face response. ⚠ Each response is a ratio — ⛔ state **what of what** before assigning a dimension,
and if a coefficient vanishes identically say its dimension is undetermined. Show each route and label it
**independent** or **definitional** — a route whose asserted equation *defines* the symbol under test is
definitional (⚠ a check reducing to `(X − 2Y) + 2Y == X` is worthless).

⭐⭐ **A HOMOGENEITY CHECK ON EVERY FINAL EQUATION, WITH UNITS RESTORED.** For each equation you report — the
in-plane equation, the thickness equation, the mass balance, the affinity `𝒜`, the closure, the face
response, the two-port power identity, and the dispersion determinant — verify that every additive term
carries the same `[L,T,M]` dimensions, and report the outcome per equation. ⛔⛔ **THE CHECK MUST BE ABLE TO
FAIL** — ⭐ **demonstrate it:** corrupt one term's dimensions deliberately, confirm the check reports a
failure, restore it, and report that demonstration.
⇒ `S11B_DIM_<name>`, `S11B_DIM_ROUTE_KIND_<name>`, `S11B_HOMOGENEITY_<equation>`,
  `S11B_HOMOGENEITY_ABLATION_DEMO`

### B8 · Controls — ⛔ FORM controls; a coefficient rescaling tests none of them

Projection / parity controls — ⭐ run each as a **controlled one-parameter family**, ⛔ not an arbitrary
asymmetric combination (an asymmetric component can lie orthogonal to the tested current, or an asymmetric
interval outside the window's support, moving nothing for reasons unrelated to the physics):

```
window     Ω_b(w) = sech²((w − b)/a)        b = 0 even; b ≠ 0 is the control
interval   [−L, L + c]                       c = 0 symmetric; c ≠ 0 is the control
```

- **control W** — `b ≠ 0`, `c = 0` (window shift only): recompute B0a's parity result as a function of `b`.
- **control I** — `b = 0`, `c ≠ 0` (interval asymmetry only): recompute it as a function of `c`.
⭐ For each, state whether the parity result is identically independent of that parameter — ⭐ this
distinguishes a parity **selection rule** from a domain-symmetry **artifact**. ⛔ Do not vary both at once.

Assembly controls (recompute the named tasks):
- **A — remove the thickness channel** (`δW = 0`): recompute B2c and B5. ⚠⚠ **A IS CONFOUNDED — report it as
  such.** Its exact cuts are `δW = 0 ⇒ e_W = 0, V_± = 0, Λ_V(ω)V_± = 0, δp_±V_± = 0, Λ_X(ω)𝒜_±V_± = 0`; but
  with the slab-side affinity active, `V_± = 0` does **not** imply `J_± = 0` or `δp_± = 0` — ⛔ do not delete
  that remaining transfer channel. ⇒ ⛔ Do not attribute any change under A to the thickness DOF alone; for
  each moving quantity report **which of the simultaneously-removed channels** it could be attributed to,
  and **say plainly if the attribution cannot be separated**.
- **B — remove the gradient stiffness** (`κ_W = 0`): recompute B2b and B5.
- **C — impermeable faces** (`Λ_A⁰ = Λ_V⁰ = 0`, `Λ_X⁰` symbolic): recompute B5.
- **D — remove the cross term** (`C = 0`): recompute B2c and B5.
- **E — remove the SLAB-SIDE part of the affinity `𝒜`** (bulk-pressure and velocity couplings untouched):
  recompute B5 and B2d. ⭐ This isolates the **new** transfer channel; ⛔ it is **not** the same cut as C.
- **F — remove the reciprocal traction** (`Λ_X⁰ = 0`, `Λ_A⁰`, `Λ_V⁰` symbolic): recompute B2a–B5 and both
  B2d questions. FORM cut — its recomputed mechanical operator must equal the symbolic `Λ_X⁰ = 0` operator
  of B2a and its power identity the B2d identity after the same substitution; report both symbolic
  differences. **Wrong derivation caught:** dropping `Λ_X` from the full system as well as the control makes
  F pass without testing the traction channel.

⭐⭐ **Recompute B6 under every one of A–D and report what moves.** ⛔ Do not assume a control cannot affect
B6, and ⛔ do not discard a dependence on the grounds it "must be" predetermined. ⚠ If none of A–D changes
B6, **state that as a finding and say why**. ⭐ For each control report which reported quantities move and
which do not; ⛔ report what each control does, ⛔ not what it was expected to do. ⭐ Keep `τ_A, τ_V, τ_X`
independent wherever the coefficient remains, and report which (if any) becomes irrelevant because its whole
channel was cut.
⇒ `S11B_CONTROL_WINDOW_PARITY`, `S11B_CONTROL_INTERVAL_SYMMETRY`, `S11B_CONTROL_NO_THICKNESS`,
  `S11B_CONTROL_A_ATTRIBUTION`, `S11B_CONTROL_NO_GRADIENT_STIFFNESS`, `S11B_CONTROL_IMPERMEABLE`,
  `S11B_CONTROL_NO_CROSS_TERM`, `S11B_CONTROL_NO_MU_COUPLING`, `S11B_CONTROL_NO_RECIPROCAL_TRACTION`,
  `S11B_CONTROLS_ON_TRANSVERSE`

### B9 · Validity — including the DERIVED convective correction

Report the conditions under which this step's linearisations hold, and **where in `(ω,k)` any fail**. State
separately any validity condition per response kernel in terms of `ωτ_A, ωτ_V, ωτ_X`; ⛔ do not replace them
by a common `ωτ` unless you explicitly impose the equal-time specialization.

⭐⭐ **THE DISCARDED CONVECTIVE CORRECTION — DERIVE it, ⛔ do not presuppose its form.** §2's acoustics are
linearised about rest. With a steady background normal drain `v_dr` the bulk operator becomes the **convective
wave equation** `(∂_t + v_dr ∂_w)² φ = c_s0² ∇₄² φ`. From this and the §1 ansatz, **derive** the leading
relative correction that the rest-frame linearisation of §2 discards, and **report where in `(ω,k)` it is
NOT small** — ⛔ do not state its order or its value in advance, and ⛔ do not presuppose any functional form
for it. ⚠ For complex `ω` and complex `q_out` a ratio of complex numbers has no
order relation — **state what modulus or norm you compare and whether the comparison is meaningful there;
if it is not, say so** rather than reporting a region computed from a measure you cannot justify. ⭐ Derive
**both** smallness conditions and emit them separately — (1) the background must not change appreciably
during a wave period (**timescale**); (2) the flow must be slow compared with `c_s0` (**speed**) — ⛔
condition 1 does not imply condition 2. ⚠ Solving the full convective problem is **out of scope**; naming
the discarded correction is the deliverable, and the rest-frame treatment is a **stated scope limit**.
⇒ `S11B_VALIDITY_CONDITIONS`, `S11B_VALIDITY_TIMESCALE`, `S11B_VALIDITY_FLOW_SPEED`,
  `S11B_DISCARDED_CONVECTIVE_CORRECTION`, `S11B_VALIDITY_FAILURE_REGION`

---

## 10 · Tag grammar — ⛔ both engines must produce PARALLEL tag sets

```
<ENGINE>_S11B_<QUANTITY>
```

- `<ENGINE>` is `WL` or `PY`. ⛔ Both prefixes are mandatory.
- Locus base-tags expand to the §8 suffixes (`_EQUATIONS`, `_SOLUTION`, …).
- ⭐ Emit **one line per tag**: `TAG: <payload>`, payload fully re-parseable (`InputForm` in Wolfram;
  `sympy.srepr`-safe string or `str()` in SymPy).

### ⛔⛔ ONE TAG PER NAMED OBJECT

⚠⚠ **This is the defect that made the two historical engine pairs uncomparable.** ⭐⭐ **Emit one tag per
named object.** ⛔ Do not bundle two named objects into one payload, and ⛔ do not split one named object
across several payloads. Where you emit an object for which this file supplies no tag name, choose one and
list it in the §13 report.

### ⭐ Engine-local tags — declare them, so parity is meaningful

Some tags **cannot** exist in both engines — §11's cross-step comparison has no counterpart in the engine
that does not import the chain's `LEDGER`, and each CAS emits its own solver-condition tags. ⛔ Those are
**not** disagreements. ⭐ Give every such tag the infix `_LOCAL_` immediately after the engine prefix
(`PY_S11B_LOCAL_…`), and emit one tag per engine listing every `_LOCAL_` name it produced.

---

## 11 · Inherited inputs, and the engine-local cross-step comparison

⭐ This step consumes objects the export chain already carries. ⛔ **Neither engine re-declares them under a
new identity.** ⚠ **Only the importing engine (the one wired to `research/pde_ledger_v3/scripts/S11_exports.py`)
performs the lookup and comparison below** — it is engine-local (`_LOCAL_` infix, §10). The blind engine
re-derives every dimension from §1–§6 and emits **no** import tag; that is by construction, ⛔ not a
denylist.

**The inherited objects and their standard names:**

| this file's symbol | ledger name | role |
|---|---|---|
| `c_s0` | `c_s0` | bulk sound speed — inherited |
| `μ_R` | `mu_R` | curl/twist modulus — inherited |
| `ρ_br⁰` | `rho_br` | ⭐ **settled here:** `ρ_br⁰ ≡ ρ_4D⁰ W₀` **is** the imported `rho_br` (S11's in-plane inertia coefficient; ⭐ S11's wall width is frozen, so `rho_br` already **is** the slab-integrated 3-inertia). ⛔ Bind to the imported object under the same name; ⛔ do not ask the engines to decide the identity, and ⛔ do not mint a second inertia knob. The cross-step residual (below) reports **consistency**, ⛔ not the identity |

⭐ `ρ_m` (bulk mass density) and `v_dr` (`v_bulk_normal_0`, the bulk normal drain) **originate in this step** —
they have **no** upstream row; ⛔ the importing engine must not manufacture one.

For each inherited object, the importing engine performs the live lookup
`LEDGER[name]`, and where the row exposes a `dimension_key`, the live
`LEDGER[name]['dimension_key']` → `LEDGER[that_key]['value']` lookup; it then emits the **derived** dimension
vector, the **imported** vector, their **difference**, and provenance (`class`, `step` of the coefficient
row; `class`, `step`, `corroborated_steps` of the resolved dimension row). ⚠ A row with **no**
`dimension_key` (as `c_s0` currently has none) resolves as **unresolved** — emit it in the unresolved list,
⛔ do not manufacture a placeholder vector. Emit `S11B_LOCAL_Q_RESOLVED_COEFFICIENTS` and
`S11B_LOCAL_Q_UNRESOLVED_COEFFICIENTS` from those live lookups, and `S11B_LOCAL_Q_RESIDUAL_SCOPE` with the
single token `CROSS_STEP_CONSISTENCY_ONLY` — ⛔ the imported vector is the predecessor's by construction, so
the residual tests cross-step consistency only and is **not** an independent operand.
⛔ If a derived and an imported vector disagree, **emit the disagreement and continue** — ⛔ do not adjust
the derivation to match, and ⛔ do not exit non-zero.

⚠ **Downstream, not this artifact.** The SymPy `S11b_exports.py` (flat `LEDGER` rows, `BUILD_INPUT_DIGESTS`
pinning `{the S11b audit, S11_exports.py, S11b_SHARED_PHYSICS.md}`, the three-valued `F9` collision, `D1`
iterate-the-emitted-collection, the `D3` in-run round-trip, the `_RELATIONALS` block) and the **separate**
cross-engine comparator (joining by name to the **frozen `T7` contract** —
`S11_C17_C18_spec_repair_decisions_v2.md:53-60`: paired-payload residuals, ⛔ **native-boolean rejection**,
three-valued undecided, plus a repoint ablation) are specified in the per-engine build directives and the
comparator artifact **per decision list `S11b_unified_decisions.md` G7–G8** — ⛔ not restated here. The
`F1`–`F9` chain rules (`S11_export_chain_decisions_v2.md`) and the blind-Wolfram control
(`S9_export_chain_rebuild_directive.md:16-18`) are **inherited** — a build directive points at them.

---

## 12 · What is supplied vs. what this build tests

| supplied — ⛔ NOT tested here | tested here |
|---|---|
| §2 bulk acoustics; the §1b branch prescription; the §4 closure, affinity `𝒜` and its normalization; the §5 term list (⛔ but its **closedness** is tested); the §6 balance-law route and virtual-displacement rule | everything in B0–B9 |
| `Λ_X⁰`, `τ_X` — a supplied reciprocal-traction channel (its **consequences** are derived, its **existence** is not) | the face response, the impedance and its regimes, the constraint, every EOM and response |
| the `μ_s = 0` equal-time reduction's **reference response** (the `S11B_ZPERM_SLICE` target of B2c) — ⭐ **held by the orchestrator, diffed on our side, ⛔ not printed** (⭐ the mapping `g_p` / `Λ_p⁰↔Λ_A⁰` is **computed from the supplied affinity**, ⛔ not withheld) | the longitudinal fate, the grazing classification, the transverse coupling, the passivity/reciprocity/admissibility regions, the stability boundary |
| `c_s0`, `μ_R` (inherited, §11); `ρ_br⁰` = `rho_br` **settled** (§11, frozen wall width) | every dimension and homogeneity outcome; the cross-step **consistency** residual (engine-local) |
| `v_dr` is a scope-limit parameter — the rest-frame linearisation is a stated limit | B9's derived discarded convective correction and both smallness conditions |

⭐ Where one engine can compute an object and the other cannot decide, the two emit typed, differing objects
and ⛔ **no computation in this build resolves that difference** — it is a finding for the orchestrator.

### ⭐⭐ Premise inventory — one named object

Each engine emits **one** engine-local tag, `<ENGINE>_S11B_LOCAL_PREMISE_INVENTORY`, listing the supplied
premises its run carried, in whatever form it holds them. It is **one named object**, so emitting it as one
tag is not the §10 bundling. Several premises assert an **absence** and have no live CAS object — ⛔ do not
manufacture one; list them as the engine holds them.

---

## 13 · Report back — ⛔ under 25 lines

1. The file you wrote, its line count, and the **total number of emitted tags**.
2. The tasks actually run, any skipped, with the runtime.
3. ⭐ Anything in this specification you found **ambiguous, under-determined, ill-posed, or tautological** —
   including any requested object that **cannot be computed** from what §1–§6 supplies, and any
   `<QUANTITY>` you could not emit under §10's one-tag-per-object rule. ⭐ This is wanted and is more
   valuable than a clean build.
4. ⛔ **Do NOT report what any value came out to be, do not interpret any result, and do not say whether
   anything "worked".** Your job ends at compute-and-print.
