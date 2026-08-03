# S11b-B — SHARED PHYSICS SPECIFICATION (rev 2)

⚠ **Inserted BYTE-IDENTICALLY into both engine directives.** It is the only part they share.

⚠⚠ **Revision 2.** Rev 1 was rejected by two independent reviews before any build ran: it mandated a
**non-uniform** background while fixing every perturbation to a plane wave, and those cannot both hold —
position-dependent coefficients mix wavevectors, so a global dispersion relation does not exist.
⇒ ⭐ **The non-uniform problem is now a separate step (S11b-C).** This step is the **homogeneous
assembly**, which is a clean Fourier problem and answers a question that needs no gradients at all.

---

## 0 · What this step is

Assemble the brane's in-plane sector, the slab's thickness degree of freedom, and the bulk's response into
one linear system on a **uniform** background, and determine the **longitudinal mode's fate**: does it
propagate freely, decay, or fail to exist as a mode.

⛔ **OUT of scope.** Any statement about whether light's confinement is **unconditional** — that requires a
non-uniform background and is **S11b-C's**. ⭐ **Task B6 computes the transverse coupling on this uniform background**; whatever it returns, ⛔ that
result does **not** settle the unconditional question either way.

⛔ **Do not consult any file for the answers.** Everything needed is below.

## ⛔⛔ 0b · WHAT NEITHER ENGINE MAY READ

`research/pde_ledger_v3/scripts/` · `mathematica/` · `steps/` · `reduction/` · `V3_STEP_PLAN.md` ·
`directives/S11b_*` and `S11bA_*` · `research/pde_audit/` · anything named `PREREGISTERED`/`PREREG` ·
**the other engine's deliverable**.
⚠ `reduction/` is barred from **both** engines this step; registry insertion is a **separate later pass**,
so neither engine can identify a symbol the other is treating as independent. ⛔ Do not go looking.

## 1 · Geometry, fields, conventions

Four spatial dimensions `(x¹,x²,x³,w)`, `x = (x¹,x²,x³)`. A slab of thickness `W(x,t) = W₀ + δW` centred on
`w = 0`, faces at `w = ±W₀/2`, `D_brane = 3`. Bulk on both sides. **All background quantities are uniform
constants.**

| symbol | meaning |
|---|---|
| `u(x,t)` | in-plane displacement, **3 components**, ⛔ no `w`-component |
| `δW(x,t)` | thickness perturbation; faces at `ζ_± = ±δW/2` |
| `θ(x,t)` | **densification** — the Eulerian fractional perturbation of the brane material's local 4D density, `ρ_4D = ρ_4D⁰(1 + θ)`. ⚠ Eulerian, ⛔ not material |
| `ρ_br⁰ ≡ ρ_4D⁰ W₀` | background slab-integrated inertia density |
| `μ_R` | curl-type (twist) modulus, established input |
| `B_ρ` | local 4D compression modulus · `B_ρ⁽³⁾ ≡ B_ρ W₀` |
| `μ_W`, `k_W`, `κ_W` | thickness inertia, restoring stiffness, gradient stiffness — **symbolic** |
| `ρ_m`, `c_s0` | bulk mass density and sound speed |

**Conventions.** Every perturbation `∝ exp[i(k·x − ωt)]`. Face displacements measured along global `+w`;
outward normals `+ŵ` upper, `−ŵ` lower; outward face velocity `V_± ≡ (∂_tζ_±)(±1)`;
`Z ≡ (pressure at a face)/(OUTWARD normal velocity of that face)`.

⭐⭐ **COMPLEX FREQUENCY — fixed here, because the deliverable is an imaginary part.**
Treat `ω` as complex with the **retarded** prescription: all quantities are the analytic continuation from
`Im ω > 0` of their real-`ω` values. Fix the branch of `q` by continuation from the real axis, where it is
the outgoing (`q² > 0`) or decaying (`q² < 0`) root. ⛔ **State explicitly which sheet you are on**, and
⛔ do not switch sheets mid-calculation.

## 2 · Established input from S11b-A — ⛔ do NOT re-derive

```
q² = ω²/c_s0² − k²            Z = ω ρ_m / q_out
q² > 0  q_out = q     Z real          radiation resistance
q² < 0  q_out = iα    Z = −iωρ_m/α    inertial loading ρ_m/α per face
q² = 0  Z singular
```

**Permeable faces:** `J_± = Λ_p(ω)δp + Λ_V(ω)V_±`, `Λ(ω) = Λ⁰/(1 − iωτ)`, with `Λ_p⁰`, `Λ_V⁰`, `τ` real,
free, and `τ ≥ 0`. ⭐ **Assume passivity**: the closure may dissipate but ⛔ may not inject energy; report
the inequality on the coefficients this requires, and whether your results respect it.

```
Z_perm = (ρ_m r + Λ_V⁰)/(y r − Λ_p⁰) ,    r = 1 − iωτ ,    y = q_out/ω
```

⭐ **Use `Z_perm` throughout as the physical face response**; `Z` (impermeable) is the `Λ_p⁰ = Λ_V⁰ = 0`
special case and is the subject of control **C**.

⛔⛔ **TWO TRAPS, both measured in S11b-A:**
1. The per-face inertial loading is `+ρ_m/α` **against the outward acceleration on BOTH faces**. The signed
   pair `(ρ_m/α, −ρ_m/α)` is an artifact of the global-`w` convention. ⛔ Do not consume it as an inertia.
2. ⛔ **Propagating `Re Z` is radiation resistance and exists with impermeable faces.** It is ⛔ not
   evidence of transfer through the interface. Only the **evanescent** `Re Z` is created by the closure.

⚠ **A third, on validity:** the bulk's rest-frame linearisation discards a relative correction of order
`v₀|q_n|/ω`, where `v₀` is a steady background normal flow speed. ⛔ Not `v₀/c_s0`. Report where in `(ω,k)`
it is not small.

## 3 · The brane's stored energy

Per unit `x`-3-volume, with `e_W ≡ δW/W₀`:

```
U = ½ μ_R |∇×u|²  +  ½ B_ρ⁽³⁾ θ²  +  ½ k_W W₀² e_W²  +  ½ κ_W W₀² |∇(δW)|²
```

and kinetic energy `½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)²`.

⚠ **`κ_W` is included because "restoring stiffness" alone is ambiguous** — a thickness stiffness may act on
`δW` or on `∇δW`, and the two give different `k`-dependence. ⭐ **Report how each result depends on `κ_W`,
and what changes if it vanishes.** ⛔ Do not set it to zero by default.

⛔ **No single in-plane compression modulus is supplied.** Compression is carried by `θ` and by `e_W`, and
how they combine is task **B4**. ⚠ Where a modulus measured with the thickness held fixed would sit is an
**output**, ⛔ not an input.

---

## TASKS

⛔⛔ **RULES.** (1) Every value must be computed, ⛔ not asserted. (2) If a task cannot be completed, emit
`NOT_ESTABLISHED` and name the missing input — **a refusal is valid and valuable**. (3) ⛔ Never silently
choose a branch, closure, path, or expansion; introduce a free symbol and report the dependence.
(4) ⛔⛔ **No task states the form of its answer.** If a requested object turns out not to be a scalar — if
the response is operator-valued, or the polarizations split, or the effect is a spatial attenuation rather
than a frequency shift — ⭐ **report that instead**, and say so explicitly.

**B1 · The constraint.** Impose conservation of slab material per unit `x`-3-volume with a source from
face flux, and linearise. Report the relation among `θ`, `e_W`, `u`, and the flux. ⭐ State how many
independent internal degrees of freedom survive, and why.
⇒ `S11BB_CONSTRAINT`, `S11BB_INTERNAL_DOF_COUNT`

**B2 · The equations of motion.** From §3, derive the in-plane equation and the thickness equation,
including the force the bulk exerts on **both** faces via §2. Report both operators.
⇒ `S11BB_INPLANE_EOM`, `S11BB_THICKNESS_EOM`, `S11BB_BULK_FORCE_ON_THICKNESS`

**B3 · Thickness response.** Solve the thickness equation for its response to whatever drives it. Report
the response function and the **effective inertia** in each regime of §2.
⇒ `S11BB_THICKNESS_RESPONSE`, `S11BB_EFFECTIVE_INERTIA_BY_REGIME`

**B4 · The compressional response.** Eliminate the thickness degree of freedom and report the in-plane
compressional stress response. ⭐ Report its behaviour in the limits `ω → 0` and `ω → ∞` **along a stated
path in `(ω,k)`** — ⛔ the limits need not commute, so name the path and report whether another path gives
a different answer. Then report **where a modulus measured with the thickness held fixed would sit**, or
that no consistent identification exists.
⇒ `S11BB_COMPRESSIONAL_RESPONSE`, `S11BB_LIMITS_AND_PATH`, `S11BB_FROZEN_THICKNESS_IDENTIFICATION`

**B5 · The longitudinal mode.** Assemble and report the dispersion relation. Report whether it admits a
closed-form `ω(k)`, whether roots are real, and for any complex root its imaginary part and **which
physical ingredient makes it nonzero** — ⛔ distinguishing the two mechanisms of §2 trap 2. ⭐ If a root
fails to exist as a normal mode, report that.
⇒ `S11BB_LONGITUDINAL_DISPERSION`, `S11BB_ROOTS`, `S11BB_IMAGINARY_PART`, `S11BB_DISSIPATION_ORIGIN`

**B6 · The transverse mode, computed.** On this uniform background, compute the coupling between the
transverse in-plane mode and the thickness degree of freedom **from B1's constraint and §3's energy**,
⛔ not by asserting a divergence-free argument. Report the coefficient and any modification to the
transverse dispersion. ⭐ Report whether the transverse mode acquires an imaginary part, and its dependence
on `Λ_p⁰`, `Λ_V⁰` and `ωτ` across the full range.
⚠ ⛔ **Whatever this returns, it does NOT settle whether confinement is unconditional** — that is a
non-uniform question and out of scope here. ⛔ Do not phrase it as if it does.
⇒ `S11BB_TRANSVERSE_COUPLING`, `S11BB_TRANSVERSE_DISPERSION`, `S11BB_TRANSVERSE_DISSIPATION`

**B7 · Dimensions.** Derive from the equations above, ⛔ not from any table, the `[L,T,M]` exponents of
`B_ρ`, `B_ρ⁽³⁾`, `μ_W`, `k_W`, `κ_W`, B3's response, B4's response, B6's coefficient. Show each route and
label it **independent** or **definitional** — a route whose asserted equation *defines* the symbol under
test is definitional.
⇒ `S11BB_DIM_<name>`, `S11BB_DIM_ROUTE_KIND_<name>`

**B8 · Controls.** ⛔ FORM controls; a coefficient rescaling tests none of them.
- **A — remove the thickness channel** (hold `δW = 0`) and recompute B4, B5, B6.
- **B — remove the gradient stiffness** (`κ_W = 0`) and recompute B3, B5.
- **C — impermeable faces** (`Λ_p⁰ = Λ_V⁰ = 0`) and recompute B5, B6.
⭐ For each, report which reported quantities move and which do not. ⛔ Report what each control does,
⛔ not what it was expected to do.
⇒ `S11BB_CONTROL_NO_THICKNESS`, `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`, `S11BB_CONTROL_IMPERMEABLE`

**B9 · Validity.** Report the conditions under which this step's linearisations hold, including §2's
background-flow condition, and **where in `(ω,k)` any fail**.
⇒ `S11BB_VALIDITY_CONDITIONS`, `S11BB_VALIDITY_FAILURE_REGION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`; explicit expressions wherever mathematical. End with a single `VERDICT:`
line reporting whether the script's own internal consistency checks contradicted each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ Not a verdict on the
physics.
