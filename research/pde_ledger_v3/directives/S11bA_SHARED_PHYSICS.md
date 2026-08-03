# S11b-A — SHARED PHYSICS SPECIFICATION

⚠ **Inserted BYTE-IDENTICALLY into both the Mathematica and the SymPy directive.** It is the only part
they share, and it exists because on an earlier step the two engines were handed *different task lists*,
leaving the one surprising result with no second engine.

⚠⚠ **This is a SUB-STEP, deliberately narrow.** Two attempts to specify the whole brane–bulk interface in
one pass were rejected by review before any build ran — the second introduced new defects *caused by* its
own fixes. ⇒ ⭐ **The assembly is deferred to S11b-B and will be written against the objects THIS step
produces, not against a description of them.**

---

## 0 · Scope

**IN scope.** What the bulk does when the slab's faces move, and the identity that projects a
four-dimensional balance law onto the slab.

⛔⛔ **OUT of scope — do NOT attempt, and do NOT let any task drift into them.** The brane's in-plane
elastic sector; any modulus of the brane; the coupled brane+thickness+bulk dispersion; the transverse
mode's coupling; any statement about whether a mode radiates, is bound, is confined, or is observable.
⚠ **Those are S11b-B's.** ⭐ A tag here that answers one of them is a defect, not a bonus.

⛔ **Do not consult any file for the answers.** Everything needed is below. If something needed is
missing, emit the tag as `NOT_ESTABLISHED` and name what is missing.

## ⛔⛔ 0b · WHAT NEITHER ENGINE MAY READ — identical for both, and that is the point

⚠ **This list lives in the SHARED block deliberately.** A bar maintained separately in two headers drifts,
and the engine with the shorter list becomes a transcriber while the other derives — their agreement then
certifies a shared source rather than independent physics.

⛔ **Neither engine may read, open, grep, `cat`, `git show`, or `git cat-file`:**

- `research/pde_ledger_v3/scripts/` — **any** `.py`
- `research/pde_ledger_v3/mathematica/` — **any** `.wl`
- `research/pde_ledger_v3/reduction/` — any file, including the registry and the gates
- `research/pde_ledger_v3/steps/` — any file
- `research/pde_ledger_v3/V3_STEP_PLAN.md`
- `research/pde_ledger_v3/directives/S11b_*` — the superseded whole-interface directives
- `research/pde_audit/` — any file. ⚠ Prior work on adjacent physics whose current validity is **not
  established**.
- any file whose name contains `PREREGISTERED`
- **the other engine's deliverable for this sub-step**, by any route including git

⭐ **You do not need any of them.** Every input is in this specification. ⛔ If you believe an input is
missing, emit that tag as `NOT_ESTABLISHED` and name what is missing — ⛔ do not go looking.
⚠ `git status` and `git diff` are fine.

## 1 · Geometry, notation, and ⭐ CONVENTIONS

Four spatial dimensions `(x¹,x²,x³,w)` and time `t`; `x = (x¹,x²,x³)`. A slab of thickness `W₀` centred on
`w = 0`, with faces at `w = ±W₀/2`. Bulk on **both** sides. `w` is normal to the slab.

⭐⭐ **Conventions, fixed here so the two engines produce comparable signs.** ⚠ Reported real and imaginary
parts are meaningless for comparison unless these are shared.

| | convention |
|---|---|
| harmonic dependence | every perturbation `∝ exp[i(k·x − ωt)]`, with `ω > 0` and `k` real |
| face displacements | `ζ₊(x,t)`, `ζ₋(x,t)` = displacement of the **upper / lower** face, **both measured along global `+w`** |
| parity combinations | **thickness** `δW ≡ ζ₊ − ζ₋` · **centre shift** `ζ_c ≡ (ζ₊ + ζ₋)/2` |
| outward normals | upper face `n̂ = +ŵ` · lower face `n̂ = −ŵ` |
| outward face velocity | `V_± ≡ (∂_t ζ_±)·(±1)` — ⭐ **face-odd**; use this, ⛔ not the global `∂_t ζ_±`, wherever a quantity is defined along the outward normal |
| response ratio | `Z ≡ (pressure perturbation at a face) / (OUTWARD normal velocity of that face)` |
| ⭐⭐ radiation condition | in each half-space retain **only** waves carrying energy **away** from the slab (real normal wavenumber) or **decaying** away from it (imaginary). ⛔ **There are NO incoming waves from infinity.** State explicitly the sign of the normal wavenumber selected in each half-space under the harmonic convention above. ⚠ Boundedness alone does **not** select a branch when the normal wavenumber is real — both progressive waves are bounded, and 4 amplitudes against 2 face conditions is under-determined. |

⛔ **`ζ_c` is NOT the `h`-branon.** The ledger's `h` is a **dimensionless** field (`ξ_w = ℓh`); `ζ_c` here
is a **length**. ⚠ They are different objects and must never be identified without the normalisation that
relates them.

⛔ `x` is a position, ⛔ **never** a displacement. `Ω` (below) is a window, ⛔ never a thickness.

## 2 · The bulk sector

A scalar superfluid linearised to acoustics, with **no shear modulus**:

```
v = ∇₄ φ ,     δp = −ρ_m ∂_t φ ,     ∂_t² φ = c_s0² ∇₄² φ
```

`ρ_m` is the bulk mass density, `c_s0` the bulk sound speed. Both half-spaces `|w| > W₀/2` are bulk.

⛔ **Nothing about the slab's interior is specified and nothing about it is needed.** The slab enters this
step **only** through the motion of its faces.

## 3 · The window and the current

`Ω(w)` is a smooth **window function**, `≈1` inside the slab and `→0` outside, used to project
four-dimensional equations onto a three-dimensional description. Unless a task says otherwise take `Ω(w)`
**even** in `w`.

`j^w` is the `w`-component of the four-dimensional **MASS** current. ⛔ Not a number current — fixed here,
because the two differ by a conversion factor that would change A1's reported dimensions.

⭐ **Relative flux.** Where material crosses a **moving** face, the physically meaningful flux is measured
**relative to that face and along its outward normal**:

```
J_± ≡ ρ_m ( v_w − ∂_t ζ_± ) · (±1)        evaluated at the upper (+) / lower (−) face
```

⚠ **Note the plain minus.** The outward-normal sign appears **once**, in the trailing `(±1)`; applying it
again to the face velocity would give a face co-moving with the bulk a nonzero relative flux.

⚠ Report explicitly which signed combination of `J₊` and `J₋` corresponds to **net accretion by the slab**
and which to **through-flow**. ⛔ Do not use a sum of global-`w` fluxes in place of relative ones; the two
differ, and they carry different physics.

## 4 · Background state

There is a **steady background transfer** of material across the interface. Let `v₀` be the resulting
steady normal flow speed in the bulk near a face (**symbolic**).

⚠⚠ **§2's acoustics are linearised about REST, and that is an approximation, ⛔ not an identity.** A
nonzero `v₀` adds convective terms and a first-order flux `δρ·v₀`. ⭐⭐ **TWO independent smallness
conditions are involved and they are NOT the same:**

1. the background must not change appreciably during a wave period — a **timescale** condition;
2. the background flow must be slow compared with the bulk sound speed — a **speed** condition.

⭐ **Derive BOTH**, state each as an inequality between named quantities, and **report the order in
`v₀/c_s0` of the leading correction** that §2's rest-frame linearisation discards. ⛔ Do not assume either
condition holds, and ⛔ do not treat condition 1 as implying condition 2.
⚠ Solving the full convective problem is **out of scope**; naming the discarded order is the deliverable,
and the rest-frame treatment is a **stated scope limit**.

---

## TASKS

⛔⛔ **RULES THAT OVERRIDE EVERY TASK.**
1. **Every reported value must come from a computation in the script.** ⛔ A printed assertion is not a
   result.
2. ⛔ **If a task cannot be completed, emit its tag as `NOT_ESTABLISHED` and say which input is missing.**
   **A refusal is a valid and valuable output.**
3. ⛔ **Never silently choose a branch, a closure, or a convention.** If a choice is required and is not
   fixed above, introduce a **free symbol**, say so, and report the answer as a function of it.

---

**A1 · Projection identity.**
Integrate the four-dimensional continuity equation against `Ω(w)` over `w`, integrating by parts to
isolate the term carrying `j^w`. Report the resulting source term for a **finite** interval `[w₁, w₂]` and
for an **infinite** one.
⇒ `S11BA_PROJECTION_FINITE`, `S11BA_PROJECTION_INFINITE`

**A2 · Parity.**
With `Ω(w)` even, evaluate A1's source term for `j^w(w)` (a) **even** in `w`, (b) **odd** in `w`.
⭐ **Report the interval used and whether it is symmetric about `w = 0`**, and state for each result whether
it is **exact**, and **on what interval** — ⛔ an oddness argument does not by itself fix an integral over
an asymmetric interval.
⇒ `S11BA_PARITY_EVEN_JW`, `S11BA_PARITY_ODD_JW`, `S11BA_PARITY_INTERVAL`

**A3 · Dynamical window.**
Repeat A1 with `Ω = Ω(w; x, t)`. Enumerate **every** term present here and absent from A1.
⇒ `S11BA_DYNAMIC_WINDOW_EXTRA_TERMS`

**A4 · Bulk response to moving faces — impermeable.**
Solve §2 with both faces displaced as `ζ₊`, `ζ₋`, imposing that the bulk normal velocity at each face
equals that face's normal velocity. Treat both half-spaces, applying §1's radiation condition. Report `Z`
(as defined in §1) and its real and imaginary parts, **separately for the `δW` and `ζ_c` combinations**
of §1.
⭐ **Inertial loading is defined here so both engines report the same object:** in a regime where `Z` is
purely imaginary, `m_add` is the coefficient satisfying `δp|_face = m_add · ∂_t² ζ_±`, reported **PER
FACE**. ⛔ Do not report a two-face sum, and ⛔ do not fold in any half-amplitude convention from `δW` —
those combinations are assembly data and belong to a later step.
⭐⭐ **There are THREE regimes** — the bulk normal wavenumber squared may be positive, negative, **or
zero**. Report all three, including the behaviour of every reported quantity **at** the zero (grazing)
case, where some may be singular. ⚠ Omitting that third case is a known prior defect of this corpus.
⇒ `S11BA_Z_IMPERMEABLE`, `S11BA_Z_BY_REGIME`, `S11BA_Z_BY_PARITY`, `S11BA_ADDED_MASS`,
  `S11BA_GRAZING_BEHAVIOUR`

**A5 · Bulk response — permeable faces.**
⭐ **Closure, fixed here so the family is finite.** Take the relative flux of §3 to obey a linear law at
each face with **one relaxation time**, since the transfer is a conversion process occurring at a finite
rate:

```
J_±  =  Λ_p(ω) · δp|_face  +  Λ_V(ω) · V_±

Λ_p(ω) = Λ_p⁰ / (1 − iωτ)          Λ_V(ω) = Λ_V⁰ / (1 − iωτ)
```

⛔⛔ **DO NOT set `τ = 0`, and do not treat an instantaneous law as the default.** ⚠ `τ → 0` is the
**memoryless limit**, and it asserts that conversion across the interface is *instantaneous* — which is
close to the quantity this programme exists to determine. Imposing it would answer a rate question by
assumption, and a rational response looks equally healthy either way.
⭐ **Report the `τ → 0` limit as a special case**, ⛔ not as the premise.

⚠ **Report the behaviour across the whole range of `ωτ`** — small, order unity, and large — for every
quantity whose value depends on it. ⛔ Do not report only one regime.

⚠ **Scope limit to record:** a **single** `τ` shared by both coefficients. A more general law would admit
separate relaxation times or a full memory kernel; that generalisation is ⛔ not attempted here.

⭐ **`V_±` is the OUTWARD face velocity of §1, ⛔ not the global `∂_t ζ_±`.** ⚠ `J_±` is face-odd by
construction, so pairing it with a face-even quantity through a single scalar coefficient would give
reflection-related faces opposite constitutive terms and spuriously mix `δW` into `ζ_c`.

⭐ **`Λ_p⁰`, `Λ_V⁰` and `τ` are free but REAL symbols, with `τ ≥ 0`.** ⛔ Do not carry them as generic
complex quantities — a dissipation classification asks which terms are real and in phase with velocity,
and that question is vacuous if the underlying constants may themselves be complex. ⚠ The frequency
dependence above is the **only** source of complexity in the closure.

Report their dimensions. Recompute A4 under this condition and report the modified `Z` **as a function of
`Λ_p⁰`, `Λ_V⁰` and `ωτ`**. Then report, **for each of A4's three regimes and each parity combination**,
whether a dissipative part (real, in phase with velocity) is present, **on which coefficient it depends**,
and **how it varies with `ωτ`** — including both limits.

⛔⛔ **REPORT THE DEGENERATE LOCI.** The closure is only *generically* solvable: there are coefficient
values at which a face equation loses its dependence on the bulk amplitude, so a driven face has **no**
solution or an undriven one has a **free** amplitude. ⚠ A generic rational `Z` conceals both and still
prints `PASS`. ⇒ solve for those loci explicitly and report them.

⛔ Do not select a value for either coefficient. ⛔ Do not report a single "the" permeable response.
⚠ **Record the restriction to a memoryless law as a stated scope limit** — a more general law would admit
memory, and that generalisation is not attempted here.
⇒ additional tag `S11BA_DEGENERATE_LOCI`
⇒ `S11BA_PERMEABLE_COEFF_DIMS`, `S11BA_Z_PERMEABLE`, `S11BA_DISSIPATIVE_BY_REGIME_AND_PARITY`,
  `S11BA_DISSIPATION_VS_OMEGA_TAU`, `S11BA_TAU_ZERO_LIMIT`, `S11BA_CLOSURE_SCOPE_LIMIT`

**A6 · Dimensions.**
Derive, from the equations above and ⛔ **not** from any external table or registry, the `[L,T,M]`
exponents of: `Z`, `m_add`, `ρ_m`, `c_s0`, `v₀`, `Λ_p⁰`, `Λ_V⁰`, `τ`, and A1's source term. Show the route
for each.
⚠ **A check reducing to an identity of the form `(X − 2Y) + 2Y == X` is worthless.** Label each route
**independent** or **definitional**.
⇒ `S11BA_DIM_<name>` per entry, plus `S11BA_DIM_ROUTE_KIND` per entry

**A7 · Controls.** ⛔ FORM controls; ⛔ do not substitute a coefficient rescaling.

⭐⭐ **Use these CONTROLLED one-parameter families, ⛔ not an arbitrary member of a class.** A control
that merely says "make it asymmetric" has no determinate effect — an asymmetric component can be chosen
orthogonal to the tested current, or an asymmetric interval can lie outside the window's support, and the
control then moves nothing for reasons that have nothing to do with the physics.

```
window     Ω_b(w) = sech²((w − b)/a)        b = 0 is even; b ≠ 0 is the control
interval   [−L, L + c]                       c = 0 is symmetric; c ≠ 0 is the control
```

- **A** — `b ≠ 0`, `c = 0`; recompute A2.
- **B** — `b = 0`, `c ≠ 0`; recompute A2.

Report A2's result **as a function of `b` and of `c`**, and state for each whether it is identically
independent of that parameter or not. ⭐ That is what distinguishes a parity selection rule from a
domain-symmetry artifact. ⛔ Report what each control does, ⛔ not what it was expected to do.
⇒ `S11BA_CONTROL_WINDOW_PARITY`, `S11BA_CONTROL_INTERVAL_SYMMETRY`

**A8 · Validity conditions.** Emit **both** of §4's derived conditions — the timescale one and the flow-speed
one — as separate tags, plus the order in `v₀/c_s0` of the leading correction discarded by §2's rest-frame
linearisation. ⛔ Do not emit one condition in place of two.
⇒ `S11BA_VALIDITY_TIMESCALE`, `S11BA_VALIDITY_FLOW_SPEED`, `S11BA_DISCARDED_CONVECTIVE_ORDER`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`. Values must be explicit expressions, not prose, wherever mathematical. End
with a single `VERDICT:` line reporting whether the script's own internal consistency checks contradicted
each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ It is not a verdict
on the physics and must not be worded as one.
