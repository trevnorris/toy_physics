# DIRECTIVE — S11b-B SymPy audit and registry insertion

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`

⛔ **One deliverable only.** ⚠ There is **no registry insertion and no gate run in this step** — both live
under the barred `reduction/`, and an earlier revision listed them here while barring them below, which
relocated the asymmetry instead of removing it.

Run the script. Iterate until it exits 0. Then stop and exit — ⛔ do not
write a report or a summary.

## ⛔⛔ ONE OF TWO INDEPENDENT ENGINES

A blind Mathematica audit of the same physics exists. **The disagreement between the engines is the test**,
so an agreement reached through a shared source certifies nothing.
⇒ **The read-bar list is §0b of the shared specification below**, byte-identical to the other engine's.

## ⛔⛔ NO REGISTRY WORK THIS STEP

⚠ Rev 1 promised registry access "after the physics". ⛔ **That ordering is not enforceable** — neither at
read time nor at import time — and the registry contains a relation this step is asked to reinterpret. ⇒ a
reviewer showed one engine could use it to choose an identification while the other stayed blind, making
apparent agreement non-independent.

⭐ **`reduction/` is therefore barred from BOTH engines this step** (§0b), and registry insertion is a
**separate pass after the physics is banked**. ⛔ Do not read, import, or modify it.

⚠ **B7's dimensions must be DERIVED from the specification's equations.** ⛔ Not read from any table.

## Script conventions

- One tag per line, `TAG: value`, no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.

---

# S11b-B — SHARED PHYSICS SPECIFICATION (rev 5)

⚠ **Inserted BYTE-IDENTICALLY into both engine directives.** It is the only part they share.

⚠⚠ **Revision 5.** Four revisions were rejected by independent review **before any build ran**. Rev 1
mandated a **non-uniform** background while fixing every perturbation to a plane wave, and those cannot
both hold — position-dependent coefficients mix wavevectors, so a global dispersion relation does not
exist. ⇒ ⭐ **The non-uniform problem is now a separate step (S11b-C).** This step is the **homogeneous
assembly**, which is a clean Fourier problem and answers a question that needs no gradients at all.

⭐ **Two things earlier revisions got wrong have since been settled by independent calculation and are now
SUPPLIED rather than left to you: the complex-frequency continuation (§1b) and the derivation route (§3b).**
⚠ Both were previously either asserted wrongly or left under-determined, and each produced a blocker.

---

## 0 · What this step is

Assemble the brane's in-plane sector, the slab's thickness degree of freedom, and the bulk's response into
one linear system on a **uniform** background, and determine the **longitudinal mode's fate**: does it
propagate freely, decay, **grow**, or fail to exist as a mode.

⭐⭐ **GROWTH IS AN ADMISSIBLE OUTCOME AND MUST BE REPORTED AS ONE.** The moduli in §3 are free symbols and
⛔ **no boundedness condition is imposed on them or on the cross term `C`**, so the quadratic form is not
assumed positive-definite and a root with `Im ω > 0` is a possible result of this calculation.
⛔ **If you find one, do not discard it, do not re-branch to remove it, and do not add a stability
assumption.** Report it, and report the condition on the moduli that would exclude it. ⚠ **A growing mode
is a first-class finding here, not an error to be cleaned up.**

⚠⚠ **BUT A GROWING ROOT CAN ALSO BE MANUFACTURED BY A MISTAKE**, and two specific mistakes do it:
a derivation-route error (§3b) and re-imposing the radiation condition at complex `ω` (§1). Each has a
**named, mechanical diagnostic** given in its section. ⭐ **Run both diagnostics before reporting any
growing root, and report their outcome alongside it.** ⛔ **This is not permission to discard growth** — it
is the requirement to **distinguish a real growing mode from an artifact**, which is what makes the
finding worth anything.

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

## ⭐⭐ 1b · COMPLEX FREQUENCY — the deliverable is an imaginary part, so the branch decides it

⚠⚠ **This section is supplied, and it is load-bearing.** `q_out` has branch points at `ω = ±c_s0|k|` sitting
**on the real axis**, and B5's roots are **below** it. Two descent paths that wind differently around a
branch point reach **different sheets**, where `q_out` differs by a factor of `−1`. That flips a decaying
**normal mode** into a growing **leaky resonance** at the *same* `ω`. ⛔ **The physical requirements below
do NOT by themselves fix this** — an earlier revision asked engines to derive the rule from them and it
was under-determined; a still earlier one asserted `Im q_out ≥ 0`, which is **non-analytic and wrong**.

**On the REAL axis, `q_out` is fixed by:**
1. `q² > 0` (propagating): the bulk solution carries energy **away** from the slab on both sides.
   ⚠ Read this as an **energy-flux** condition, valid for **both signs of `ω`**. ⛔ It is **not** a
   phase-velocity condition — that reading breaks for `ω < 0`.
2. `q² < 0` (evanescent): the solution **decays** as `|w| → ∞`.

**At COMPLEX `ω`, `q_out` is DEFINED as follows.** ⛔ Not derived, not re-selected — defined:

> `q_out(ω,k)` is the analytic continuation of its real-axis value reached from `ω + i0⁺` by moving
> **downward along the ray of fixed `Re ω`**. Equivalently: deform the inverse-Fourier-in-time contour
> downward from above the real axis while the branch points `ω = ±c_s0|k|` stay fixed and their cuts are
> dragged **vertically downward**, so that `q_out` is single-valued on the lower half-plane cut along
> `Re ω = ±c_s0|k|`.

⛔⛔ **REQUIREMENTS 1–2 MUST NOT BE RE-IMPOSED AT COMPLEX `ω`.** Whatever `|w| → ∞` behaviour the
continuation produces there is a **RESULT to report**, ⛔ never a criterion for re-selecting the root.
⚠⚠ **THIS IS THE DIAGNOSTIC §0 REFERS TO:** an engine that re-applies "must decay" at a complex pole lands
on the opposite sheet and turns a **damped resonance into an apparent instability**. ⇒ ⭐ **if you report a
growing root, state explicitly that you did not re-impose 1–2 to obtain it.**

⭐ **Verify, and report:** that this definition reproduces requirements 1–2 on the real axis. ⛔ If it does
not, report the disagreement rather than adjusting either side.
⭐ **Report the degenerate point** `ω = ±c_s0|k|`: there `q_out = 0`, the two bulk solutions **coalesce**
(the second going linear in `w`), and ⛔ neither requirement selects anything — continuity supplies it.
⭐ **If a root's trajectory crosses `Re ω = ±c_s0|k|` under parameter variation, report that it has LEFT
this sheet** — ⛔ do not re-select it onto one — and say whether the object is a normal mode or a resonance.

⭐⭐ **MAKE THE DEPENDENCE MEASURABLE:** report `S11BB_IMAGINARY_PART` **also on the opposite sheet**, and
report the ratio. ⚠ Expect it to be large — the deliverable rests on the continuation, not only on the
algebra, and that must be visible in the output rather than buried in a convention.
⇒ `S11BB_BRANCH_REALAXIS_CHECK`, `S11BB_BRANCH_DEGENERATE_POINT`, `S11BB_BRANCH_SENSITIVITY`,
  `S11BB_SHEET_OF_EACH_ROOT`

## 2 · Established input from S11b-A — ⛔ do NOT re-derive

```
q² = ω²/c_s0² − k²            Z = ω ρ_m / q_out
q² > 0  q_out = q     Z real          radiation resistance
q² < 0  q_out = iα    Z = −iωρ_m/α    inertial loading ρ_m/α per face
q² = 0  Z singular
```

⚠ **`Z` above is the BULK's OWN impedance:** the ratio of the bulk pressure perturbation at a face to the
**bulk material's** outward normal velocity there. ⛔ For an impermeable face that equals the face's own
velocity, but **with permeation they differ** — see below.

## ⭐⭐ 2b · PERMEABLE FACES — the closure, and the face response you must DERIVE

⭐ **Signs fixed here, because two engines can otherwise obtain opposite damping:** `J_±` is the **mass flux
per unit face area leaving the slab through that face**, measured along that face's **outward** normal (so
`J_± > 0` removes material from the slab). `δp` is the **bulk** pressure perturbation at that face.

```
J_± = Λ_p(ω) δp  −  Λ_μ(ω) μ_θ  +  Λ_V(ω) V_± ,        Λ(ω) = Λ⁰/(1 − iωτ)
```

with `Λ_p⁰`, `Λ_μ⁰`, `Λ_V⁰`, `τ` real, free, `τ ≥ 0`, and `μ_θ ≡ ∂U/∂θ|_{e_W}` the slab's own internal
chemical potential (`U` is §3; `μ_θ` is defined precisely in §3b).

⚠⚠ **WHY `Λ_μ` IS HERE, and it matters for how you read B5.** The thermodynamic conjugate force for mass
transfer across an interface is the **chemical-potential jump**, ⛔ not the bulk pressure alone. A closure
driven by `δp` only is **not Onsager-admissible**, and with it a face at rest carrying an internal density
perturbation would have `δp = V_± = J_± = 0` — ⇒ conversion driven by the slab's own state would be
**excluded by the form of the closure rather than by physics.** ⭐ `Λ_μ⁰` is carried **free**: `Λ_μ⁰ = 0`
recovers the earlier closure exactly, and ⛔ **nothing here asserts what value it takes.**

### ⭐ DERIVE the face response — ⛔ it is not supplied

Two relations close the face:
- **The bulk**, as above: `δp = Z · v_bulk,±` where `v_bulk,±` is the **bulk material's** outward normal
  velocity at that face.
- **Interfacial mass balance** (kinematics, ⛔ not a result): material leaving the slab enters the bulk, so
  the bulk's outward normal velocity exceeds the face's own by the volume flux of that material:
  `v_bulk,± = V_± + J_±/ρ_m`.

⭐ **Solve these together with the closure** and report the face pressure `δp` as a function of `V_±` and
`μ_θ`. ⚠ **Report it as what it is** — with `Λ_μ⁰ ≠ 0` the face response is ⛔ **not** a pure
velocity-impedance, because it acquires a term driven by a **brane** variable. ⭐ Report both coefficients
separately, and ⛔ do not force the result into the shape of a single impedance if it does not have it.
⇒ `S11BB_FACE_RESPONSE`, `S11BB_FACE_RESPONSE_MU_COEFF`

### ⭐⭐ ACCEPTANCE CHECK — an independently derived value your assembly must reproduce

At `Λ_μ⁰ = 0` your face response **must** reduce to the following, which was derived in a **separate step by
a different pair of engines**:

```
Z_perm = (ρ_m r + Λ_V⁰)/(y r − Λ_p⁰) ,    r = 1 − iωτ ,    y = q_out/ω
```

⭐ **Report the comparison explicitly.** ⛔ **If it does not reduce, report the disagreement — do not adjust
your derivation to match.** ⚠ A disagreement here is a real finding about one of the two steps.
⇒ `S11BB_ZPERM_REDUCTION_CHECK`

### Passivity — ⭐ TWO SEPARATE QUESTIONS. ⛔ Assert neither; compute both

1. **Face-port dissipativity.** Whether the face, as a port, absorbs rather than injects energy is an
   explicit inequality on the supplied symbols. ⭐ **Derive it and report it.**
2. **Thermodynamic admissibility of the closure.** Write the entropy production for interfacial transfer as
   a product of fluxes and their conjugate forces, and report the condition on `Λ_p⁰`, `Λ_μ⁰`, `Λ_V⁰`, `τ`
   for the corresponding **Onsager matrix to be positive semi-definite** — including any **reciprocity**
   relation it forces among them. ⭐ **Report whether the information given determines this**; if it does
   not, emit `NOT_ESTABLISHED` and name exactly what is missing.

⛔⛔ **DO NOT IMPOSE EITHER CONDITION TO REMOVE A GROWING ROOT.** ⭐⭐ **Instead, report the question that
matters: does a growing root survive INSIDE the thermodynamically admissible region?** ⚠ A growing root
that exists **only** outside it says nothing about the model; one that survives **inside** it is a
first-class result.
⇒ `S11BB_PORT_DISSIPATIVITY`, `S11BB_ONSAGER_CONDITION`, `S11BB_ONSAGER_DETERMINABLE`,
  `S11BB_GROWTH_INSIDE_ADMISSIBLE`

⛔⛔ **TWO TRAPS, both measured in S11b-A:**
1. The per-face inertial loading is `+ρ_m/α` **against the outward acceleration on BOTH faces**. The signed
   pair `(ρ_m/α, −ρ_m/α)` is an artifact of the global-`w` convention. ⛔ Do not consume it as an inertia.
2. ⛔ **Propagating `Re Z` is radiation resistance and exists with impermeable faces.** It is ⛔ not
   evidence of transfer through the interface. Only the **evanescent** `Re Z` is created by the closure.

⚠ **A third, on validity:** the bulk's rest-frame linearisation discards a relative correction of order
`v₀ q_n/ω`, where `v₀` is a steady background normal flow speed and **`q_n` is the bulk normal wavenumber —
the same object as `q_out` under the branch prescription of §1.** ⛔ The correction is **not** of order
`v₀/c_s0`; that is its value only on the sound cone. ⭐ Report where in `(ω,k)` the correction is **not**
small. ⚠ For complex `ω` and complex `q_n` a ratio of complex numbers has no order relation — **state what
modulus or norm you compare, and whether the comparison is meaningful there. If it is not, say so** rather
than reporting a region computed from a measure you cannot justify.

## 3 · The brane's stored energy

Per unit `x`-3-volume, with `e_W ≡ δW/W₀`:

```
U = ½ μ_R |∇×u|²  +  ½ B_ρ⁽³⁾ θ²  +  C W₀ θ e_W  +  ½ k_W W₀² e_W²  +  ½ κ_W W₀² |∇(δW)|²
```

and kinetic energy `½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)²`.

⭐⭐ **`C` is the symmetry-allowed CROSS term between densification and thickening, and it is included
deliberately.** ⚠ A diagonal energy would already impose that the two channels are separable — which is
part of what B4 is meant to determine, ⛔ not an input. **Report how every result depends on `C`, and what
changes when `C = 0`.** ⛔ Do not set it to zero by default.
⛔⛔ **THE LIST ABOVE IS THE SET OF TERMS CARRIED. IT IS NOT ASSERTED TO BE A CLOSED BASIS — CONSTRUCT THE
BASIS YOURSELF AND CHECK IT.** ⚠ An earlier revision gave this list and merely asked engines to "report
omissions"; that is too weak, because an omitted allowed term can change the dispersion while a control
that removes a *listed* term still reports cleanly.

⭐ **Do this as a task, before deriving anything:** enumerate the fields and their first gradients
(`u`, `∇u`, `θ`, `∇θ`, `e_W`, `∇e_W`) and construct **every** scalar quadratic in them allowed by the
symmetry group below. Compare that basis against the list above.

⛔⛔ **THE SYMMETRY GROUP, STATED IN FULL** — an earlier revision gave only isotropy and `w → −w`, which is
**not enough to define the basis**: under those alone a **pinning term `½K|u|²` is allowed**, and an engine
that includes it **gaps the modes** while one that does not finds them gapless. Both would be obeying the
text. The group is:
- **In-plane translation invariance** ⇒ `u` may enter **only through its gradients**, never undifferentiated.
- **In-plane isotropy _and_ parity** (the full `O(3)` acting on the three in-plane directions).
- **Reflection `w → −w`.**
- **Equivalence modulo total divergences** — two densities differing by a total in-plane divergence are the
  same term; ⛔ do not count both.
- ⛔ **Time-reversal is NOT assumed**, and no positivity or boundedness is assumed (§0).

⭐⭐ **JUDGE INDEPENDENCE AS FIELD BILINEARS, with `B1`'s constraint NOT applied.** ⚠ An earlier revision
asked for redundancy *"modulo B1's constraint"*; that has **no well-defined meaning** — `B1` is sourced,
carries memory, and changes rank at `ω = 0`, so "modulo B1" differs between the impermeable reduction and
the flux-on case, and two engines eliminating different terms would both be obeying the text.
⇒ **Carry EVERY independent invariant with a free symbolic coefficient**, and report its effect on B4/B5.
⭐ **Separately**, report which basis elements become redundant **once the constraint is applied**, and
⚠ **whether that set differs between the impermeable and the flux-on case** — ⛔ report the difference
rather than choosing one.
⇒ `S11BB_ENERGY_BASIS`, `S11BB_ENERGY_BASIS_OMISSIONS`, `S11BB_ENERGY_BASIS_INDEPENDENT_TERMS`,
  `S11BB_BASIS_REDUNDANCY_UNDER_CONSTRAINT`

⚠ **`κ_W` is included because "restoring stiffness" alone is ambiguous** — a thickness stiffness may act on
`δW` or on `∇δW`, and the two give different `k`-dependence. ⭐ **Report how each result depends on `κ_W`,
and what changes if it vanishes.** ⛔ Do not set it to zero by default.

## ⭐⭐ 3b · HOW TO DERIVE THE EQUATIONS — **BALANCE LAWS**. ⛔ NOT an action principle

⚠⚠ **Earlier revisions mandated a Lagrangian with a multiplier. That was the wrong tool, and it produced a
blocker twice.** The reason is structural: **a retarded kernel cannot be varied.** In any bilinear
`∫ λ 𝒴[∂_t δW] dt`, variation transposes the operator, whose symbol is `Y(−ω)` — the **advanced** kernel,
pole moved from `ω = −i/τ` to `ω = +i/τ`. ⇒ **stationarity of an action can never yield a retarded
response**; it symmetrises and keeps only the reactive part. `B1` is **not a constraint** at all — it is a
**balance law** whose source transfers mass irreversibly to the bulk.

⛔⛔ **THREE FORBIDDEN ROUTES.** Each is wrong, ⛔ not merely disfavoured:
- **(i) Substituting `B1` into `U` and then varying.** With `Λ⁰ ≠ 0`, `θ` is a **history functional**, not a
  state relation, so it may not be substituted into a stored energy. ⚠ Signature: an extra root, raising the
  dispersion determinant's **degree by one** — an invented mode.
- **(ii) Varying `J_±` inside a multiplier term.** Manufactures an anti-causal, **energy-generating**
  kernel. ⚠ Signature: a root at `ω ≈ +i/τ`.
- **(iii) Any route in which a response kernel is differentiated with respect to a field.**

⭐ **USE THIS AND NO OTHER:**
1. **State variables `u`, `δW`, `θ`.** ⛔ `θ` is **never eliminated**.
2. **Every constitutive quantity is a partial derivative of `U` at fixed remaining state variables:**
   `μ_θ ≡ ∂U/∂θ|_{e_W}` and `p_W ≡ ∂U/∂e_W|_θ`. ⛔ **No chain rules through the mass balance.**
3. **In-plane momentum balance** for `u`.
4. **Thickness equation: a force balance on the faces.** ⭐ **Every term in `U` that changes when the faces
   move contributes — INCLUDING the slab's own internal pressure pushing the faces apart.** ⛔ Do not omit a
   contribution on the grounds that it is "internal"; omitting it makes this route disagree with every other.
5. **`B1`'s mass balance as an INDEPENDENT EVOLUTION EQUATION**, ⛔ not a constraint to be substituted and
   ⛔ not adjoined with a multiplier.
6. **`J_±` and `δp` are PRESCRIBED RESPONSES.** Substitute them only **after** 3–5 are written. ⛔ They
   appear in no variation and are differentiated with respect to no field.

### ⭐⭐ THE CAUSALITY DIAGNOSTIC — mechanical, gradeable, ⛔ and it does NOT foreclose instability

Every response kernel in the final equations must be **retarded**: only `Λ(ω)`, `Z(ω)` and functions of them
may appear. ⛔⛔ **If `Λ(−ω)`, `Z(−ω)`, `Y(−ω)` or `Λ*(ω)` appears anywhere in your final equations, the
derivation is WRONG.** ⭐ **Run this check and report its outcome explicitly.**
⚠ **If you find a root at or near `ω = +i/τ`, treat it as a signal to re-check step 6 before reporting it**
— and then report **both** the root **and** the outcome of that re-check. ⛔ Do not simply delete it.
⇒ `S11BB_CAUSALITY_CHECK`, `S11BB_KERNEL_ARGUMENTS_PRESENT`

### ⭐ ENERGY ACCOUNTING — the discriminator that can actually fail

⚠ An earlier revision asked engines to enumerate loss channels **from the equations their own prescribed
route produced**, and to check power balance against it. ⛔ That is an identity for any system built that
way — **it could not fail.** Replace it with this:

⭐ **Compute `d/dt(T + U)` from YOUR equations.** Report **every** sink term, and **name the transport
process each one corresponds to**. ⛔⛔ **If a term corresponds to no transport process, REPORT IT** — that
is a **defect in the derivation**, ⛔ not a result about the physics.
⇒ `S11BB_ENERGY_SINKS`, `S11BB_UNATTRIBUTED_SINK_TERMS`

⛔⛔ **DO NOT USE THE IMPERMEABLE LIMIT AS A CHECK ON THE ROUTE.** All candidate routes **coincide** at
`Λ⁰ → 0`, so that limit **cannot discriminate between a correct derivation and any of the forbidden ones.**
⚠ Note also it is **not** the lossless limit — §2 trap 2: radiation resistance survives it. ⛔ Do not
describe `Λ⁰ → 0` as "reversible".

⛔ **The §3 list supplies no single in-plane compression modulus.** In it, compression is carried by `θ` and
by `e_W`, and how they combine is task **B4**. ⚠ Where a modulus measured with the thickness held fixed
would sit is an **output**, ⛔ not an input.

⚠⚠ **BUT THAT IS A PROPERTY OF THE LIST, ⛔ NOT AN ASSUMPTION YOU MAY IMPOSE.** If §3's basis construction
yields an **independent** invariant built from `∇·u`, ⭐ **carry it with a free coefficient as §3 instructs**
— ⛔ do **not** drop it to preserve the property above — and report how it changes B4's identification.
⚠ Both readings have been proposed; ⭐ **report which one your basis construction actually produces.**

---

## TASKS

⛔⛔ **RULES.** (1) Every value must be computed, ⛔ not asserted. (2) If a task cannot be completed, emit
`NOT_ESTABLISHED` and name the missing input — **a refusal is valid and valuable**. (3) ⛔ Never silently
choose a branch, closure, path, or expansion; introduce a free symbol and report the dependence.
(4) ⛔⛔ **No task states the form of its answer.** If a requested object turns out not to be a scalar — if
the response is operator-valued, or the polarizations split, or the effect is a spatial attenuation rather
than a frequency shift — ⭐ **report that instead**, and say so explicitly.

**B1 · The constraint.** ⛔⛔ **The exact balance is GIVEN — it is kinematics, not a result — and your job
is to linearise it.** ⚠ It is written out because three previous revisions left it verbal, and a verbal
statement lets the in-plane term be omitted, which produces a longitudinal mode with no restoring force
that **both engines would agree on**.

Let `Σ(x,t) ≡ ρ_4D(x,t) · W(x,t)` be the slab-integrated mass per unit `x`-3-volume, and let the slab
material move with in-plane velocity `v(x,t) = ∂_t u`. Conservation of slab material with a source from
face flux is

```
∂_t Σ  +  ∇_x·( Σ v )  =  −( J₊ + J₋ )
```

with `J_±` as defined in §2 (outward mass flux per unit face area, `J_± > 0` removes material) and `∇_x·`
the divergence in the **three** in-plane directions.

⭐ **Linearise this** to first order in `θ`, `e_W`, `u` and the flux, and report the resulting relation.
⭐⭐ **Report every term in your linearised constraint together with the term of the exact balance it came
from**, so that nothing is dropped without an account of why. ⛔ Do not add terms that are not consequences
of the balance above.

⭐ Then state how many **independent internal degrees of freedom** survive, and why. ⚠ **Define what you
are counting**: give the field list before the constraint is imposed and after, at fixed `(k, ω)`, and say
whether you are counting fields, amplitudes, or independent initial data.
⇒ `S11BB_CONSTRAINT`, `S11BB_CONSTRAINT_TERM_ORIGINS`, `S11BB_INTERNAL_DOF_COUNT`,
  `S11BB_DOF_COUNTING_CONVENTION`

**B2 · The equations of motion.** Following the route of §3b, derive the in-plane equation and the
thickness equation, including the force the bulk exerts on **both** faces via the face response you derived
in §2b. Report both operators.
⇒ `S11BB_INPLANE_EOM`, `S11BB_THICKNESS_EOM`, `S11BB_BULK_FORCE_ON_THICKNESS`

**B3 · Thickness response.** Solve the thickness equation for its response. ⭐ **State explicitly what the
response is a ratio of** (which output field to which driving quantity) and give its dimensions.
Then report **the bulk's contribution to the thickness operator in each regime of §2**, decomposed into a
part in phase with `∂_t²δW` and a part in phase with `∂_tδW`.
⛔⛔ **Do NOT report it as an "effective inertia".** ⚠ In the propagating regime the bulk load is radiation
**resistance**, which is velocity-like; calling it an inertia would smuggle damping into a mass and
collapse the distinction §2 trap 2 exists to preserve. ⭐ Report which regimes, if any, admit a mass
interpretation at all.
⇒ `S11BB_THICKNESS_RESPONSE`, `S11BB_RESPONSE_NORMALIZATION`, `S11BB_BULK_OPERATOR_BY_REGIME`,
  `S11BB_MASS_INTERPRETATION_VALID_WHERE`

**B4 · The compressional response.** Eliminate the thickness degree of freedom and report the in-plane
compressional stress response. ⭐ **State explicitly what this response is a ratio of** — which stress
component to which measure of deformation — **before** reporting it; ⚠ two engines can otherwise emit
incomparable objects under the same tag. ⚠ **Check whether B1's constraint changes rank at exactly `ω = 0`** — if an
integration constant survives there, say what fixes it, and report whether dividing by `ω` before or after
taking the limit changes the answer. ⭐ Report its behaviour in the limits `ω → 0` and `ω → ∞` **along a
stated path in `(ω,k)`** — ⛔ the limits need not commute, so name the path and report whether another path gives
a different answer. Then report **where a modulus measured with the thickness held fixed would sit**, or
that no consistent identification exists.
⇒ `S11BB_COMPRESSIONAL_RESPONSE`, `S11BB_LIMITS_AND_PATH`, `S11BB_FROZEN_THICKNESS_IDENTIFICATION`

**B5 · The longitudinal mode.** Assemble and report the dispersion relation. Report whether it admits a
closed-form `ω(k)`, whether roots are real, and for any complex root its imaginary part and **which
physical ingredient makes it nonzero** — ⛔ distinguishing the two mechanisms of §2 trap 2. ⭐ If a root
fails to exist as a normal mode, report that.
⭐⭐ **Report the SIGN of every imaginary part and classify each root as decaying or GROWING**, and report
the condition on the moduli and on `C` that separates the two. ⛔ **A growing root is a reportable result**
(§0) — ⛔ do not suppress it, re-branch to remove it, or assume the quadratic form is positive-definite.

⭐⭐ **FOR ANY GROWING ROOT, REPORT ALL THREE OF THESE — it is what separates a finding from an artifact:**
1. **The causality diagnostic of §3b** — did `Λ(−ω)`, `Z(−ω)`, `Y(−ω)` or `Λ*(ω)` appear anywhere? Is the
   root at or near `ω = +i/τ`?
2. **The sheet the root sits on**, and confirmation that requirements 1–2 of §1 were ⛔ **not** re-imposed
   at complex `ω` to obtain it.
3. ⭐⭐ **Whether it survives INSIDE the thermodynamically admissible region of §2b.** ⚠ A growing root
   existing **only outside** that region says nothing about the model; one surviving **inside** it is a
   first-class result. ⛔ Report which it is; ⛔ do not impose admissibility to make it disappear.
⇒ `S11BB_LONGITUDINAL_DISPERSION`, `S11BB_ROOTS`, `S11BB_IMAGINARY_PART`, `S11BB_DISSIPATION_ORIGIN`,
  `S11BB_ROOT_STABILITY_CLASS`, `S11BB_STABILITY_CONDITION`, `S11BB_GROWTH_ARTIFACT_DIAGNOSTICS`

**B6 · The transverse mode, computed.** On this uniform background, compute the coupling between the
transverse in-plane mode and the thickness degree of freedom **from B1's constraint and §3's energy**,
⛔ not by asserting a divergence-free argument. Report the coefficient and any modification to the
transverse dispersion. ⭐ **State explicitly what the coefficient couples to what**, and its normalization,
⛔ before assigning it a value or a dimension; ⚠ if it vanishes identically, say so and say that its
normalization is then undetermined. ⭐ Report whether the transverse mode acquires an imaginary part, and
its dependence on `Λ_p⁰`, `Λ_μ⁰`, `Λ_V⁰` and `ωτ` across the full range.
⚠ ⛔ **Whatever this returns, it does NOT settle whether confinement is unconditional** — that is a
non-uniform question and out of scope here. ⛔ Do not phrase it as if it does.
⇒ `S11BB_TRANSVERSE_COUPLING`, `S11BB_TRANSVERSE_DISPERSION`, `S11BB_TRANSVERSE_DISSIPATION`

**B7 · Dimensions.** Derive from the equations above, ⛔ not from any table, the `[L,T,M]` exponents of
`B_ρ`, `B_ρ⁽³⁾`, `μ_W`, `k_W`, `κ_W`, `C`, B3's response, B4's response, B6's coefficient, **and the
coefficient of any additional independent invariant you carried under §3**, plus `Λ_p⁰`, `Λ_μ⁰`, `Λ_V⁰`,
`μ_θ` and §2b's face response. ⚠ Each of those responses is a ratio; ⛔ state what of what before assigning
a dimension, and if a coefficient vanishes identically say that its dimension is undetermined. Show each route and
label it **independent** or **definitional** — a route whose asserted equation *defines* the symbol under
test is definitional.
⇒ `S11BB_DIM_<name>`, `S11BB_DIM_ROUTE_KIND_<name>`

**B8 · Controls.** ⛔ FORM controls; a coefficient rescaling tests none of them.
- **A — remove the thickness channel** (hold `δW = 0`) and recompute B4 and B5.
  ⚠⚠ **A IS CONFOUNDED AND YOU MUST REPORT IT AS SUCH.** Holding `δW = 0` also holds the faces still, which
  removes the face velocity, hence the bulk pressure load **and** the velocity-coupled permeability, all at
  once. ⇒ ⛔ **Do not attribute any change under A to the thickness degree of freedom alone.** ⭐ For each
  quantity that moves, report **which of the simultaneously-removed channels it could be attributed to**,
  and **say so plainly if the attribution cannot be separated** by this control.
- **B — remove the gradient stiffness** (`κ_W = 0`) and recompute B3 and B5.
- **C — impermeable faces** (`Λ_p⁰ = Λ_μ⁰ = Λ_V⁰ = 0`) and recompute B5.
- **D — remove the cross term** (`C = 0`) and recompute B4 and B5.
- **E — remove the chemical-potential coupling** (`Λ_μ⁰ = 0`, other couplings untouched) and recompute B5
  and the passivity results of §2b. ⭐ This is the control that isolates the **new** transfer channel; it is
  ⛔ **not** the same cut as **C**, which removes all three couplings at once.
⭐⭐ **Recompute B6 under every one of A–D as well, and report what moves.** ⛔ Do not assume in advance
that a control cannot affect B6, and ⛔ do not discard a dependence you find on the grounds that it "must
be" algebraically predetermined. ⚠ If none of A–D changes B6's reported quantities, **state that as a
finding and say why** — that is a result about the structure of the coupling, and it must be **discovered
here, not assumed**.
⭐ For each, report which reported quantities move and which do not. ⛔ Report what each control does,
⛔ not what it was expected to do.
⇒ `S11BB_CONTROL_NO_THICKNESS`, `S11BB_CONTROL_A_ATTRIBUTION`, `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`,
  `S11BB_CONTROL_IMPERMEABLE`, `S11BB_CONTROL_NO_CROSS_TERM`, `S11BB_CONTROL_NO_MU_COUPLING`,
  `S11BB_CONTROLS_ON_TRANSVERSE`

**B9 · Validity.** Report the conditions under which this step's linearisations hold, including §2's
background-flow condition, and **where in `(ω,k)` any fail**.
⇒ `S11BB_VALIDITY_CONDITIONS`, `S11BB_VALIDITY_FAILURE_REGION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`; explicit expressions wherever mathematical. End with a single `VERDICT:`
line reporting whether the script's own internal consistency checks contradicted each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ Not a verdict on the
physics.
