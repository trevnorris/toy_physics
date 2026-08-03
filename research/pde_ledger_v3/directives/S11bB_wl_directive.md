# DIRECTIVE — S11b-B blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bB_interface_assembly_mathematica_audit.wl`

Run with `math -script <path>`. Iterate until it completes with no errors and no unevaluated output. Then
stop and exit — ⛔ do not write a report, a summary, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It is an **independent** check on a SymPy audit that does not yet exist. If it agrees with that audit
because it copied from something, the check is worthless and the step is not verified.
⇒ **The read-bar list is §0b of the shared specification below**, byte-identical to the other engine's.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads.**
- Prefix every output tag with `WL_`.
- Strip `ConditionalExpression[0, …]` when checking that something vanishes; test poles with
  `1/expr == 0`.
- Keep total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically rather
  than raising a limit.

---

# S11b-B — SHARED PHYSICS SPECIFICATION (rev 4)

⚠ **Inserted BYTE-IDENTICALLY into both engine directives.** It is the only part they share.

⚠⚠ **Revision 4.** Three revisions were rejected by independent review **before any build ran**. Rev 1
mandated a **non-uniform** background while fixing every perturbation to a plane wave, and those cannot
both hold — position-dependent coefficients mix wavevectors, so a global dispersion relation does not
exist. ⇒ ⭐ **The non-uniform problem is now a separate step (S11b-C).** This step is the **homogeneous
assembly**, which is a clean Fourier problem and answers a question that needs no gradients at all.

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

⭐⭐ **COMPLEX FREQUENCY — the deliverable is an imaginary part, so the branch matters.**
Treat `ω` as complex. ⛔⛔ **NO BRANCH RULE IS SUPPLIED HERE. DERIVE IT.** An earlier revision supplied one
(`Im q_out ≥ 0`); it was **wrong** — it is not the retarded continuation, and both engines would have
implemented it and agreed on a deliverable of the wrong sign.

⭐ **The physical requirements, which are what you must satisfy:**
1. For **real** `ω` in the propagating regime (`q² > 0`), the bulk solution must carry energy **away** from
   the slab on both sides.
2. For **real** `ω` in the evanescent regime (`q² < 0`), it must **decay** as `|w| → ∞`.
3. The **retarded** response is analytic in the **upper half** complex-`ω` plane. The value at complex `ω`
   is the analytic continuation of the real-axis choice, reached from `Im ω > 0`.

⇒ ⭐⭐ **Derive the resulting prescription for `q_out(ω,k)` from 1–3, state it explicitly, and justify it.**
⛔ Do not assert it. ⛔ Do not switch sheets mid-calculation.

⭐⭐ **AND MAKE THE AMBIGUITY MEASURABLE:** report whether `S11BB_IMAGINARY_PART` **changes sign or
magnitude** under the opposite branch choice. ⚠ If it does, that dependence is itself a reportable result —
it says the deliverable rests on the continuation and not only on the algebra.
⇒ `S11BB_BRANCH_PRESCRIPTION`, `S11BB_BRANCH_JUSTIFICATION`, `S11BB_BRANCH_SENSITIVITY`

⭐ If a root of B5 requires continuation onto a second sheet, ⛔ do not silently follow it — report that it
leaves the physical sheet, and say whether the object found there is a normal mode or a resonance.

## 2 · Established input from S11b-A — ⛔ do NOT re-derive

```
q² = ω²/c_s0² − k²            Z = ω ρ_m / q_out
q² > 0  q_out = q     Z real          radiation resistance
q² < 0  q_out = iα    Z = −iωρ_m/α    inertial loading ρ_m/α per face
q² = 0  Z singular
```

**Permeable faces.** ⭐ **Signs fixed here, because two engines can otherwise obtain opposite damping:**
`J_±` is the **mass flux per unit face area leaving the slab through that face**, measured along that
face's **outward** normal (so `J_± > 0` removes material from the slab). `δp` is the **bulk** pressure
perturbation evaluated at that face. Then

`J_± = Λ_p(ω)δp + Λ_V(ω)V_±`, `Λ(ω) = Λ⁰/(1 − iωτ)`, with `Λ_p⁰`, `Λ_V⁰`, `τ` real, free, and `τ ≥ 0`.

⛔ **PASSIVITY IS NOT ASSUMED, because from this closure alone it is not computable.** Deciding whether the
interface can only dissipate requires the brane–bulk affinity conjugate to `J_±` and the reciprocal
traction law, ⛔ **neither of which is supplied**. ⇒ ⭐ **Report the inequality on `Λ_p⁰`, `Λ_V⁰`, `τ` that
would be required for the interface not to inject energy, and state whether the information given is
sufficient to determine it.** ⚠ If it is not, emit `NOT_ESTABLISHED` and name exactly what is missing —
**a refusal here is the correct answer if it is the true one.** ⛔ Do not impose passivity to remove a
growing root; see §0.
⇒ `S11BB_PASSIVITY_CONDITION`, `S11BB_PASSIVITY_DETERMINABLE`

⚠⚠ **A declared limitation of this closure, which you must report the consequences of.** `δp` is the
**bulk** pressure perturbation at the face. A face at rest with an internal brane density perturbation has
`V_± = 0`, and with an unperturbed bulk `δp = 0`, hence `J_± = 0`. ⇒ **pressure-driven conversion by an
internal brane perturbation is excluded by the form of this closure, not by physics.** ⭐ Report what this
excludes and whether it affects your B5 and B6 conclusions.
⇒ `S11BB_CLOSURE_EXCLUDES`

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
(`u`, `∇u`, `θ`, `∇θ`, `e_W`, `∇e_W`), and construct **every** scalar quadratic in them allowed by the
symmetries of a uniform slab — isotropy in the three in-plane directions, and reflection `w → −w`. Compare
that basis against the list above. **For each invariant the list omits, state whether it is independent of
the listed terms or redundant with them *modulo B1's constraint*** — ⚠ the constraint relates these fields,
so a term can be symmetry-allowed and still carry no new physics. ⭐ **If any omitted invariant is
genuinely independent, carry it with a free symbolic coefficient and report its effect on B4 and B5.**
⇒ `S11BB_ENERGY_BASIS`, `S11BB_ENERGY_BASIS_OMISSIONS`, `S11BB_ENERGY_BASIS_INDEPENDENT_TERMS`

⚠ **`κ_W` is included because "restoring stiffness" alone is ambiguous** — a thickness stiffness may act on
`δW` or on `∇δW`, and the two give different `k`-dependence. ⭐ **Report how each result depends on `κ_W`,
and what changes if it vanishes.** ⛔ Do not set it to zero by default.

## ⭐⭐ 3b · HOW TO DERIVE THE EQUATIONS — the route is fixed, because inequivalent routes exist

⚠ B1's constraint carries a **source** (face flux) and **memory** (`τ`). ⇒ it must be **adjoined**,
⛔ **not substituted into `U` to eliminate a field** — substitution, and local stress plus continuity, give
**inequivalent** equations and different dispersions, and both would look like obeying "derive from §3".
⭐ **Report whether the constraint is holonomic** once `J_±` is treated as determined (see the accounting
rule below), and whether that changes what substitution is legitimate. ⛔ Report what you find; do not
assume either answer.

⭐ **Use this prescription and no other:**
1. Form the Lagrangian `L = T − U` with `T` the kinetic energy above, treating `u`, `δW` and `θ` as
   **independent** fields.
2. Adjoin B1's constraint with a **Lagrange multiplier** field `λ(x,t)`. ⭐ Report what `λ` turns out to be
   physically.
3. ⛔⛔ **NO RAYLEIGH DISSIPATION FUNCTION. NO DISSIPATIVE TERM IN `U`.** An earlier revision mandated a
   Rayleigh function built from `J_±` **on top of** `Z_perm`, which counts the same closure twice — see the
   accounting rule below.
4. The **bulk** enters as an **external generalized force on the faces**, obtained from the face pressure
   through `Z_perm` of §2. ⭐ **State explicitly which of the variations this force contributes to, and
   why** — ⛔ do not assume it appears in only one equation, and ⛔ do not assume it appears in all of them.
5. Vary independently with respect to `u`, `δW`, `θ` and `λ`, and report **all four** resulting equations,
   with the external face force added to whichever of them step 4 determines.

### ⭐⭐ THE SINGLE ACCOUNTING RULE — count every loss channel exactly once

⭐ **`J_±` is NOT an independent field.** It is **determined** by the closure of §2 evaluated on the bulk
solution: the face motion fixes `V_±`, the bulk solution fixes `δp`, and the closure then fixes `J_±`.
⇒ `Z_perm` **already** contains `Λ_p⁰`, `Λ_V⁰` and `τ`. The flux appearing in B1's constraint is **the same
determined object**, ⛔ not a second independent source.

⭐⭐ **THE CLAIM TO BE CHECKED — ⛔ it is a hypothesis here, not a premise — is that `Z_perm` is then the
ONLY place the interface enters the mechanical equations.** ⛔ Do not assume it; test it, and ⭐ **if it is
false, say so and say what the missing term is.** Do it explicitly:
- **Enumerate every channel by which energy leaves the slab** in this model.
- For each, show **where it appears** in your equations, and that it appears **exactly once**.
- ⚠ **If a channel exists that the face force through `Z_perm` does not capture — for instance one carried
  by the material crossing the face rather than by the work done on it — report it and say what would be
  needed to include it.** ⛔ Do not silently fold it into `Z_perm`, and ⛔ do not silently drop it.
- **Cross-check:** compute the time-averaged power delivered to the bulk from the face force, and compare
  it against the decay rate of the B5 root. ⭐ **Report whether they agree.** ⛔ If they do not, that is a
  result — report the discrepancy, do not tune either side to match.
⇒ `S11BB_LOSS_CHANNELS`, `S11BB_LOSS_COUNTED_ONCE`, `S11BB_POWER_BALANCE_CHECK`

⭐ **Report whether the IMPERMEABLE limit (`Λ_p⁰ = Λ_V⁰ = 0`) of this prescription agrees with direct
substitution of the constraint into `U`.** ⛔ If it does not, that is a result, not an error to smooth over.
⚠ **Note that this limit is NOT the lossless limit** — §2 trap 2 says radiation resistance survives it with
impermeable faces. ⛔ Do not describe `Λ⁰ → 0` as "reversible".

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

**B2 · The equations of motion.** From §3, derive the in-plane equation and the thickness equation,
including the force the bulk exerts on **both** faces via §2. Report both operators.
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
⇒ `S11BB_LONGITUDINAL_DISPERSION`, `S11BB_ROOTS`, `S11BB_IMAGINARY_PART`, `S11BB_DISSIPATION_ORIGIN`,
  `S11BB_ROOT_STABILITY_CLASS`, `S11BB_STABILITY_CONDITION`

**B6 · The transverse mode, computed.** On this uniform background, compute the coupling between the
transverse in-plane mode and the thickness degree of freedom **from B1's constraint and §3's energy**,
⛔ not by asserting a divergence-free argument. Report the coefficient and any modification to the
transverse dispersion. ⭐ **State explicitly what the coefficient couples to what**, and its normalization,
⛔ before assigning it a value or a dimension; ⚠ if it vanishes identically, say so and say that its
normalization is then undetermined. ⭐ Report whether the transverse mode acquires an imaginary part, and
its dependence on `Λ_p⁰`, `Λ_V⁰` and `ωτ` across the full range.
⚠ ⛔ **Whatever this returns, it does NOT settle whether confinement is unconditional** — that is a
non-uniform question and out of scope here. ⛔ Do not phrase it as if it does.
⇒ `S11BB_TRANSVERSE_COUPLING`, `S11BB_TRANSVERSE_DISPERSION`, `S11BB_TRANSVERSE_DISSIPATION`

**B7 · Dimensions.** Derive from the equations above, ⛔ not from any table, the `[L,T,M]` exponents of
`B_ρ`, `B_ρ⁽³⁾`, `μ_W`, `k_W`, `κ_W`, `C`, B3's response, B4's response, B6's coefficient, **and the
coefficient of any additional independent invariant you carried under §3**. ⚠ Each of those three responses is a ratio; ⛔ state what of what before assigning a dimension, and if a coefficient vanishes identically say that its dimension is undetermined. Show each route and
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
- **C — impermeable faces** (`Λ_p⁰ = Λ_V⁰ = 0`) and recompute B5.
- **D — remove the cross term** (`C = 0`) and recompute B4 and B5.
⭐⭐ **Recompute B6 under every one of A–D as well, and report what moves.** ⛔ Do not assume in advance
that a control cannot affect B6, and ⛔ do not discard a dependence you find on the grounds that it "must
be" algebraically predetermined. ⚠ If none of A–D changes B6's reported quantities, **state that as a
finding and say why** — that is a result about the structure of the coupling, and it must be **discovered
here, not assumed**.
⭐ For each, report which reported quantities move and which do not. ⛔ Report what each control does,
⛔ not what it was expected to do.
⇒ `S11BB_CONTROL_NO_THICKNESS`, `S11BB_CONTROL_A_ATTRIBUTION`, `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`,
  `S11BB_CONTROL_IMPERMEABLE`, `S11BB_CONTROL_NO_CROSS_TERM`, `S11BB_CONTROLS_ON_TRANSVERSE`

**B9 · Validity.** Report the conditions under which this step's linearisations hold, including §2's
background-flow condition, and **where in `(ω,k)` any fail**.
⇒ `S11BB_VALIDITY_CONDITIONS`, `S11BB_VALIDITY_FAILURE_REGION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`; explicit expressions wherever mathematical. End with a single `VERDICT:`
line reporting whether the script's own internal consistency checks contradicted each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ Not a verdict on the
physics.
