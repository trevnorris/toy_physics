# DIRECTIVE — S11b blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_brane_bulk_interface_mathematica_audit.wl`

Run it with `math -script <path>`. Iterate until it runs to completion with no errors and no unevaluated
output. Then stop and exit — ⛔ do not write a report, a summary document, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It exists to be an **independent** check on a SymPy audit that does not yet exist. If it agrees with that
audit because it copied from something, the check is worthless and the step is not verified.

⛔ **DO NOT READ, open, grep, `cat`, `git show`, or otherwise inspect:**
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/` — any file, especially `quantities.yaml`
  and `relations.yaml`
- any `.py` file anywhere under `/var/projects/toy_physics/research/pde_ledger_v3/`
- `/var/projects/toy_physics/research/pde_audit/` — any file
- any file whose name contains `PREREGISTERED`

⭐ **You do not need any of them.** Every input is in the specification below. If you believe an input is
missing, emit the relevant tag as `NOT_ESTABLISHED` and name what is missing. ⛔ Do not go looking.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads.** This matches the 44
  Mathematica scripts already in this corpus, and their reliability comes from exactly that shape.
- Prefix every output tag with `WL_` — so `S11B_PROJECTION_IDENTITY_FINITE` is printed as
  `WL_S11B_PROJECTION_IDENTITY_FINITE`.
- Strip `ConditionalExpression[0, …]` wrappers when checking that something vanishes.
- To test whether a quantity has a pole, test `1/expr == 0` rather than `expr == Infinity`.
- Keep the total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically
  rather than raising a limit.

---

# S11b — SHARED PHYSICS SPECIFICATION (rev 2)

⚠ **This block is inserted BYTE-IDENTICALLY into both the Mathematica directive and the SymPy
directive.** It is the only part they share, and it exists because on the previous step the two engines
were handed *different task lists*, leaving the one surprising result with no second engine.

⚠ **Revision 2** — rev 1 was rejected by its own review before any build ran. It under-specified the
system in ways that would have made **both** engines return the same wrong answer. Fixes are marked ⭐.

---

## 0 · What this step is

Derive the **linear coupling between the brane and the bulk at their interface**, and from it the physical
status of the brane's longitudinal mode.

⛔ **Do not consult any file for the answers.** Everything needed is in this specification. If something
needed is missing, emit the relevant tag as `NOT_ESTABLISHED` and name what is missing, rather than
supplying it from elsewhere or inventing it.

## 1 · Geometry and notation

Four spatial dimensions with coordinates `(x¹, x², x³, w)` and time `t`. Write `x = (x¹,x²,x³)`.

The **brane** is a slab centred on `w = 0` of thickness `W`, with **two faces** near `w = ±W/2`. Its own
spatial dimension is `D_brane = 3`. The **bulk** is the region outside the slab. `w` is normal to the brane.

⭐ **Notation, fixed here to prevent collisions:**

| symbol | meaning |
|---|---|
| `x` | in-plane position, ⛔ **never** a displacement |
| `W(x,t) = W₀ + δW` | the slab **thickness** |
| `Ω(w)` | the **window function** used for projection — ⛔ a different object from `W` |
| `ζ₊(t)`, `ζ₋(t)` | **normal displacement of the upper / lower face** |
| `u(x,t)` | in-plane brane displacement, 3 components |
| `h(x,t)` | normal displacement of the brane as a whole |

⭐ **The two faces are independent, and their parity is the physics.** Define
`ζ₊ − ζ₋` and `ζ₊ + ζ₋` and state which of `δW` and `h` each corresponds to. ⛔ Do not model a single face;
a task that drives only one face cannot represent either mode correctly.

## 2 · The brane sector (established input, from earlier steps)

An in-plane displacement `u(x,t)` with **three** components `u_i`, `i = 1..3`.
⛔ **`u` has NO component along `w`.** Normal displacement is `h`, a different field; never combine them.

```
ρ_br ω² a_i = M_ij a_j ,      M_ij = μ_R (k·k) δ_ij + (B_comp⁰ − μ_R) k_i k_j
```

`ρ_br` is an inertia density, `μ_R` a curl-type modulus, `B_comp⁰` a compression modulus. These are
**slab-integrated**: `ρ_br = ρ_4D·W₀`, and likewise for the moduli, where the `4D` quantity is local.

⭐⭐ **`B_comp⁰` carries a superscript zero because it was measured in a calculation where the thickness
was NOT a degree of freedom.** This step restores that degree of freedom. ⛔ **Do not assume whether
`B_comp⁰` already contains the thickness channel or does not** — determining that is task **W10**, and
assuming either way double-counts or omits a stiffness that neither engine can then see.

## 3 · The bulk sector

A scalar superfluid, linearised to acoustics. It has **no shear modulus**.

```
v = ∇₄ φ ,     δp = −ρ_m ∂_t φ ,     ∂_t² φ = c_s0² ∇₄² φ
```

`ρ_m` is the bulk mass density, `c_s0` the bulk sound speed. The bulk occupies **both** sides of the slab.

## 4 · The thickness degree of freedom

`δW(x,t)` is a dynamical degree of freedom obeying, per unit `x`-3-volume,

```
μ_W ∂_t²(δW)  =  − k_W δW  +  F_face  +  F_int
```

- `μ_W` — inertia, **symbolic**
- `k_W` — restoring stiffness, **symbolic**. ⛔ **DO NOT assume a functional form or a `k`-dependence.**
- `F_face` — force per unit 3-area exerted by the bulk on the faces, from §3. ⭐ Both faces contribute.
- `F_int` — the internal driving from the brane sector, determined by §5's constraint.

⭐ This equation is given so that `[μ_W]` and `[k_W]` are **defined**; ⛔ it is not a claim about the
values of any term, and `F_face`/`F_int` must be **derived**, not posited.

## 5 · ⭐⭐ THE CONSTRAINT THAT COUPLES THE SECTORS

⚠ **rev 1 omitted this entirely, and without it the linear system is block-diagonal and the step returns a
false null.**

The slab's mass per unit `x`-3-volume is `ρ_br = ρ_4D·W`. Conservation of slab material, per unit
`x`-3-volume, with a source from material crossing the faces:

```
∂_t(ρ_4D W)  +  ∇·(ρ_4D W ∂_t u)  =  𝒮
```

where `𝒮` is the flux of material across the two faces per unit `x`-3-volume, and `∇` is the in-plane
divergence. ⭐ **Linearise this.** It relates `∇·u`, `δW`, `δρ_4D` and `𝒮`, and it is the **only** linear
channel connecting the brane sector to the thickness degree of freedom. ⛔ Do not introduce any other
coupling by hand.

**Energy storage.** Compressing the 4D brane material costs energy through its equation of state, with a
local modulus `B_ρ` (**symbolic**). Changing the thickness costs energy through the wall, with the
stiffness `k_W` of §4. ⛔ **Do not assume how a given compression divides between the two channels** —
that division is a consequence of the constraint above and is part of **W10**.

## 6 · The window and the current

`Ω(w)` is a smooth **window function**, `≈1` inside the slab and `→0` outside, used to project
four-dimensional equations onto the three-dimensional brane description. Unless a task says otherwise take
`Ω(w)` **even** in `w`. ⛔ `Ω` is **not** the thickness `W`.

`j^w` denotes the `w`-component of the four-dimensional mass (or number) current.

## 7 · Background state

There is a steady background transfer of material across the interface. Treat the background as **constant
over one wave period** and work with the **perturbation** about it.

⭐ **Derive the condition under which that treatment is valid**, state it as an inequality between named
timescales, and emit it as a tag. ⛔ Do not assume it holds.

---

## TASKS

⛔⛔ **RULES THAT OVERRIDE EVERY TASK BELOW.**
1. **Every reported value must come from a computation performed in the script.** ⛔ A `Print` of an
   asserted string is not a result.
2. ⛔ **If a task cannot be completed, emit its tag with the literal value `NOT_ESTABLISHED` and state
   which input is missing.** ⛔ Never emit a conclusion the computation does not earn. **A refusal is a
   valid and valuable output**, and several tasks below may honestly end that way.
3. **Do not assume the answer to one task while doing another.** Several are deliberately able to
   contradict each other.
4. ⛔ **Do not silently choose a closure, a branch, or a division between channels.** Where a choice is
   required and not determined by this specification, introduce a **free symbol**, say so, and report the
   answer as a function of it.

---

**W1 · Projection identity.**
Integrate the four-dimensional continuity equation against `Ω(w)` over `w`, integrating by parts to
isolate the term carrying `j^w`. Report the resulting source term for a **finite** integration interval
`[w₁, w₂]` and for an **infinite** one.
⇒ `S11B_PROJECTION_IDENTITY_FINITE`, `S11B_PROJECTION_IDENTITY_INFINITE`

**W2 · Parity.** ⭐ Report the interval you use and whether it is symmetric about `w = 0`.
With `Ω(w)` even, evaluate the W1 source term for `j^w(w)` (a) **even** in `w`, (b) **odd** in `w`. Report
each. State for each whether the result is **exact** or leading-order, and **on what interval** it holds —
⛔ an evenness/oddness argument does not by itself fix an integral over an asymmetric interval.
⇒ `S11B_PARITY_EVEN_JW`, `S11B_PARITY_ODD_JW`, `S11B_PARITY_INTERVAL`

**W3 · Dynamical window.**
Repeat W1 with `Ω = Ω(w; x, t)`. Enumerate **every** term present here and absent from W1.
⇒ `S11B_DYNAMIC_WINDOW_EXTRA_TERMS`

**W4 · Bulk response to oscillating faces — impermeable.**
Solve §3 for the **two** faces displaced as `ζ₊(t)`, `ζ₋(t)` with in-plane wavevector `k` and frequency
`ω`, imposing that bulk normal velocity at each face equals that face's normal velocity. Treat the bulk on
both sides. Select the admissible branch in each regime and **state the selection criterion used**.
Report the pressure-to-velocity response, its real and imaginary parts, and any effective inertial loading
per unit area, **separately for the `ζ₊ − ζ₋` and `ζ₊ + ζ₋` combinations**.
⭐⭐ **There are THREE regimes, not two** — the bulk normal wavenumber squared may be positive, negative,
**or zero**. ⛔ Report all three, including the behaviour of every reported quantity **at** the zero
(grazing) case, where some may be singular. ⚠ Omitting this third case is a known prior defect.
⇒ `S11B_IMPEDANCE_IMPERMEABLE`, `S11B_IMPEDANCE_BY_REGIME`, `S11B_ADDED_MASS`,
  `S11B_GRAZING_BEHAVIOUR`, `S11B_IMPEDANCE_BY_PARITY`

**W5 · Bulk response — permeable face.** ⭐ Restructured; rev 1 was not closable.
An acoustic half-space is closed by one condition at the face plus a radiation/boundedness condition at
infinity. Allowing material to cross the face therefore requires an **additional constitutive law**, which
this specification does **not** supply.
⇒ **Write down the most general LINEAR, LOCAL relation** between the normal mass flux across a face and
the interface perturbation fields available (face velocity, pressure perturbation, and any thickness
perturbation). Introduce its coefficient(s) as **free symbols**, state their dimensions, and report the
modified response **as a function of them**. Then report, **for each of W4's three regimes**, whether a
dissipative (real, in-phase-with-velocity) part is present, and **on which coefficients that depends**.
⛔ Do not select a value for any coefficient. ⛔ Do not report a single "the" permeable impedance.
⇒ `S11B_PERMEABLE_CLOSURE_FORM`, `S11B_PERMEABLE_COEFFS`, `S11B_IMPEDANCE_PERMEABLE`,
  `S11B_PERMEABLE_DISSIPATIVE_BY_REGIME`

**W6 · Stress balance.** ⭐ New; rev 1 left this implicit.
State and derive the **stress matching condition** at a face: the bulk pressure perturbation against the
brane-side normal stress. ⚠ The brane moduli of §2 are **slab-integrated** while the bulk pressure is a
local 4D stress, so report explicitly what converts between them. This condition supplies `F_face` in §4.
⇒ `S11B_STRESS_BALANCE`, `S11B_F_FACE`

**W7 · The coupled system.**
Assemble §2, §4, §5's constraint, and W4/W6 into one linear system. ⭐ State which of W4 (impermeable) and
W5 (permeable) you take as the physical assembly and which as a comparison, and **why**. Report:
- the dispersion relation;
- whether the longitudinal speed inside the bulk branch condition is the same quantity the system solves
  for, and if so the resulting self-consistency condition;
- whether a self-consistent solution exists and whether it is unique;
- the long-wavelength behaviour **as a function of the assumed scaling of `k_W`** — if `k_W ∝ k^p`, report
  the exponent in `ω(k)` as a formula in `p`, in **all three** regimes of W4.
⇒ `S11B_COUPLED_DISPERSION`, `S11B_ASSEMBLY_CHOICE`, `S11B_SELF_CONSISTENCY`, `S11B_FIXED_POINT`,
  `S11B_LONGWAVE_EXPONENT_MAP`

**W8 · The transverse channel.**
Using §5's constraint as the only available coupling route, compute the coupling coefficient between the
**transverse** brane mode of §2 and `δW`, in a homogeneous brane, at linear order. ⛔ **Compute it from
the constraint — do not infer it from a divergence-free condition.** Report the coefficient. Then report
whether the transverse mode acquires any dissipative term from W5, and **on which of W5's coefficients
that depends**.
⇒ `S11B_TRANSVERSE_COUPLING`, `S11B_TRANSVERSE_DISSIPATIVE`

**W9 · Dimensions.**
Derive, from the equations of this specification and ⛔ **not** from any external table or registry, the
`[L, T, M]` exponents of: the W4 response, the W4 inertial loading, `μ_W`, `k_W`, `B_ρ`, `ρ_m`, W5's
coefficients, and the W1 source term. Show the route for each.
⚠ **A check reducing to an identity of the form `(X − 2Y) + 2Y == X` is worthless.** Label each route
**independent** or **definitional**.
⇒ `S11B_DIM_<name>` per entry, plus `S11B_DIM_ROUTE_KIND` per entry

**W10 · ⭐⭐ Reconciling `B_comp⁰`.**
With the thickness restored as a degree of freedom, derive the **effective** in-plane compression modulus
as a function of `ω` and `k`, from §5's constraint and the two energy channels of §5.
Report its limits as `ω → 0` and `ω → ∞`, and state **which limit, if either, corresponds to a
calculation in which the thickness is not a degree of freedom** — i.e. where `B_comp⁰` of §2 sits.
⛔ Do not assume the answer; ⛔ if the two limits do not bracket a single consistent identification, report
that as the result.
⇒ `S11B_B_EFF_OF_OMEGA`, `S11B_B_EFF_LIMITS`, `S11B_B_COMP0_IDENTIFICATION`

**W11 · Controls.** ⛔ All are FORM controls; ⛔ do not substitute a coefficient rescaling.
- **A** — give the bulk a nonzero shear modulus; recompute W4 and W8. Report what changes.
- **B** — make `Ω(w)` **asymmetric** in `w`; recompute W2. ⭐ Then, **separately**, keep `Ω` even and make
  the **interval** asymmetric. Report which of the two changes moves W2's result, so that a parity
  selection rule is distinguished from a domain-symmetry artifact.
- **C** — remove §5's constraint (set the coupling to zero by hand) and recompute W7. Report what
  collapses. ⚠ This exists because a system with no coupling returns a clean-looking dispersion.
For each, report which reported quantities move and which do not.
⇒ `S11B_FORM_CONTROL_BULK_SHEAR`, `S11B_FORM_CONTROL_WINDOW_PARITY`,
  `S11B_FORM_CONTROL_INTERVAL_SYMMETRY`, `S11B_FORM_CONTROL_NO_CONSTRAINT`

**W12 · Validity condition.** Emit §7's derived condition as a tag.
⇒ `S11B_BACKGROUND_VALIDITY_CONDITION`

---

## OUTPUT FORMAT

Print one line per tag, `TAG: value`. Values must be explicit expressions, not prose, wherever the value is
mathematical. End with a single `VERDICT:` line reporting whether the script's own internal consistency
checks contradicted each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ It is not a verdict
on the physics and must not be worded as one.
