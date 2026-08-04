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

- **Standalone, print-only, no arguments, no exports.**
- **No external file reads.**
- Prefix every output tag with `WL_`.
- Strip `ConditionalExpression[0, …]` when checking that something vanishes; test poles with
  `1/expr == 0`.
- Keep total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically rather
  than raising a limit.

---

# S11b-B — SHARED PHYSICS SPECIFICATION (rev 10)

⚠ **Inserted BYTE-IDENTICALLY into both engine directives.** It is the only part they share.

⚠⚠ **Revision 10.** This step is the **homogeneous assembly** on a uniform background. A non-uniform
background together with a global plane-wave ansatz would mix wavevectors, so a single global dispersion
relation would not exist. ⇒ ⭐ **The non-uniform problem is a separate step (S11b-C).**

⭐ **The complex-frequency continuation (§1b), affinity (§2b), and corrected balance-law route (§3b) are
SUPPLIED rather than left to engine choice.** ⚠ These load-bearing inputs are binding.
**This revision is awaiting independent review; it is not presented as final or fully reviewed.**

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

`/var/projects/toy_physics/research/pde_ledger_v3/scripts/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/steps/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/paper/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/NEXT_SESSION.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/SESSION_2026-08-01.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/SESSION_REASONING.md` ·
`/var/projects/toy_physics/research/pde_ledger_v3/TECHNIQUES_THAT_WORKED.md` ·
anything matching `/var/projects/toy_physics/research/pde_ledger_v3/_review_prompt*.md` ·
**every other file in** `/var/projects/toy_physics/research/pde_ledger_v3/directives/` **except the one
assembled directive supplied to that engine** · `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/` ·
`/var/projects/toy_physics/research/pde_ledger_v3/LAUNCH_PROMPT.md` ·
`/var/projects/toy_physics/research/pde_audit/` · anything named `PREREGISTERED`/`PREREG` ·
**the other engine's deliverable**.
⛔ **The repository's entire version-control history and metadata are also barred:** the complete Git
object store, `/var/projects/toy_physics/.git/`, refs, reflogs, index, commit messages, patches and prior
file contents. Do not run `git` or any other command/API that reads any of them, including history, object,
blame, status or diff operations.
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

## ⭐⭐ 1b · COMPLEX FREQUENCY — if a root is complex, the branch fixes its continuation

⚠⚠ **This section is supplied, and it is load-bearing.** `q_out` has branch points at `ω = ±c_s0|k|` sitting
**on the real axis**. If B5 yields a complex root, ⛔ **nothing here says on which side of the real axis it
lies; that is B5's to determine, and §0 admits either.** Two continuation paths that wind differently
around a branch point reach **different sheets**, where `q_out` differs by a factor of `−1`. That exchanges
a **normal mode** for a **leaky resonance** at the *same* `ω`. ⛔ **The physical requirements below do NOT
by themselves fix the continuation.**

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

> For `Im ω > 0`, use the analytic continuation **upward from the same upper-rim values `ω + i0⁺`**. No
> dragged cut enters the upper half-plane. This supplies `q_out` there without re-imposing either real-axis
> requirement.

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

⭐⭐ **MAKE ANY DEPENDENCE MEASURABLE:** if B5 yields a complex root, report
`S11BB_IMAGINARY_PART` **also on the opposite sheet**, and report the ratio where it is defined. If there is
no complex root or the ratio is undefined, report that instead. ⛔ **No expectation is supplied for the
ratio** — report whatever the calculation gives. Any dependence on the continuation must be **visible in
the output** rather than buried in a convention.
**Wrong derivations caught:** leaving the upper-half-plane sheet undefined, re-selecting it by spatial
decay, mishandling the coalescence point, or, when a complex root exists, hiding its dependence on the
opposite sheet.
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
J_± = Λ_A(ω) 𝒜_±  +  Λ_V(ω) V_± ,
Λ_I(ω) = Λ_I⁰/(1 − iωτ_I) ,               I ∈ {A,V,X}
```

with `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰`, `τ_A`, `τ_V`, and `τ_X` real and free, and each `τ_I ≥ 0`. Keep the three
relaxation times independent throughout. The affinity and its normalization are **SUPPLIED**:

```
𝒜_± ≡ μ_s − δp_±/ρ_m ,
μ_s ≡ μ_θ/ρ_br⁰ ,       μ_θ ≡ (δU/δθ)|_{u,e_W and all other fields fixed} .
```

⛔ **Use this affinity exactly; do not derive, re-normalize, or re-select it. `μ_θ` and `μ_s` are not
interchangeable.** Derive the dimensions of every term in B7 before combining them in an equation.

### ⭐ DERIVE the face response — ⛔ it is not supplied

Two relations close the face:
- **The bulk**, as above: `δp = Z · v_bulk,±` where `v_bulk,±` is the **bulk material's** outward normal
  velocity at that face.
- **Interfacial mass balance** (kinematics, ⛔ not a result): material leaving the slab enters the bulk, so
  the bulk's outward normal velocity exceeds the face's own by the volume flux of that material:
  `v_bulk,± = V_± + J_±/ρ_m`.

⭐ **Solve these together with the closure**, keeping `τ_A` and `τ_V` independent, and report the face pressure `δp` as a function of `V_±` and
`μ_θ`. ⚠ **Report it as what it is** — with a nonzero slab-side coupling the face response is ⛔ **not** a pure
velocity-impedance, because it acquires a term driven by a **brane** variable. ⭐ Report both coefficients
separately, and ⛔ do not force the result into the shape of a single impedance if it does not have it.
⇒ `S11BB_FACE_RESPONSE`, `S11BB_FACE_RESPONSE_MU_COEFF`

### ⭐⭐ ACCEPTANCE CHECK — an independently derived value your assembly must reproduce

**With the slab-side contribution to `𝒜_±` set to zero**, first write the coefficient fixed by the supplied
affinity and define the coefficient mapping **before** solving the face equations:

```
𝒜_±|_{μ_s=0} ≡ g_p δp_± ,          Λ_p⁰ ≡ g_p Λ_A⁰ .
```

⛔ `g_p` is the sign and normalization fixed by the supplied affinity; do not fit or change it in
this reduction. With that fixed mapping, the face response **must** reduce to the following result from a
separate step:

```
Z_perm = (ρ_m r + Λ_V⁰)/(y r − Λ_p⁰) ,    r = 1 − iωτ ,    y = q_out/ω
```

Before this comparison, report the general response with `r_A ≡ 1 − iωτ_A` and
`r_V ≡ 1 − iωτ_V` still distinct. For the comparison only, apply the explicit specialization
`τ_A = τ_V = τ`, so `r_A = r_V = r`.
⭐ **Report the comparison explicitly.** ⛔ **If it does not reduce, report the disagreement — do not adjust
your derivation to match.** ⚠ A disagreement here is a real finding about one of the two steps.
⚠ **Limit of this standard:** it applies only after the stated specialization. No independently supplied
acceptance value in this directive checks the general response with `r_A` and `r_V` distinct; report that
general response, but state plainly that this acceptance check leaves that part unchecked.
⚠ **Scope:** this reduction checks the algebra combining the bulk relation, interfacial mass balance and
closure; it does **not** validate `g_p`, which is an input fixed by the supplied affinity.
**Wrong derivation caught:** a sign or algebra error that survives the stated reduction makes the
specialized face response disagree with the supplied `Z_perm` value.
⇒ `S11BB_ZPERM_REDUCTION_CHECK`

### Passivity — ⭐ TWO SEPARATE QUESTIONS. ⛔ Assert neither; compute both

⛔⛔ **THE TRANSFERRED-MASS PRESSURE WORK MUST NOT BE ABSENT.** With time-average convention
`overline(P) = ½ Re(force × velocity*)`, define positive `overline(P_out,±)` as energy leaving the slab. The
complete two-port identity is

```
overline(P_out,±)
  = ½ Re[ (δp_± + Λ_X 𝒜_±) V_±* + μ_s J_±* ]
  = ½ Re[ δp_± v_bulk,±* + 𝒜_± J_±* + Λ_X 𝒜_± V_±* ] ,
v_bulk,± = V_± + J_±/ρ_m ,                 Λ_I(ω) = Λ_I⁰/(1 − iωτ_I) .
```

The first line is the slab-side pairing; the second separates outgoing bulk transport from interfacial
conversion. ⛔ The mixed expression `½ Re(δp_±V_±* + 𝒜_±J_±*)` is incomplete by
`½ Re(δp_±J_±*/ρ_m + Λ_X𝒜_±V_±*)`. The `Λ_X⁰ = 0` comparison standard is obtained solely by making that
substitution in the displayed identity, with every other coefficient left symbolic.

1. **Port dissipativity.** On real `ω`, use the derived face response with the fixed independent input vector

   ```
   a_± ≡ (V_±, μ_s)^T ,       overline(P_out,±) ≡ ½ a_±† H_±(ω,k) a_± ,       H_± = H_±† .
   ```

   Construct `H_±`, including `Λ_X`, and compute the necessary and sufficient condition `H_± ⪰ 0` for every
   `a_±`, using its eigenvalues or principal minors. Report whether that condition depends on parameters
   alone or also on `(ω,k)` through `Z`. The supplied face and bulk laws determine this quadratic form;
   ⛔ do not emit `NOT_ESTABLISHED` for this part.
   **Wrong derivation caught:** a one-port test, or a two-port test missing either transferred-mass pressure
   work or reciprocal-traction work, can give the wrong sign for the full Hermitian power form.
2. **Thermodynamic admissibility and reciprocity.** The supplied reciprocal traction in §3b defines
   `f_X,± ≡ Λ_X(ω)𝒜_±` and the interfacial part of the power is
   `½ Re(𝒜_±J_±* + f_X,±V_±*)`. Thus the supplied mixed representation is

   ```
   ( J_±   )   ( Λ_A  Λ_V ) ( 𝒜_± )
   ( f_X,± ) = ( Λ_X   0  ) ( V_± ) .
   ```

   On real `ω`, construct its Hermitian power form for arbitrary complex `(𝒜_±,V_±)`, and report the
   necessary and sufficient condition on `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰`, `τ_A`, `τ_V`, and `τ_X` for nonnegative interfacial
   entropy production. ⚠ This is a mixed force/flux representation. Time-reversal invariance remains
   **not assumed** (§3).
   For the conditional Onsager–Casimir test, the definitions make `𝒜_±` and `f_X,±` even under time
   reversal and the rates `J_±` and `V_±` odd; this parity assignment does not assume that time reversal is
   a symmetry. Compute and report the relation between `Λ_X` and `Λ_V` **if one is forced**, and distinguish
   a conditional reciprocity relation from unconditional admissibility. Separately for admissibility and
   for conditional reciprocity, report whether any relation among `τ_A`, `τ_V`, and `τ_X` is forced.
   Finally report the necessary-and-sufficient coefficient condition **in whatever form it takes**, both
   without and, if applicable, with reciprocity. The supplied laws determine the admissibility form; ⛔ do
   not emit `NOT_ESTABLISHED` for it. ⛔ Assert no outcome in advance, and ⛔ treat no coefficient as the
   expected one.
   **Wrong derivation caught:** omitting the supplied second row, changing the displayed power pairing, or
   identifying relaxation times before computing can give a false reciprocity relation or admissibility
   condition.

⛔⛔ **DO NOT IMPOSE EITHER CONDITION TO REMOVE A GROWING ROOT.** ⭐⭐ **Instead, report whether a growing root
lies inside each region you computed, distinguishing unconditional admissibility from any smaller region
that also imposes reciprocity.** Still report every root as a first-class outcome; admissibility is a
classification, ⛔ not an acceptance gate.
⇒ `S11BB_TWO_PORT_POWER_IDENTITY`, `S11BB_PORT_DISSIPATIVITY`, `S11BB_PORT_CONDITION_KIND`,
  `S11BB_ONSAGER_CONDITION`, `S11BB_ONSAGER_RECIPROCITY`, `S11BB_ONSAGER_DETERMINABLE`,
  `S11BB_RELAXATION_TIME_RELATIONS`, `S11BB_COEFFICIENT_ADMISSIBILITY`,
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
BASIS YOURSELF AND CHECK IT.** ⚠ An omitted allowed term can change the dispersion while a control that
removes a *listed* term still reports cleanly.

⭐ **Do this as a task, before deriving anything:** enumerate the fields and their first gradients
(`u`, `∇u`, `θ`, `∇θ`, `e_W`, `∇e_W`) and construct **every** scalar quadratic in them allowed by the
symmetry group below. Compare that basis against the list above.

⛔⛔ **THE SYMMETRY GROUP, STATED IN FULL.** Isotropy and `w → −w` alone are **not enough to define the
basis**: under those alone a **pinning term `½K|u|²` is allowed**, and an engine
that includes it **gaps the modes** while one that does not finds them gapless. Both would be obeying the
text. The group is:
- **In-plane translation invariance** ⇒ `u` may enter **only through its gradients**, never undifferentiated.
- **In-plane isotropy _and_ parity** (the full `O(3)` acting on the three in-plane directions).
- **Reflection `w → −w`.**
- **Equivalence modulo total divergences** — two densities differing by a total in-plane divergence are the
  same term; ⛔ do not count both.
- ⛔ **Time-reversal is NOT assumed**, and no positivity or boundedness is assumed (§0).

⭐⭐ **JUDGE INDEPENDENCE AS FIELD BILINEARS, with `B1`'s constraint NOT applied.** ⚠ Redundancy *"modulo
B1's constraint"* has **no single meaning** — `B1` is sourced, carries memory, and changes rank at
`ω = 0`, so "modulo B1" differs between the impermeable reduction and the flux-on case.
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

⚠ **Under ORDINARY single-copy variation a retarded kernel cannot be varied.** In a bilinear
`∫ λ 𝒴[∂_t δW] dt`, variation transposes the operator, whose symbol is `Y(−ω)` — the **advanced** kernel,
moving its finite response pole from the lower to the upper half of the complex `ω` plane. ⇒ **an
irreversible flux may not be placed inside the varied functional.** ⚠⚠ **This is a statement about
ORDINARY variation only** — ⛔ doubled-variable (in-in / Galley)
constructions **do** yield genuinely retarded kernels, and ⭐ **using one as an INDEPENDENT CROSS-CHECK on
the route below is allowed and encouraged.** ⛔ Do not treat the signatures below as universal laws; they
are signatures **under ordinary single-copy variation**.

⛔⛔ **THREE FORBIDDEN ROUTES**, each wrong under ordinary variation, ⛔ not merely disfavoured:
- **(i) Substituting the mass balance WITH ITS FLUX SOURCE into `U` and then varying.** With
  `(Λ_A⁰,Λ_V⁰) ≠ (0,0)`, the
  relation carries a history functional and may not enter a stored energy. ⚠ Typical signature: an extra
  root raising the dispersion determinant's **degree by one** — an invented mode.
- **(ii) Varying `J_±` or `δp_±` inside the action.** Manufactures an anti-causal, **energy-generating**
  kernel. ⚠ Typical signature: a finite pole inherited from a response kernel lies in the upper half of
  the complex `ω` plane.
- **(iii) Any route in which a response kernel is differentiated with respect to a field.**

### ⭐⭐ THE VIRTUAL-DISPLACEMENT RULE — binding, and it is where two engines would otherwise diverge

Use `δ_v` only for a virtual variation, to distinguish it from the perturbation named `δW`. Let `X` label
an in-plane material point and define the current in-plane map, its Jacobian, the Eulerian slab density,
and the material-area mass by

```
x(X,t) = X + u(X,t) ,                         𝒥_x ≡ det(∂x/∂X) ,
u(X,t) = u(x,t) + O(2) ,
Σ_E(x,t) ≡ Σ(x,t) ≡ ρ_4D(x,t) W(x,t) ,
Σ_mat(X,t) ≡ Σ_E(x(X,t),t) 𝒥_x(X,t) .
```

**A virtual variation is taken at one instant and transfers no mass through either face.** Its binding
constraint is the material, ⛔ not Eulerian, equation

```
δ_v Σ_mat(X,t) = 0 .
```

On the uniform background, the first-order identities and their virtual variation are

```
Σ_E = ρ_br⁰(1 + θ + e_W) + O(2) ,            𝒥_x = 1 + ∇_x·u + O(2) ,
δ_v e_W = δ_v(δW)/W₀ ,
δ_v θ + δ_v e_W + ∇_x·δ_v u = 0 .            (VIRTUAL CONSTRAINT)
```

⛔⛔ **The term `∇_x·δ_v u` MUST NOT BE ABSENT.** The equation `δ_vΣ_E = 0` is not the material constraint
and must not be used. Thus `δ_vθ`, `δ_v(δW)` and `δ_vu` are not independent.

⛔ **Do NOT vary `U` with `θ` held fixed.** Impose the displayed VIRTUAL CONSTRAINT either by eliminating
one virtual variation or by a Lagrange multiplier which you then eliminate. ⭐ **Report which you used and
report the multiplier's physical identity.** ⚠ The **same** multiplier supplies the in-plane restoring
force and the thickness term; ⛔ you may not keep one and drop the other.

⭐ **USE THIS AND NO OTHER:**
1. **State variables `u`, `δW`, `θ`**, with
   `δ_v θ + δ_v e_W + ∇_x·δ_v u = 0` binding their **variations** as displayed above.
2. **Constitutive quantities are VARIATIONAL derivatives of `U`:** `μ_θ ≡ δU/δθ` and `p_W ≡ δU/δe_W`.
   ⛔⛔ **Functional, not ordinary partial** — if §3's basis carries gradient terms such as `|∇θ|²`, an
   ordinary partial drops their contribution while the variational derivative keeps a `k²` term, and two
   engines would obtain different physics. ⭐ **State, for every derivative you write, what is held fixed.**
3. **In-plane momentum balance** for `u`.
4. **Thickness equation**, obtained from the variation under rule above — ⛔ **not** from a verbal
   description of which forces "push the faces apart".
5. **The mass balance WITH its `J_±` source is a separate EVOLUTION equation**, restored **after** the
   variation. ⛔ It is not the same object as `δ_vΣ_mat = 0`.
6. **`δp_±` and the affinity traction proportional to `Λ_X𝒜_±` are PRESCRIBED mechanical responses**
   entering through the external virtual work displayed below. **`J_±` enters the mechanics through the
   separate mass evolution and through its contribution to the prescribed `δp_±`; it has no direct
   mechanical generalized-force term:**

   ```
   Q_J^direct = 0 .
   ```

   ⛔ Substitute the prescribed responses only after 3–5 are complete; ⛔ vary or differentiate them with
   respect to nothing; ⛔ add no term proportional to `J_± δ_v(δW)` or `J_± ∇·δ_v u`, and no mechanical
   face term beyond the two terms supplied below.

### ⛔⛔ GEOMETRY AND SIGN CONVENTIONS — fixed here, because a sign error here is an ENERGY-SOURCE error

The face geometry, traction and complete pressure virtual work are

```
ζ_± = ±δW/2 ,               n̂_± = ±ŵ ,
V_± ≡ (∂_tζ_±)(±1) = ½∂_t(δW) ,
δ_vx_± = n̂_± δ_v(δW)/2 ,   t_± = −(δp_± + Λ_X(ω)𝒜_±) n̂_± ,
δ_v𝒲_bulk ≡ Σ_{s=±} t_s·δ_vx_s
            = −½[δp_+ + δp_− + Λ_X(ω)(𝒜_+ + 𝒜_−)] δ_v(δW) .
```

⭐ Both faces therefore move outward under positive `δ_v(δW)`. ⛔ **Use the displayed virtual-work
equation; do not guess a pressure force and do not append another flux-force term.**
⚠⚠ **Taking this sign the other way reverses the pressure-work term and can manufacture an instability.**
The independent, explicitly bounded pressure-work check below is the diagnostic for that error.

### ⭐⭐ THE CAUSALITY DIAGNOSTIC — mechanical, gradeable, ⛔ and it does NOT foreclose instability

In the **equations of motion, closure, face response and dispersion determinant**, every response kernel
must retain the **retarded analytic structure** supplied here. Reduce algebraically equivalent forms far
enough to cancel removable factors, then locate every finite pole inherited from a response kernel and
report its position relative to the real `ω` axis. Such inherited poles must remain in the lower half-plane;
an upper-half-plane inherited pole is an advanced response. ⛔ **Judge the analytic structure, not the
literal appearance of `Λ_I(−ω)`, `Z(−ω)`, `Y(−ω)` or conjugated factors:** rationalizing a denominator may
display such a factor without changing the object's poles. Conjugated kernels in the time-averaged and
Hermitian power forms required by §2b are outside this diagnostic. Zeros of the dispersion determinant,
including growing-mode roots, are not response-kernel poles and are also outside this pole test. ⭐ **Run
this scoped check and report its outcome explicitly.**
**Wrong derivation caught:** transposing or varying a retarded kernel leaves an uncancelled inherited
response pole above the real axis in an in-scope dynamical object.
⚠ **If this pole test finds an upper-half-plane inherited response pole, re-check step 6 and report both
the pole and the outcome.** ⛔ Do not use this diagnostic to delete or reclassify a dispersion root.
⇒ `S11BB_CAUSALITY_CHECK`, `S11BB_KERNEL_POLE_LOCATIONS`

### ⭐⭐ TWO MANDATORY CONVENTION CROSS-CHECKS — ⛔ run them, ⛔ they must be able to fail

⚠ These exist because the variational convention has **two candidate readings that give different
dispersion relations**, and a wrong one manufactures growth. ⛔ Report both outcomes explicitly.

**(a) The in-plane equation your variation produces must carry the restoring force `−∇(δU/δθ)`.**
⭐ This single check selects the convention uniquely. ⛔ If your in-plane equation does not have it, your
variational rule is wrong — ⛔ fix the rule, do not patch the equation.
**Wrong derivation caught:** omitting `∇·δ_vu` from the VIRTUAL CONSTRAINT, or varying at fixed `θ`, removes
this contribution from the in-plane balance.

**(b) With NO bulk, `Λ_A⁰ = Λ_V⁰ = Λ_X⁰ = 0`, `κ_W = 0` and `k = 0`:** report `ω²` for the thickness mode, and
report **the inequality on `B_ρ⁽³⁾`, `C`, `k_W` under which it is positive.** ⭐ It must be positive for
**every** `U` that is positive-definite — ⛔ not merely on a convenient sub-region. ⭐ **State explicitly
whether `B_ρ⁽³⁾` appears in that `ω²` at all.**
The kinematic reduction used in this check is explicitly

```
J_+ = J_- = 0 , k = 0  ⇒  ∂_t(θ + e_W) = 0  ⇒  θ + e_W = constant ,
K_check ≡ d² U(θ = constant − e_W, e_W)/de_W² .
```

⭐ Report whether the thickness stiffness from the equation equals `K_check`.
**Wrong derivation caught:** a fixed-`θ` thickness variation or omission of the common constraint
multiplier gives a conservative stiffness inconsistent with the independently reduced stored energy.

⛔⛔ **SCOPE OF CHECK (b) — read this before applying it anywhere else.** It is confined to the case with
**no bulk, no permeation, no reciprocal traction, and positive-definite `U`**, where the system is
conservative and growth would therefore be a **derivation error**. ⛔ **It says NOTHING about the full
problem** and ⛔ **must not be used to reject a growing root of B5**, where the bulk, the interface and
possibly indefinite moduli are all in play. ⚠ Confusing the two would re-close the falsification channel
of §0.
⇒ `S11BB_CONVENTION_CHECK_INPLANE`, `S11BB_CONVENTION_CHECK_CONSERVATIVE`,
  `S11BB_CONSERVATIVE_POSITIVITY_INEQUALITY`

### ⭐ ENERGY ACCOUNTING — three discriminators that can fail

⚠ Enumerating loss channels **from the same equations used to derive them** and then checking power balance
against that enumeration is an identity for any system built that way — **it cannot fail.** Use these
independent discriminators:

1. ⭐ **Compute `d/dt(T + U)` from YOUR equations.** Report every signed external exchange term — sink or
   source — and name the transport process each term corresponds to. ⛔⛔ **If a term corresponds to no
   transport process, REPORT IT** as a defect in the derivation.
   **Wrong derivation caught:** a response-kernel or flux-force term not among the supplied pressure,
   reciprocal-traction and transfer terms creates an energy exchange with no transport counterpart.

2. ⭐⭐ **INDEPENDENT PRESSURE-WORK SIGN CHECK.** Restrict only this diagnostic to real `ω`, propagating
   `q² > 0`, impermeable faces, and the form cut `Λ_X⁰ = 0`. Then the equations supplied independently by
   §2 give

   ```
   J_± = 0 ,                 v_bulk,± = V_± ,
   δp_± = Z V_± ,           overline(P_bulk,±) = ½ Re(δp_± V_±*) ,
   overline(d(T+U)/dt)|_pressure = −Σ_{s=±} overline(P_bulk,s) .
   ```

   Compute the left side **off shell** by contracting the slab equations with their corresponding
   velocities and substituting their pressure-force terms while leaving the harmonic amplitude free. Do
   not impose the homogeneous thickness equation or dispersion relation, and do not first replace the
   left side by the literal period average of an exact total derivative. Compute the right side from the
   outgoing bulk flux, then report their symbolic difference without changing either derivation.
   **Wrong derivation caught:** reversing the traction or face-displacement sign changes the slab pressure
   work into a source while the independently computed outgoing bulk power keeps its sign.

⛔⛔ **SCOPE OF CHECK 2:** it checks only the pressure-work sign in the displayed real-frequency,
propagating, impermeable, zero-reciprocal-traction sub-case. It does not classify any B5 root and must not
be used to discard or re-branch a growing root of the full problem.

3. ⭐⭐ **FULL TWO-PORT BALANCE CHECK.** At real `ω`, keep `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰` and all three `τ_I`
   symbolic. This is an **off-shell coefficient check**: contract each slab equation with its corresponding
   velocity, substitute the equation's force terms, and keep the harmonic amplitudes and the face
   quantities `δp_s`, `J_s` and `𝒜_s` algebraically free. ⛔ Do not impose the homogeneous equations,
   determinant or any on-shell amplitude relation, and do not replace the result by the literal period
   average of an exact total derivative. Decompose that slab power expression by order in `Λ_X⁰`, and at
   every order compare it face by face and channel by channel with the independently supplied slab-side
   exchange

   ```
   −½ Σ_{s=±} Re[ (δp_s + Λ_X(ω)𝒜_s) V_s* + μ_s J_s* ] .
   ```

   Report every symbolic difference without changing either derivation. This check classifies no root and
   makes no statement about the sign of any imaginary part.
   ⚠ **Verified blind spot:** because the face response, closure and the `δp/ρ_m` normalization inside the
   affinity enter both sides through the same opaque face quantities, they cancel from this comparison.
   This check does not validate any of those three objects; it tests only the slab traction,
   constraint/multiplier bookkeeping and the slab-side `μ_s` conjugate normalization.
   **Wrong derivation caught:** a traction sign or factor error, an omitted face, or an elimination against
   B1 that substitutes the wrong slab-side conjugate leaves a nonzero difference in at least one channel.
⇒ `S11BB_ENERGY_SINKS`, `S11BB_ENERGY_SOURCES`, `S11BB_UNATTRIBUTED_SINK_TERMS`,
  `S11BB_UNATTRIBUTED_EXCHANGE_TERMS`,
  `S11BB_PRESSURE_WORK_SIGN_CHECK`, `S11BB_FULL_TWO_PORT_BALANCE_CHECK`

⛔⛔ **DO NOT USE THE ZERO-INTERFACE-COEFFICIENT LIMIT AS A CHECK ON THE ROUTE.** All candidate routes
**coincide** at `Λ_A⁰ = Λ_V⁰ = Λ_X⁰ = 0`, so that limit **cannot discriminate between a correct derivation
and any of the forbidden ones.** ⚠ It is also **not** the lossless limit — §2 trap 2: radiation resistance
survives it. ⛔ Do not describe this limit as "reversible".

⛔ **The §3 list supplies no single in-plane compression modulus.** In it, compression is carried by `θ` and
by `e_W`, and how they combine is task **B4**. ⚠ Where a modulus measured with the thickness held fixed
would sit is an **output**, ⛔ not an input.

⚠⚠ **BUT THAT IS A PROPERTY OF THE LIST, ⛔ NOT AN ASSUMPTION YOU MAY IMPOSE.** If §3's basis construction
yields **one or more independent invariants containing `∇·u`**, ⭐ **carry every one with its own free
coefficient as §3 instructs** — ⛔ do **not** drop any to preserve the property above — and report how each
changes B4's identification.
⚠ Both readings have been proposed; ⭐ **report which one your basis construction actually produces.**
**Wrong derivation caught:** omitting an allowed independent invariant — including carrying only one
divergence-containing invariant — leaves a mismatch in the required full-basis comparison and changes
downstream response calculations.

---

## TASKS

⛔⛔ **RULES.** (1) Every value not explicitly supplied as a route or convention must be computed, ⛔ not
asserted. (2) If a task cannot be completed, emit
`NOT_ESTABLISHED` and name the missing input — **a refusal is valid and valuable**. (3) ⛔ Never silently
choose a branch, closure, path, or expansion; introduce a free symbol and report the dependence.
(4) Except for equations explicitly marked **SUPPLIED**/**GIVEN**, route or normalization disambiguations,
and exact FORM-control cuts, ⛔⛔ **no task states the form of its answer.** If a requested object turns out
not to be a scalar — if the response is operator-valued, or the polarizations split, or the effect is a
spatial attenuation rather than a frequency shift — ⭐ **report that instead**, and say so explicitly.

**B1 · The constraint.** ⛔⛔ **The exact balance is GIVEN — it is kinematics, not a result — and your job
is to linearise it.** ⚠ The written equation makes every source term gradeable; a merely verbal statement
would allow the in-plane term to be omitted, producing a longitudinal mode with no restoring force that
**both engines could agree on**.

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
in §2b and the reciprocal traction proportional to `Λ_X`. Report both operators, isolate the `Λ_X` term in
the thickness operator, and report the symbolic difference from the `Λ_X⁰ = 0` operator.
⇒ `S11BB_INPLANE_EOM`, `S11BB_THICKNESS_EOM`, `S11BB_BULK_FORCE_ON_THICKNESS`,
  `S11BB_RECIPROCAL_TRACTION_THICKNESS_EFFECT`

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
**Wrong derivation caught:** dividing B1 by `ω` before checking `ω = 0` can erase an integration constant
and change the reported static response.
⇒ `S11BB_COMPRESSIONAL_RESPONSE`, `S11BB_LIMITS_AND_PATH`, `S11BB_FROZEN_THICKNESS_IDENTIFICATION`

**B5 · The longitudinal mode.** Assemble and report the dispersion relation. Report whether it admits a
closed-form `ω(k)`, whether roots are real, and for any complex root its imaginary part and **which
physical ingredient makes it nonzero** — ⛔ distinguishing the two mechanisms of §2 trap 2. ⭐ If a root
fails to exist as a normal mode, report that.
⭐⭐ **Report the SIGN of every imaginary part and classify each root as decaying or GROWING**, and report
the condition on the moduli and on `C` that separates the two. ⛔ **A growing root is a reportable result**
(§0) — ⛔ do not suppress it, re-branch to remove it, or assume the quadratic form is positive-definite.

⭐ Carry `Λ_X` symbolically through the determinant and roots. Compare the result with the form cut
`Λ_X⁰ = 0`, and report which roots, multiplicities, and classifications change. ⛔ Assert no outcome in
advance.

⭐⭐ **FOR ANY GROWING ROOT, REPORT ALL THREE OF THESE — it is what separates a finding from an artifact:**
1. **The causality diagnostic of §3b** — after removable factors are cancelled, where is every finite pole
   inherited from a response kernel relative to the real `ω` axis? Do not count the required conjugations
   in time-averaged or Hermitian power forms, and do not count the dispersion root itself as a kernel pole.
2. **The sheet the root sits on**, and confirmation that requirements 1–2 of §1 were ⛔ **not** re-imposed
   at complex `ω` to obtain it.
3. ⭐⭐ **Whether it lies INSIDE each thermodynamically admissible region computed in §2b**, distinguishing
   unconditional admissibility from any sub-region that additionally imposes reciprocity. ⛔ Do not use
   either as a gate that removes, suppresses or re-branches any root.
   **Wrong derivation caught:** suppressing a root with an admissibility condition, or conflating the
   unconditional and reciprocity-conditioned regions.
⇒ `S11BB_LONGITUDINAL_DISPERSION`, `S11BB_ROOTS`, `S11BB_IMAGINARY_PART`, `S11BB_DISSIPATION_ORIGIN`,
  `S11BB_ROOT_STABILITY_CLASS`, `S11BB_STABILITY_CONDITION`, `S11BB_GROWTH_ARTIFACT_DIAGNOSTICS`,
  `S11BB_RECIPROCAL_TRACTION_ROOT_EFFECT`

**B6 · The transverse mode, computed.** On this uniform background, compute the coupling between the
transverse in-plane mode and the thickness degree of freedom **from B1's constraint and §3's energy**,
⛔ not by asserting a divergence-free argument. Report the coefficient and any modification to the
transverse dispersion. ⭐ **State explicitly what the coefficient couples to what**, and its normalization,
⛔ before assigning it a value or a dimension; ⚠ if it vanishes identically, say so and say that its
normalization is then undetermined. ⭐ Report whether the transverse mode acquires an imaginary part, and
its dependence on `Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰`, `ωτ_A`, `ωτ_V`, `ωτ_X` and the slab-side part of `𝒜` across the full range.
⚠ ⛔ **Whatever this returns, it does NOT settle whether confinement is unconditional** — that is a
non-uniform question and out of scope here. ⛔ Do not phrase it as if it does.
⇒ `S11BB_TRANSVERSE_COUPLING`, `S11BB_TRANSVERSE_DISPERSION`, `S11BB_TRANSVERSE_DISSIPATION`

**B7 · Dimensions.** Derive from the equations above, ⛔ not from any table, the `[L,T,M]` exponents of
`B_ρ`, `B_ρ⁽³⁾`, `μ_W`, `k_W`, `κ_W`, `C`, B3's response, B4's response, B6's coefficient, **and the
coefficient of any additional independent invariant you carried under §3**, plus `Λ_A⁰`, `Λ_V⁰`,
`Λ_X⁰`, `τ_A`, `τ_V`, `τ_X`, `𝒜`, `μ_θ`, `μ_s` and §2b's face response. ⚠ Each of those responses is a ratio; ⛔ state what of what before assigning
a dimension, and if a coefficient vanishes identically say that its dimension is undetermined. Show each route and
label it **independent** or **definitional** — a route whose asserted equation *defines* the symbol under
test is definitional.

⭐⭐ **AND A HOMOGENEITY CHECK ON EVERY FINAL EQUATION, WITH UNITS RESTORED.** For each equation you report
— the in-plane equation, the thickness equation, the mass balance, the affinity `𝒜`, the closure, the face
response, the two-port power identity, and the dispersion determinant — ⭐ **verify that every additive term carries the same `[L,T,M]`
dimensions**, and report the outcome per equation.
⛔⛔ **THE CHECK MUST BE ABLE TO FAIL.** ⭐ **Demonstrate that it is:** corrupt one term's dimensions
deliberately, confirm the check reports a failure, and restore it. ⚠ Report that demonstration.
**Wrong derivation caught:** combining `μ_θ` directly with a chemical potential per mass, or omitting a
factor of `ρ_br⁰` or `ρ_m`, makes the affinity or two-port power identity dimensionally inhomogeneous;
assigning `Λ_X⁰` the dimensions of `Λ_A⁰` instead of deriving them makes the traction term inhomogeneous.
A check confined to the final determinant would not expose which normalization failed.
⇒ `S11BB_DIM_<name>`, `S11BB_DIM_ROUTE_KIND_<name>`, `S11BB_HOMOGENEITY_<equation>`,
  `S11BB_HOMOGENEITY_ABLATION_DEMO`

**B8 · Controls.** ⛔ FORM controls; a coefficient rescaling tests none of them.
- **A — remove the thickness channel** (hold `δW = 0`) and recompute B4 and B5.
  Its exact cuts are

  ```
  δW = 0  ⇒  e_W = 0 ,  V_± = 0 ,  Λ_V(ω)V_± = 0 ,
              δp_±V_± = 0 ,  Λ_X(ω)𝒜_±V_± = 0 .
  ```

  ⚠⚠ **A IS CONFOUNDED AND YOU MUST REPORT IT AS SUCH.** It removes the thickness field, face-motion drive,
  velocity-coupled permeability, mechanical pressure work and reciprocal-traction work together. ⛔ But
  with the slab-side affinity active, `V_± = 0` does **not** imply `J_± = 0` or `δp_± = 0`; do not delete
  that remaining transfer channel. ⇒ ⛔ **Do not attribute any change under A to the thickness degree of
  freedom alone.** ⭐ For each quantity that moves, report **which of the simultaneously-removed channels
  it could be attributed to**, and **say so plainly if the attribution cannot be separated** by this control.
  **Wrong derivation caught:** setting `J_± = δp_± = 0` merely because `V_± = 0`, which silently removes
  the slab-state-driven transfer channel from control A.
- **B — remove the gradient stiffness** (`κ_W = 0`) and recompute B3 and B5.
- **C — impermeable faces** (`Λ_A⁰ = Λ_V⁰ = 0`, with `Λ_X⁰` left symbolic) and recompute B5.
- **D — remove the cross term** (`C = 0`) and recompute B4 and B5.
- **E — remove the SLAB-SIDE part of the affinity `𝒜`** (bulk-pressure and velocity couplings untouched) and recompute B5
  and the passivity results of §2b. ⭐ This is the control that isolates the **new** transfer channel; it is
  ⛔ **not** the same cut as **C**, which removes the mass-flux closure row while leaving the reciprocal
  traction symbolic.
- **F — remove the reciprocal traction** (`Λ_X⁰ = 0`, with `Λ_A⁰` and `Λ_V⁰` left symbolic) and recompute
  B3–B5 and both passivity questions of §2b. This is a FORM cut: its recomputed mechanical operator must
  equal the symbolic `Λ_X⁰ = 0` operator reported under B2, and its power identity must equal the §2b
  identity after the same substitution; report both symbolic differences. **Wrong derivation caught:**
  dropping `Λ_X` from the full system as well as the control makes F pass without testing the added
  traction channel.
⭐⭐ **Recompute B6 under every one of A–D as well, and report what moves.** ⛔ Do not assume in advance
that a control cannot affect B6, and ⛔ do not discard a dependence you find on the grounds that it "must
be" algebraically predetermined. ⚠ If none of A–D changes B6's reported quantities, **state that as a
finding and say why** — that is a result about the structure of the coupling, and it must be **discovered
here, not assumed**.
⭐ For each, report which reported quantities move and which do not. ⛔ Report what each control does,
⛔ not what it was expected to do.
⭐ In every control keep `τ_A`, `τ_V`, and `τ_X` independent wherever the corresponding coefficient remains,
and report which, if any, becomes irrelevant because its entire channel was cut.
⇒ `S11BB_CONTROL_NO_THICKNESS`, `S11BB_CONTROL_A_ATTRIBUTION`, `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`,
  `S11BB_CONTROL_IMPERMEABLE`, `S11BB_CONTROL_NO_CROSS_TERM`, `S11BB_CONTROL_NO_MU_COUPLING`,
  `S11BB_CONTROL_NO_RECIPROCAL_TRACTION`, `S11BB_CONTROLS_ON_TRANSVERSE`

**B9 · Validity.** Report the conditions under which this step's linearisations hold, including §2's
background-flow condition, and **where in `(ω,k)` any fail**. State separately any validity condition you
use for each response kernel in terms of `ωτ_A`, `ωτ_V`, and `ωτ_X`; do not replace them by a common
`ωτ` unless you explicitly impose the equal-time specialization.
⇒ `S11BB_VALIDITY_CONDITIONS`, `S11BB_VALIDITY_FAILURE_REGION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`; explicit expressions wherever mathematical. End with a single `VERDICT:`
line reporting whether the script's own internal consistency checks contradicted each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ Not a verdict on the
physics.
