# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)

**Version:** v10, repair 2 (2026-10-08).

**Objects:** the v8 optical objects in Parts A–C and the neutral, linked steady in-plane
brane balance and its conditional optical requirements in Part D, retaining every unsupplied O2 input.

**Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions
are `directives/S9b_repair_decision_list.md` at `79055918` (amendment 1). v9's Part D is not an input. Paths
below are relative to `research/pde_ledger_v3/`. The additional source abbreviations are:

- **O2-R:** `steps/O2_steady_brane_balance.md`, accepted at `72866fcf`, especially §§4 and 8.
- **O2-S:** `directives/O2_SHARED_PHYSICS.md` at `4680e251`.
- **O2-C:** `directives/O2_input_contract.md` at `217a92e9`.
- **O2-P:** `directives/O2_premise_decision_list.md` at `77d2c39a`.

Equations labelled **supplied** are conditional inputs that the computation cannot test. Adopted
premises retain that status and their provenance. Every
dependent result is flagged with the supplied identification or adopted premise it uses. Reference
objects are comparison inputs, rather than premises for deriving the brane's observables. The spec
names computed objects and solution conditions; it supplies no outcome or expected-value acceptance test.

## Why

S9 states what light needs from the medium in order to exist and to stay confined. It does not state what light
needs in order to be deflected and delayed by a mass. This step states that requirement as conditions on the
local wave equation that brane light sees near a mass. It does not derive that equation from the substrate:
the nonuniform brane operator is unfinished S11c work, and the substrate comes at the knit.

## Governing object (supplied; the computation cannot test the input)

Light is the transverse branch of the brane's in-plane displacement `u`, and `u` is the material displacement
of the stuff whose density is `ρ_br` (`steps/S11_stray_longitudinal.md:32–35`; S10). Near the mass, at
leading eikonal order, that branch obeys the supplied local dispersion relation

```
(ω − V^i k_i)² = c_γ(x)² · g^{ij}(x) k_i k_j ,      g_ij = δ_ij + ∂_iξ_w ∂_jξ_w .
```

Each piece is a supplied identification. Flag any result that depends on one.
- **Local light speed `c_γ(x)`.**

  ```
  c_γ(x)² ≡ μ_⊥(x)/ρ_br(x) .
  ```

  Here `μ_⊥` is the basis-invariant transverse stiffness (`steps/S11b_interface_coupling_law.md:74–87`);
  S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
  remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
  for both polarizations, measured relative to the shear-carrying material; the supplied polarization
  identification is `c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`.
- **Advection (supplied).** `V(x)` is the in-plane velocity of the material that `u` displaces. Its
  supplied kinetic symbol is

  ```
  K_adv(ω,k;x) ≡ (ω − V^i(x) k_i)² .
  ```
- **Steady brane mass balance.** The supplied balance is

  ```
  ∇·(ρ_br(x)V(x)) = −j_n(x) .
  ```

  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
  coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as
  `ρ_br`. This input supplies no induced-measure or finite-slab replacement mass law.
  S9b computes and claims on the coordinate `d³x` measure only. No induced-measure object,
  re-expression or correction is computed or printed. The induced-measure readings and
  O2-R L478–480's recorded qualification for those readings are outside this step's
  computation and claims (D4).

  `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
  owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object
  (`directives/S11b_SHARED_PHYSICS.md:99–111`); the relation between `j_n` and `v_dr` is open.
- **Embedding (supplied).** `ξ_w(x)` is the brane's displacement into the bulk direction, as a length.
  The ledger's `h` is dimensionless. The supplied equations are

  ```
  w = ξ_w(x) ,      ξ_w = ℓh ,      g^{ij} = (g_ij)⁻¹ .
  ```

  `ℓ` is the fixed reduction scale (`directives/S11b_SHARED_PHYSICS.md:95–96`; O2-C §5), rather than
  a selected slab width.
- **Uniform anchor (supplied on its historical domain).** S11b's uniform transverse input is

  ```
  ρ_br⁰ω² = μ_⊥k² ,      V = 0 ,      ξ_w = 0 ,      ∂_iμ_⊥ = ∂_iρ_br = 0 .
  ```

  Source: `steps/S11b_interface_coupling_law.md:74–87`. It supplies no live flowing material law.
- **Background anchoring of `c_γ` (supplied).** This step selects LAB_HELD from the two distinct
  physical anchorings (`directives/S11c_a_SHARED_PHYSICS.md:232–242`; O2-C §6):

  ```
  Q_bg^L(x,t) ≡ Q_bg(x) ,      Q_bg^M(x,t) ≡ Q_bg(χ(x,t)) ,
  c_γ^L(x,t) ≡ c_γ(x) .
  ```

  `χ` is the inverse material map. MATERIAL_ADVECTED is not selected. LAB_HELD supplies neither a
  material-reference evolution law nor a physical holder.

**Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
- direction-dependent (radial versus tangential) stiffness;
- coupling to thickness or bulk fields;
- polarization-dependent propagation;
- mixed `ωk` content other than advection.

**Reference speed (supplied).** `c₀ ≡ lim_{r→∞} c_γ(r)`, identified with the measured light speed. The bulk
sound speed `c_s` is separate and is not identified with `c₀`. Local light speed is not an observable here,
because rulers and clocks are made of the same medium. Only the far-field observables below are compared.

**Order (supplied optical retained set and counting).**
- **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
  coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
  terms of the underlying operator, which affect amplitude and polarization transport.
- **Optical smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained optical set is every monomial

  ```
  δ^a (V/c₀)^b ((∂ξ_w)²)^c,      0 ≤ a ≤ 1,  0 ≤ b ≤ 2,  0 ≤ c ≤ 1,
  ```

  with nonnegative integer `a`, `b` and `c`. Each retained optical monomial is printed separately;
  optical terms outside this set are not computed. This is not a truncation of Part D's mechanical
  or inherited energy accounting.
- **Order counting** (supplied; the computation cannot test it; owned by the gravity sector or S12). In the far zone,
  let `ε(r) ≡ GM/(c₀² r)`, with `GM` as supplied below. The supplied counting is

  ```
  δ = O(ε),      (∂ξ_w)² = O(ε),      V/c₀ = O(ε^{1/2})   ⇒   (V/c₀)² = O(ε) .
  ```

  Within these bounds, `δ`, `V` and `ξ_w` stay free radial profiles in Parts A–C. Part B's comparison
  object is the sum of contributions of the grades `δ`, `V/c₀`, `(V/c₀)²` and `(∂ξ_w)²`, minus the
  reference. The `V/c₀` contribution is computed and printed at its supplied grade `ε^{1/2}`, with
  no outcome assumed. Every other retained monomial is printed separately with its order under the
  supplied counting, and is not compared with the first-order references.
- **Profiles.** Prefer general radial functions in Parts A–C. If an engine restricts those objects to a
  family, it says so and keeps every exponent and coefficient symbolic. Part D's unsupplied responses
  remain general live unknowns, including admissible gradients and material history; no engine-chosen
  response family or derivative/history cutoff is permitted (D3; O2-R §8).
- **Inherited limits and freezes.** S9 took the sharp zero-width-sheet limit, no dissipation,
  frequency-independent moduli, the continuum limit and amplitude `→ 0`
  (`steps/S9_light_requires_shear.md:349–350`). It also took two distinct background limits
  (`directives/S9_wl_rebuild_directive.md:359`; `steps/S10_two_transverse_photons.md:816–819`):
  - S9's in-plane background flow `v₀ → 0` removed all convective terms. This step lifts that freeze by
    keeping the in-plane material velocity `V` live.
  - S9's background-strain freeze `strain → 0` is named separately. This step lifts its isotropic part by
    keeping it live in `c_γ(x)`; its direction-dependent part is the anisotropic stiffness excluded above.

  S11b separately froze its normal background drain `v_dr`: the correction recorded there as
  `O(v₀|q_n|/ω)` is uncarried (`steps/S11b_interface_coupling_law.md:158–164`). That correction remains
  outside this step, the `v_dr` freeze is open, and `v_dr` is distinct from both the live in-plane `V` and
  the live normal-exchange profile `j_n`. The effect of lifting every other inherited limit remains outside
  this step's claim.

  Adopted P1 below revisits S9's no-dissipation and frequency-independent-moduli limits for the
  steady-load material/reference response (O2-P P1; O2-C §4). Those historical limits do not constrain
  Part D's relaxation response. Parts A–C retain their supplied optical dispersion and idealizations;
  compatibility of that optical regime with the material's relaxation remains a later question.

## The observables

**Quantifier.** The Part B conditions to compute are the profile conditions for matching the reference
for **every** far-zone `b`. Part D uses that same quantifier. No matching outcome is supplied.

**Setting.** One isolated, spherically symmetric mass at rest, with its drain flowing: steady, not frozen,
with `V` live. Far field; linear waves. Time is the lab time of the brane's far-field rest frame.
1. **Δθ(b):** the total turning angle of a full flyby with impact parameter `b`.
2. **Round-trip (radar) excess time.** An emitter at distance `Z_E` on one side of the mass and a reflector
   at `Z_R` on the other, along the line. Subtract the flat round-trip time. Here
   `r_E = √(b² + Z_E²)` and `r_R = √(b² + Z_R²)`. Part A prints the full excess time; the Part B comparison
   uses only the coefficient of `ln(1/b²)` in the regime `Z_E, Z_R ≫ b`, for every profile, including tails
   other than `1/r`. The radar claim is limited to this logarithmic component.
3. **The two one-way excess times** between the same endpoints, and their **nonreciprocal part** (half the
   difference). Print whether that part depends on the path or only on the endpoints.

**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them.

## Reference and bulk input objects

- **`GM`** (owned by the gravity sector). The far-field Newtonian mass parameter, as measured by the orbits of
  slowly moving test matter. It is an independent symbol, not identified with any profile amplitude. Every
  Part B and Part C condition is stated relative to it.
- **GR reference objects (supplied comparison inputs), in PPN form with `c₀`:**

  ```
  Δθ_ref(b;γ) ≡ (1+γ)·2GM/(b c₀²) ,
  Δt_RT,log,ref(b;γ) ≡ 2(1+γ)(GM/c₀³)·ln(4 r_E r_R/b²) ,
  γ_GR ≡ 1 .
  ```

  The logarithmic reference has domain `Z_E, Z_R ≫ b`; Part B uses only its coefficient of
  `ln(1/b²)`. These equations define the oracle objects to compare against. They do not specify
  either observable computed from the supplied brane dispersion, or an acceptance test.
  The reference family defines each effective `γ`; GR residuals and matching conditions use its
  supplied `γ_GR` member.
- **Part C bulk (supplied):**
  - `P = Kρ^n`, with `n` symbolic, `ρ` the bulk number density and `m` the particle mass (v3 convention). So
    `c_s² = nKρ^(n−1)/m`.
  - Define the fractional bulk-density change at the brane by
    `f(x) ≡ δρ(x)/ρ₀ = ρ(x)/ρ₀ − 1`. It is a symbolic radial profile owned by the gravity sector or S12.
    `ρ₀` is the asymptotic bulk number density, `c_s0` is the corresponding asymptotic sound speed, and Part C
    keeps first order in `f`.

## Deliverables (printed by each engine; no conclusions in the scripts)

- **Part A.** The branch-existence and path-traversal conditions. Then, for the LAB_HELD profile, as expressions
  or functionals of `δ`, `V` and `ξ_w`, with each retained monomial printed separately:
  - `Δθ`;
  - the round-trip excess time;
  - the two one-way excess times and their nonreciprocal part.
- **Part B.**
  - An effective `γ` from `Δθ`, and one from the coefficient of `ln(1/b²)` in the round-trip excess time for
    `Z_E, Z_R ≫ b`.
  - Their residuals against the references.
  - For each observable, the solution condition on the profiles, relative to `GM`, for the comparison
    residual of the selected grade sum (see "Order counting") at every far-zone `b`. This is an object
    to compute, with its domain and unresolved inputs; it is not a demanded residual value. Other
    retained monomials are printed but not set against the references.
  - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
  - The difference between the two `γ`s, printed as an object.
- **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
  to the local bulk density, with `V` and `ξ_w` left live in each:
  - `c_γ(x) ≡ c₀` (`δ(x) ≡ 0`);
  - pointwise fixed ratio: `c_γ(x)/c_s(x) = c₀/c_s0`;
  - symbolic power response: `c_γ(x)/c₀ = (ρ(x)/ρ₀)^s`, with `s` symbolic.

  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.

- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
  momentum balance below in force. Print the constructed balance pieces and conditional objects,
  their premises and remaining OPEN inputs, which of `δ`, `V`, `ρ_br` and `j_n` each condition
  determines relative to the independent orbital `GM`, which stay free, and the implied `j_n` for
  each condition. For each condition, also print what the OPEN pieces left in the in-plane balance
  must supply, each with its general live dependence (amended D3). Preserve the every-far-zone-`b`
  quantifier and print the domain of each condition.
  No profile, relation among these quantities, or success in matching is supplied.

## Part D inputs: linked steady in-plane brane

**Object and domain (D3; O2-R §8; O2-S §§1, 4–7).** The object is the brane's far-field steady
in-plane momentum balance: the in-plane component object of O2's `ℬ_hold^live`, conditional on the
adopted inputs below, with O2's named accounting and force/power qualifications retained. It is a
neutral-sector restriction of the supplied optical
setting. The supplied coordinate/radial restrictions are

```
r ≡ |x| > 0 ,      V^i(x) = V_r(r) x^i/r ,
Q(x,t) = Q(r)      (steady scalar profiles in this setting).
```

`V_r`, `ρ_br`, `j_n` and all unsupplied responses remain general live profiles/actions. Eulerian
steadiness supplies no material constancy or history cutoff. Coordinate momentum storage, transport,
exchange, force, load, energy and power densities use the same coordinate `d³x` measure as the mass
law. Native face-area factors and reductions remain explicit through O6. Coordinate `w` content and
graph-normal content are distinct. The centre graph fixes neither finite-thickness face nor interior
normal velocity content; any content it does not determine stays with `𝒩_br^live` and `𝒥_map`.

### Adopted equations and their scope

All entries in this subsection are **supplied adopted substrate inputs to a conditional model**,
rather than derived substrate laws. P1–P4 retain the label and date
**adopted substrate input to a conditional model (2026-10-06)** from O2-P. The new inputs carry the
dates and qualifications of D2. Each printed result carries the labels of every entry it uses.

**P1 — material reference (2026-10-06; O2-P P1; O2-C §4; O2-S §§2, 3.2, 6).** The carrier responds
elastically in the optical shear regime, represented here by the supplied dispersion and
`c_γ² ≡ μ_⊥/ρ_br`. It relaxes under steady load. The named equations identifying its unsupplied
relaxation and power inputs are

```
ℛ_ref/strain^live ≡ OPEN reference/strain evolution and relaxation response ,
𝒫_ref/relax^live ≡ OPEN power of that response .
```

These operand identifications supply no kernel, time scale, rate, sign or material-reference law.
The response keeps general live fields, gradients and material history. Its power is carried
explicitly, with any net supply accompanied by its physical supplier and budget. Optical
compatibility with that relaxation response remains a later light question.

**P2 — drain drive (2026-10-06; O2-P P2; O2-S §§2, 4–5).** Adopt O2's supplied representation:

```
F_drive,separate ≡ 0 .
```

The drive is dynamical order conversion in conserved material, represented through the material
stress, face/support loading and O3 exchange entries. Local source/controller functions and
mouth/collar/return/IR/bulk-boundary data are separate OPEN inventories. Their reaction and supply
systems stay named. The historical frozen-wall total-mass sink is not substituted for this drive.
The source-side S14a bridge and response-side S16 interface to `GM` remain OPEN; no `F_drive(GM)`
law is supplied.

**P3 — exchange momentum (2026-10-06; O2-P P3; O2-S §5).** Use the signed outward-loss convention
of the coordinate mass law. The supplied local carried-material equation is

```
(Π_n^carry)^i ≡ j_n V^i ,      i = 1,2,3 .
```

This fixes carried in-plane momentum only. Additional non-variational S12 partners and the system
carrying their reaction remain OPEN. Bulk-direction carry retains `𝒥_map`, `𝒩_br^live`, native
geometry and unsupplied face-to-material velocity identifications. P3 supplies no carried-total-energy
formula. `j_n`, `V_s`, native bulk velocity and `v_dr` are not identified with one another.

**P4 — bulk loading (2026-10-06; O2-P P4; O2-C §7; O2-S §§2, 3.3, 5).** Retain the postulated
shear-free scalar bulk and its supplied native face-normal traction restriction:

```
t_bulk,s^live ≡ 𝒯_bulk,n,s^live n̂_s .
```

`𝒯_bulk,n,s^live` is a general signed OPEN amplitude, and `n̂_s` is the native face normal.
There is no independent tangential bulk mechanical stress. Full `T_hold,s`, its support partition,
native geometry and reduction remain OPEN; the bulk amplitude qualifies its bulk part, rather
than being a second load beside it. A centre-graph restriction supplies no native-face reduction.
Mechanical loading, carried momentum and any external support remain separately identifiable.

**P5 — momentum density and carriage (2026-10-07; carriage 2026-10-08; amended D2).** Scoped to
Part D, this **adopted substrate input to a conditional model** supplies the brane's in-plane
momentum density and its carriage with the brane material at `V`:

```
𝒫_br,inplane^live ≡ ρ_br V ,
𝒥_br,carry^ij ≡ ρ_br V^i V^j .
```

`𝒥_br,carry^ij` names P5's contribution to the in-plane momentum current. Any other in-plane
momentum current stays a general OPEN action, counted beside this supplied current and P6's stress;
it contains neither again. Each Part D condition prints what the remaining OPEN pieces must supply.
The material current is distinct from P3's outward exchange momentum. If stressed brane material
carried additional momentum from its stress, P5 would change. Whether that enters at a retained
grade is OPEN and owned by S8. P5 supplies no total-energy law or normal-response map.

**P6 — adopted steady in-plane stress (2026-10-08, revised the same day; amended D2).** Scoped to
Part D, this **adopted substrate input to a conditional model** supplies isotropic pressure plus
linear viscous stress in the Cartesian Cauchy convention:

```
𝒯_br,inplane^live,ij ≡ T^{ij} ,
T^{ij} ≡ −p_br(ρ_br) δ^{ij}
         + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k)
         + ζ δ^{ij} ∂_kV^k ,
t_br,inplane^i ≡ T^{ij} n_j ,      F_br,inplane^i ≡ ∂_j T^{ij} ,
c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
```

`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
The shear viscosity `η` and bulk viscosity `ζ` are live general profiles, including their gradients.
P6 supplies the steady in-plane part of `𝒯_br^live` only. It is an adopted steady form, not a
consequence of P1; P1's relaxation response remains OPEN outside this form. The linear viscous
form excludes power-law creep. P6 supplies no normal stress or normal material law, relaxation
evolution law, or optical stiffness law. In the optical regime light continues to see `μ_⊥`.
The power this stress dissipates has no supplied identification with an O2 energy operand. It
stays in the OPEN energy accounting, with its physical supplier and budget OPEN. No order or value
is assigned to `c₀/c_comp`, `η` or `ζ`.

**Density link (2026-10-06; D2).** Part D alone adopts the supplied one-exponent stiffness response.
Writing the proportionality relative to the live asymptotic density and stiffness gives its input
equations, together with the supplied asymptotic optical identification:

```
ρ_br⁰ ≡ lim_{r→∞} ρ_br(r) ,      μ_⊥⁰ ≡ lim_{r→∞} μ_⊥(r) ,
μ_⊥(r)/μ_⊥⁰ ≡ [ρ_br(r)/ρ_br⁰]^α ,
c₀² ≡ μ_⊥⁰/ρ_br⁰ ,      c_γ(r)² ≡ μ_⊥(r)/ρ_br(r) .
```

`α` is one live symbol and `ρ_br⁰` stays live. The domain of the symbolic response is printed.
This supplies `μ_⊥` as a function of `ρ_br` alone in Part D, setting aside its other O1 `ℳ_⊥`
dependences there; every result using this restriction is flagged. Parts A–C keep the two local
profiles independent subject to their supplied optical ratio, and Part C's bulk-density route stands.
The link supplies no `ρ_br(f)` law or separate smallness grade for the brane density or stiffness.
The engines compute its implication for `δ`; no such implication is stated here.

**w-parity — neutral sector (2026-10-07; D2, D4).** For an electrically neutral mass the adopted
far-field `w → −w` symmetry restricts the graph displacement and the constructed centre-graph
material velocity. Define the supplied graph embedding and its material tangent lift by

```
X_graph(x) ≡ (x¹,x²,x³,ξ_w(x)) ,      U_graph ≡ V^i ∂_i X_graph .
```

`U_graph^w` denotes the bulk-coordinate component of this centre-graph material velocity, constructed
from the steady graph and `V` (O2-S L68–79; O2-R §§1–2). It is not an independent normal velocity
operand. The adopted neutral-sector equations are

```
ξ_w(r) ≡ 0 ,      U_graph^w(r) ≡ 0 .
```

Print **neutral-sector restriction** with Part D and each dependent result, including the reduction
of the supplied `g_ij`. Parts A–C keep `ξ_w` live, including the charged case. This restriction adds
no force term. Normal velocity content not determined by the centre graph remains OPEN in
`𝒩_br^live` and `𝒥_map`, including face and interior content and unsupplied face-to-material
identifications. The reduced bulk-direction carry `(Π_n^carry)^w` remains an OPEN action of those
operands with native geometry; it is not supplied as `j_n` times one centre-graph velocity
(O2-S §5). No native-face, thickness, exchange-map or normal constitutive law is supplied.

### O2 balance pieces that remain live

The engines compose the balance from these physical roles and the supplied equations above.
No assembled momentum balance, total transport tensor, cancellation or profile solution is supplied
here (D3; O2-S §§3–6; O2-R §8). The supplied content and remaining OPEN actions are:

| O2 content | Supplied content in Part D | Content retained as OPEN |
|---|---|---|
| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density and the material-carried contribution `ρ_br V^i V^j` to the in-plane current. | Any other in-plane momentum current, with general live dependence, counted beside the supplied current and P6 stress and containing neither again; normal/embedded responses and unresolved relations wherever unsupplied. The conservative antecedent is named, rather than added as another momentum species. |
| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the full steady in-plane stress in its adopted pressure-plus-linear-viscous form, with the stated Cauchy traction/force convention. | Normal material response, conservative antecedents, reference evolution outside the adopted form, and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; no additional in-plane stress duplicates P6. |
| Optical/material constitutive inputs; O1 `ℳ_⊥`, O7 `ℛ_br` | The density link supplies Part D's optical stiffness response only. | Brane-density response and every other unsupplied material or reduction identification. Parts A–C retain the v8 local profiles. |
| Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric, constructed `U_graph` and neutral centre-graph restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, face/interior normal velocity content not determined by the graph, normal response and material/order/projection map. No second normal equation is imposed. |
| Mechanical face/support loading; `T_hold,s`, `𝒯_bulk,n,s^live` | P4 supplies the bulk part's native direction only. | Full load/support partition, native normal-load amplitude, projections/maps and application-point velocities. No external support is selected by LAB_HELD. |
| Exchange; O3 `Π_n`, S12 partners/reactions | P3 supplies the coordinate in-plane carry `j_n V^i`. | Additional momentum partners/reactions and bulk-direction carry with O6/normal material content, native geometry and unsupplied face-to-material velocities. Centre-graph parity supplies no reduced bulk-carry formula. Native convective transfer and the same O3 carry are not counted twice. |
| Boundary/source/core inputs; O5 `ℋ_core`, `𝔅_A13` | P2 supplies the drain-drive representation only. | A13 branch; distinct local conversion/controller and boundary/domain inventories; physical core holder/mouth data and their response. Q2/S22 retain O5; S12 retains conversion and return partners. |

All OPEN entries are general unknown actions, with admissible live fields, gradients, entire material
history and native dependences (O2-S §3.2; O2-C §1; O2-R §8). A name supplies no closed argument list,
locality, finite internal-variable set, derivative order, stress split or constitutive family.
Each engine retains the named operands and complete live-object dependences supplied in this spec
in content no adopted premise supplies. This includes native/chart, core/material compatibility,
the derivative keys and component locations below, and general admissible dependences outside those
explicit keys. No finite dependence inventory exhausts the general OPEN actions. Checking retention
against either O2 engine's historical emissions belongs to the comparator and record, as stated below.
Historical homogeneous kinetic, uniform quadratic, static embedding and frozen-profile/face laws
retain their original domains (O2-S §§7–8; O2-C §§2–7); they supply no additional live law here.

**Paired energy content (P1; O2-S §6; O2-C §8; O2-R §8).** Carry the inherited force/power and
energy-accounting qualifications with each Part D condition wherever its material/load/exchange
content requires them. The named OPEN object is `ℬ_E^steady`, retaining `ℰ_br^live`, `𝒥_E^live`,
`𝒫_ref/relax^live`, `𝒫_convert/exchange^live`, `𝒫_boundary^live`, `𝒮_E,net` and `𝒫_E,supply`.
P5, P6, the density link and P3 do not supply the total energy, relaxation power, exchanged energy,
or net supplier/budget. The power dissipated by P6's adopted stress remains in this OPEN accounting,
with no supplied identification with an O2 energy operand or physical supplier/budget. Preserve
the unresolved `C_ref` energy-reference/improvement convention
where applicable. Pair material stress, face loads and generalized/couple actions with the actual
application-point velocities or generalized rates on compatible measures/maps. Identify shared
work occurrences without adding them twice; the channel names prescribe no additive split.
If a conditional closure requires net supply, it carries the physical supplier and budget as OPEN
inputs. A selected non-passive interface response would separately require its reservoir and stated
budget; none is selected. Part D closes no energy law or physical supplier.

**Counting (D3; O2-R §8; O2-S §7; O2-C §9).** The mechanical and inherited energy accounting is
untruncated. The optical grade box removes no material, force, exchange, boundary or power term.
Individual density/stiffness, inertia/stress/normal response, exchange/source/load, relaxation/power,
holder/mouth/embedding and derivative grades remain OPEN wherever unsupplied. Only the optical
stiffness/density ratio inherits the speed-change grade. First order in bulk `f` supplies no `f`–`ε`
relation. `j_n`, `p_br(ρ_br)`, `c_comp`, `η`, `ζ`, `α` and `ρ_br⁰` remain live. No order or value is
assigned to `c₀/c_comp`, `η` or `ζ`. No new scale is used to remove an OPEN operand or select a profile.

### Returned O2 interpretation and dependence locations

**D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7
obligation returned to the orchestrator under M1**. This is an unresolved input qualification,
distinct from ownership of a physical input. Each dependent Part D condition carries that status;
no physical interpretation or reconciliation is supplied. The dependence inputs are explicit below.

The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
`mu_perp`, `o2_rho_br_live` at order 1 and `xi_w` at order 2, each at the stored argument
`sqrt(x1²+x2²+x3²)` (O2-R §4). The names denote `V_r`, `δ`, `f`, `h`, `j_n`, `μ_⊥`, `ρ_br` and
`ξ_w` respectively. Each engine retains all eight as admissible live dependences at every following
component/role location in content no adopted premise supplies. These are complete dependence
inputs in this spec, rather than a request to inspect an O2 emission or reproduce its serialization.
Equivalent component/action notation is permitted, with the roles identifiable:

| O2 balance location | Roles carrying the eight keys |
|---|---|
| `hold_inplane[0]` | `OPEN_MomentumFlux_0_0`, `OPEN_MomentumFlux_0_1`, `OPEN_MomentumFlux_0_2` |
| `hold_inplane[1]` | `OPEN_MomentumFlux_1_0`, `OPEN_MomentumFlux_1_1`, `OPEN_MomentumFlux_1_2` |
| `hold_inplane[2]` | `OPEN_MomentumFlux_2_0`, `OPEN_MomentumFlux_2_1`, `OPEN_MomentumFlux_2_2` |
| `hold_bulk` | `OPEN_MomentumFlux_3_0`, `OPEN_MomentumFlux_3_1`, `OPEN_MomentumFlux_3_2` |
| `hold_normal` | All twelve `OPEN_MomentumFlux_i_j` roles, `i=0,1,2,3`, `j=0,1,2` |
| `energy_balance` | `OPEN_MaterialEnergyFlux_0`, `OPEN_MaterialEnergyFlux_1`, `OPEN_MaterialEnergyFlux_2` |

P5 supplies the in-plane momentum density and its material-carried current contribution; P6 supplies
the steady in-plane pressure-plus-linear-viscous stress; the density link supplies the optical
stiffness response, each with the dependence stated in its adopted equations. Other in-plane
momentum-current content stays OPEN beside the supplied current and P6 stress, containing neither
again. Their use in a constructed action is flagged; they do not replace unsupplied momentum-current
or energy-flux content. Keep the general OPEN objects, their full dependence declarations and the
labelled neutral-sector restrictions visible with their restricted component objects.
Neutral centre-graph geometry does not remove unsupplied native or material-history content.

**Additional engine-facing OPEN content.** Retain named/native geometry, chart/measure/map,
core/material compatibility and generalized/rotational work as general OPEN actions wherever they
bear on a condition. The explicit material-compatibility role is `OPEN_MaterialCompatibility` in
`energy_storage`, `energy_transport`, `energy_power` and `energy_balance`. Material energy-density
and energy-flux roles are `OPEN_MaterialEnergyDensity` and `OPEN_MaterialEnergyFlux_i`, `i=0,1,2`;
their content remains OPEN in energy storage/transport and power accounting. The named joint action
`OPEN_JointPowerAccounting` keeps force/power compatibility and the unresolved identification of
shared work/energy content visible, without choosing a decomposition or counting the same work twice.
These role names denote the existing material/energy operands above, not additional independently
additive energy species or a prescribed nesting/occurrence count.

**O2 emission provenance and comparison limits (O2-R §§4–5, 8).** Retention against historical
O2 emissions is comparator/record scope, outside engine construction from this spec. The historical
content includes every named operand and complete live-object dependence from either O2 engine,
the eight keys at every listed location, and PY-only native/chart and core/material-compatibility
content. The explicit engine inputs above retain general admissible dependences independently of
that finite inventory. The O2 difference ledger is unresolved wherever it bears on a condition.
It includes one-sided material-compatibility actions in all four named energy rows;
one-sided material density/flux actions in `energy_power`; four density/flux orientation differences
in `energy_balance`; and density/flux occurrences inside `OPEN_JointPowerAccounting` in
`energy_power` and `ℬ_E^steady` alongside separate `energy_balance` occurrences. Those differences
travel with the force/power-compatibility and duplicate-power questions. No limited inventory match
supplies full OPEN equality, compatibility or a duplicate-power conclusion. The returned interpretation
and unsupplied physical inputs remain OPEN.

**Requirement status (D5; O2-R §8).** An established sourced requirement is a condition on which the
conditional object rests: a material/phase identification or an expressly carried accounting/
admissibility obligation. Unsupplied closure forms remain OPEN. An unresolved engine difference
alone supplies no requirement.

## Model point and scope

- **Optical model point.** As in "Setting": leading eikonal with the retained multigraded set above.
  Part D's mechanical and inherited energy content retains its stated untruncated domain. Its results
  are conditional on P5's adopted momentum carriage and P6's adopted steady pressure-plus-linear-viscous
  form, with general live `η` and `ζ`. P6 is not derived from P1; power-law creep is outside this model
  point, and relaxation outside the adopted form and energy accounting remain OPEN. Results do not
  transfer to:
  - the strong field;
  - the throat mouth or interior;
  - a moving or rotating mass;
  - a drain that changes while light crosses;
  - polarization transport;
  - anisotropic or coupled branches.
