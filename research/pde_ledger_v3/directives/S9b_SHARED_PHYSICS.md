# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v9)

**Authors:** Codex (v5–v6, v9), with v7–v8 edits by Claude (orchestrator). **Status:** v9 authored for review;
v8 was cleared for the build on 2026-10-06. Part D is the v9 addition.

## Why

S9 states what light needs from the medium in order to exist and to stay confined. It does not state what light
needs in order to be deflected and delayed by a mass. This step states that requirement as conditions on the
local wave equation that brane light sees near a mass. It does not derive that equation from the substrate:
the nonuniform brane operator is unfinished S11c work, and the substrate comes at the knit.

## Governing object (supplied; the build cannot test it)

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
  for both polarizations, measured relative to the shear-carrying material.
- **Advection.** `V(x)` is the in-plane velocity of the material that `u` displaces. Its kinetic symbol is
  `(ω − V^i k_i)²`: the coefficient of `ω²` is `1`, the coefficient of `ω` is `−2V^i k_i`, and the `kk`
  slot contains `V^iV^j k_i k_j`.
- **Steady brane mass balance.** The supplied balance is

  ```
  ∇·(ρ_br(x)V(x)) = −j_n(x) .
  ```

  `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
  owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object
  (`directives/S11b_SHARED_PHYSICS.md:99–111`); the relation between `j_n` and `v_dr` is open.
- **Embedding.** `ξ_w(x)` is the brane's displacement into the bulk direction, as a length. The ledger's `h`
  is dimensionless, with `ξ_w = ℓh` (`directives/S11b_SHARED_PHYSICS.md:95–96`). `g^{ij}` is the inverse of
  the induced spatial metric `g_ij` above.
- **Anchor.** With `V = 0`, `ξ_w = 0` and constant coefficients, the relation reduces to S11b's uniform
  transverse branch, `ρ_br⁰ω² = μ_⊥k²` (`steps/S11b_interface_coupling_law.md:74–87`).
- **Background anchoring of `c_γ`.** The steady object in this step is the LAB_HELD profile `c_γ(x)`
  (`directives/S11c_a_SHARED_PHYSICS.md:232–242`). The S11c-a MATERIAL_ADVECTED anchoring is not used:
  for a steady radial profile, `D_t c_γ = 0` and `∂_t c_γ = 0` give `V_r ∂_r c_γ = 0`, inconsistent with
  keeping both a radial drain and a radially varying `c_γ` live.

**Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
- direction-dependent (radial versus tangential) stiffness;
- coupling to thickness or bulk fields;
- polarization-dependent propagation;
- mixed `ωk` content other than advection.

**Reference speed.** `c₀` is the asymptotic value of `c_γ`, identified with the measured light speed. The bulk
sound speed `c_s` is separate and is not identified with `c₀`. Local light speed is not an observable here,
because rulers and clocks are made of the same medium. Only the far-field observables below are compared.

**Order.**
- **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
  coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
  terms of the underlying operator, which affect amplitude and polarization transport.
- **Smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained set is every monomial

  ```
  δ^a (V/c₀)^b ((∂ξ_w)²)^c,      0 ≤ a ≤ 1,  0 ≤ b ≤ 2,  0 ≤ c ≤ 1,
  ```

  with nonnegative integer `a`, `b` and `c`. Each retained monomial is printed separately, and nothing outside
  this set is computed.
- **Order counting** (supplied; the build cannot test it; owned by the gravity sector or S12). In the far zone,
  let `ε(r) ≡ GM/(c₀² r)`, with `GM` as supplied below. The supplied counting is

  ```
  δ = O(ε),      (∂ξ_w)² = O(ε),      V/c₀ = O(ε^{1/2})   ⇒   (V/c₀)² = O(ε) .
  ```

  Under it, every other monomial in the retained box is `o(ε)`. Within these bounds, `δ`, `V` and `ξ_w` stay
  free radial profiles in Parts A–C; Part D adds the relations and OPEN operands below. Part B compares at first
  order in `ε`. For each Part B observable, it sums the
  contributions of the grades `δ`, `V/c₀`, `(V/c₀)²` and `(∂ξ_w)²` and subtracts the reference. The `V/c₀`
  contribution is computed and printed at its own order, `ε^{1/2}`, not assumed. Every other retained monomial
  is printed as its own higher-order object and is not compared with the first-order references.
- **Profiles.** Prefer general radial functions. If an engine restricts itself to a family, it says so and keeps
  every exponent and coefficient symbolic.
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

## The observables

**Quantifier.** Each Part B comparison is required to hold for **every** far-zone `b`, not only at selected
`b`.

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

## Supplied references (inputs the build cannot test)

- **`GM`** (owned by the gravity sector). The far-field Newtonian mass parameter, as measured by the orbits of
  slowly moving test matter. It is an independent symbol, not identified with any profile amplitude. Every
  Part B and Part C condition is stated relative to it.
- **GR reference, in PPN form with `c₀`:**
  - full-flyby deflection: `Δθ = (1+γ)·2GM/(b c₀²)`;
  - for `Z_E, Z_R ≫ b`, the round-trip logarithmic term is
    `2(1+γ)(GM/c₀³)·ln(4 r_E r_R/b²)`. Part B uses only its coefficient of `ln(1/b²)`;
  - `γ = 1` in GR.
- **Part C bulk:**
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
  - For each observable, the condition on the profiles, relative to `GM`, that sets the first-order sum (see
    "Order counting") minus its reference to zero. Higher-order monomials are printed but not set against the
    references.
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

### Part D — linked steady brane (v9)

The object is the Part B and Part C conditions, relative to the same independent `GM`, and their implied
`j_n`, when the relations among `δ`, `V`, `ξ_w`, `ρ_br`, `μ_⊥` and `j_n` are included. Keep `V` and `j_n`
live. Parts A–C retain the independent-profile question; Part D states which of those conditions can also
describe a steady brane. A conditional relation keeps its domain, and an OPEN relation stays an operand in
the answer. The records do not supply a closed steady system with live flow.

**Supplied relations (unfalsifiable inputs in this build).** Source paths and line ranges, including the
applicability limits, are collected in `S9b_linked_brane_sources.md`.

- **L1 — speed, stiffness and inertia (supplied identification, pointwise).**

  ```text
  c_γ(r)² ≡ μ_⊥(r)/ρ_br(r) ,          c_γ(r) ≡ c₀[1 + δ(r)] .
  ```

  This is the governing identification above, anchored by
  `steps/S11b_interface_coupling_law.md:74–87`. The uniform dispersion in that record does not supply a
  constitutive law for `μ_⊥(r)`. Do not identify variable-coefficient basis representatives by the uniform
  coefficient fold (`directives/S11c_b_SHARED_PHYSICS.md:165–174`).
- **L2 — sourced mass conservation (supplied; live steady flow).**

  ```text
  ∇·(ρ_br V) = −j_n .
  ```

  The sourced physical evolution law in the records, with its original slab notation and measures, is

  ```text
  Σ ≡ ρ_4D W ,
  ∂_tΣ + ∇_x·(Σ v) = −(J₊ + J₋)                         (flat faces),
  ∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α             (tilted faces),
  v ≡ ∂_t u ,       J_s^α ≡ ρ_m(v_bulk,s − v_face,s^α)·n̂_s^α .
  ```

  Sources: `directives/S11c_a_SHARED_PHYSICS.md:126–148,348–371`. The steady balance is already supplied
  above; the record's `v=∂_t u` refers to its displacement dynamics and is not an instruction to identify
  S9b's steady `V` with a light-wave velocity. `J_s^α` is outward relative flux per true face area and `a_s^α`
  supplies the projected-area conversion. The identification of these finite-slab quantities with S9b's
  sheet `j_n` beyond the supplied balance remains OPEN (O6). The virtual constraint `δ_vΣ_mat=0` is not a
  second steady conservation equation.
- **L3 — embedding field identity and reduction (supplied; postulated G0 sector).**

  ```text
  ξ_w = ℓh ,       h = P₀H ≡ N₀⁻¹ ∫ dw 2f₀H ,
  f₀(w) = 1/[ℓ cosh²(w/ℓ)] ,       N₀ = ∫ dw 2f₀² ,
  M_h = N₀M₄ ,       K_h = N₀K₄ = M_h c_E² .
  ```

  Sources: v2 `ledger_stage030_electric_scalar_localized_h_closure.md:46–51,63–66,88–112,127–130` and
  `ledger_stage031_puncture_deflection_field_identity_source.md:60–76`. Here `ξ_w` is a length, `h` is the
  dimensionless reduced electric scalar, and `ℓ` is the fixed reduction scale. It is not S11b's independent
  face-centre perturbation `ζ_c` (`directives/S11b_SHARED_PHYSICS.md:95–96`). `c_E` is not identified with
  `c_γ` or `c_s` (stage 030:259–261).
- **L4 — reduced embedding-sector action (supplied, conditional static closure).**

  ```text
  A_eff = ρ_br + C_J²/κ_phase ,
  S_Lh = ∫ dt d³x [ ½ A_eff (∂_t u_L)² + ½ M_h (∂_t h)²
                   − ½ B_eff |∇u_L|² − ½ K_h |∇h|² − C_hu ∇u_L·∇h ] .
  ```

  Source: stage 030:134–142. This is the coupled scalar closure with the record's constant coefficients,
  conditional on its postulated G0 action (stage 030:15–22,228–233; `V3_STEP_PLAN.md:871–893`). Its `u_L`
  is distinct from `h`. It supplies neither a background material velocity nor a live-flow momentum
  equation. Keep it as a qualified static-sector input; promotion of its coefficients or kinetic terms to
  a draining background is O4, not a supplied operation.
- **L5 — mouth source governing the static embedding (supplied, conditional on that closure).**

  ```text
  (δΩ/δh)_mouth = η_i(k_m h − g_χh s_i) ,
  k_m = K_m/ℓ² ,       g_χh = J_m/ℓ ,       Q_χ[r_Σ,s_i] = s_i .
  ```

  Sources: stage 031:68–76,96–138,147–155. `s_i` is the puncture's orientation, distinct from Part C's
  response exponent `s`; `J_m` is a mouth coupling, distinct from the normal exchange `j_n`. The orientation
  projection uses the record's postulated frozen sleeve/profile class (stage 031:27–30,243–248). Its source
  and boundary data are not an equation for a mass-driven draining background (O4–O5).
- **L6 — exterior embedding equation (supplied, static and conditional).**

  ```text
  d/dr (r² dh/dr) = 0       (source-free exterior, constant exterior stiffness κ_ext > 0),
  h(a) = h_A ,       h_A ≡ ξ_w|_A/ℓ = P₀H|_A ,       h → 0 at infinity .
  ```

  Source: stage 031:157–174. `κ_ext` names that record's generic exterior `κ`; it is not its response
  `κ=D/B_eff` and is not `μ_⊥`. Supply the governing equation and the held datum, not its solved profile.
  This static one-puncture exterior is the Q2 regime (`V3_STEP_PLAN.md:896–915`); applying it with live
  `V`, `j_n` and varying material coefficients is O4. The datum `h_A` is not a physical holder (O5).

**Reported relations that cannot be adopted as additional live-flow equations.** The finite-slab density
factorization `ρ_br=ρ_4D W` is supplied kinematics, not a density response to `f`. S11c-a's two representatives
set respectively `ρ_4D,bg⁰=ρ_4D,ref⁰` or `ρ_br,bg⁰=rho_br` while `W_bg` varies
(`directives/S11c_a_SHARED_PHYSICS.md:210–230`); neither freeze is imposed here. Its anchoring maps are
`Q_bg^L(x,t)=Q_bg(x)` and `Q_bg^M(x,t)=Q_bg(χ(x,t))` (:232–244); Part D retains v8's LAB_HELD choice.
The face law `J_s=Λ_A(ω)𝒜_s+Λ_V(ω)V_s`, `𝒜_s=μ_θ/ρ_br⁰−δp_s/ρ_m`, is a linear perturbation response
(`directives/S11b_SHARED_PHYSICS.md:194–224`), not a supplied DC law for `j_n`.

The S11b momentum balance at `directives/S11b_SHARED_PHYSICS.md:353` is for the wave displacement `u`;
the background drain is excluded from its operators (:99–111). The S11c support test compares a stationary
energy/geometry operand with `𝒮_hold⁰={f_hold⁰,t_hold,s⁰}` while
`V_s⁰=J_s⁰=𝒜_s⁰=0` (`directives/S11c_a_SHARED_PHYSICS.md:246–279`;
`directives/S11c_b_SHARED_PHYSICS.md:190–196,215–235`). These supply no balance with live flow. The S11c
step records retain rest-frame and supplied-profile limits; their current closeout supplies no holder
(`steps/S11c_PARTIAL_CLOSEOUT.md:19–21,35–37`). Carry the support balance as O2 and the holder as O5.

**OPEN premises (named operands, not invented laws).** All dependences and spatial derivatives remain live.

- **O1 — stiffness response:** `ℳ_⊥` denotes the unknown constitutive response of `μ_⊥` to brane density,
  embedding, flow and bulk state. L1 supplies its ratio to inertia, not this response. Part C's prescribed
  speed responses must be imposed together with L1; their constitutive compatibility is conditional on
  `ℳ_⊥`, as Part C is conditional on `s`.
- **O2 — live steady momentum/support balance:** `ℬ_hold^live` denotes the missing relation among the
  profiles, body force `F_drive(r)` and face/support tractions `T_hold,s(r)`. The driving force and its
  coupling to the independent `GM` are not supplied. Neither the wave balance nor the frozen S11c support
  residual determines this operand.
- **O3 — exchange momentum:** `Π_n(r)` denotes the unspecified momentum transferred with `j_n`, including
  its direction and the velocity of exchanged material. It is an operand of `ℬ_hold^live`. The absence of a
  direct generalized `J_s` force in S11b's linear perturbation model
  (`directives/S11b_SHARED_PHYSICS.md:356–363`) does not set `Π_n` to zero in the steady flowing background.
- **O4 — applicability of the embedding closure to the live steady mass:** `ℰ_h^live` denotes the missing
  embedding/longitudinal relation including flow, exchange, varying coefficients and their gradients.
  L3 identifies the field; L4–L6 retain their postulated/static qualifications. The records do not supply
  a promotion of those equations to the current background, a map from their static `u_L` to steady `V`,
  or an identification of their mouth source with the mass's drain. Print a static-sector conditional
  restriction separately if used; do not impose it on the live system without this operand.
- **O5 — physical core holder and boundary data:** `ℋ_core` denotes the unknown holder response selecting
  `h_A` or a mouth flux/Robin datum, with its relation to `GM`. Q2 leaves the holder and amplitude as debts;
  the recorded candidate mechanisms are unselected (`V3_STEP_PLAN.md:908–939`). Keep the boundary datum
  live; no nonzero amplitude or holder mechanism is supplied.
- **O6 — normal-exchange and bulk-drain map:** `𝒥_map` denotes the missing identification of the sheet
  `j_n` with finite-slab face fluxes, their measures, and the background bulk-normal `v_dr`, together with
  the bulk/return boundary data. L2 fixes only the supplied conservation balance and its source convention.
  The S11b/S11c rest-frame response supplies no live-flow exchange law.
- **O7 — brane density response and further grades:** `ℛ_br` denotes the missing response of `ρ_br` to
  Part C's bulk density change `f`, with any thickness/projection dependence left symbolic. No selected
  density representative or bulk EOS closes it. The individual density, stiffness, source, support and
  holder grades and their derivative scales, beyond the ratio and slope counting below, are also OPEN.

**Order counting for the links.** Use exactly v8's `ε=GM/(c₀²r)` and retained monomial box. L1 enters the
`δ` grade, `O(ε)`; only the ratio `μ_⊥/ρ_br` is so constrained. Do not assign separate `O(ε)` fractional
changes to `μ_⊥` and `ρ_br`, or freeze either in a derivative. L2 contains one explicit `V/c₀`, at velocity
grade `b=1`, `O(ε^{1/2})`, with the full density factor and its derivatives retained. Its induced-metric
correction keeps v8's relative `O(ε)` qualification for `j_n`. L3's fixed `ℓ` maps the live `h` slope into
`∂ξ_w`; its square enters at `O(ε)`. This assigns no order to the mouth amplitude itself. L4–L6 enter only
as qualified static governing inputs, or through O4–O5 in the live system; no grade for their sources,
coefficients, longitudinal field or holder is supplied. O1–O7 must remain named operands wherever their
grades or forms are missing. Do not discard an operand by assigning it a convenient higher order. Part C
still keeps first order in `f`; a relation between `f` and `ε` is not supplied.

**Printed Part D deliverables.** For each Part B observable and each of the three Part C responses, print
the condition on the linked profiles relative to `GM`, its residual against the same reference, and its
implied `j_n` through L2, with `ρ_br` live. Include L1 and every supplied relation applicable on the stated
domain. Show the dependence on O1–O7 as live symbols or responses, and distinguish a condition for the
static postulated sector from one for the live steady brane. Print compatibility conditions with each
Part C response and with v8's order counting, without selecting a constitutive response, support,
exchange momentum or holder. If an OPEN operand prevents elimination of a profile or establishment of
compatibility, print the remaining conditional relation and name that operand; do not replace it by a
freeze or turn it into a solved profile. No expected value, sign, cancellation or holding mechanism is
supplied as an outcome.

## Engines, review, scope

- **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).
- **Build review.** Codex-written, so a fresh Claude agent and Grok, each with a mandatory FORM ablation.
- **Model point.** As in "Setting": leading eikonal with the retained multigraded set above. Results do not
  transfer to:
  - the strong field;
  - the throat mouth or interior;
  - a moving or rotating mass;
  - a drain that changes while light crosses;
  - polarization transport;
  - anisotropic or coupled branches.
- **Deferred to the build** (implementation, not new physics; the build directive owns it):
  - the symbolic handling of the every-`b` requirement;
  - the mass balance on the induced metric. Under the supplied counting, its difference from the flat form is a
    relative `O(ε)` correction to the implied `j_n`.
- **Stop and report**, without choosing, when any of these happens:
  - a second method failure;
  - a sub-problem this spec does not name;
  - a premise this spec does not supply or explicitly carry as an OPEN operand in Part D.
- **The step record** interprets the results.
