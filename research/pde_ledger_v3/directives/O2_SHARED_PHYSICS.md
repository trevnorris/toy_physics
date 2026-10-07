# O2 — shared physics for the live steady momentum and support balance

**Author:** Codex, 2026-10-06 (versions preserved at `f18c67e8`, `a7badd69`); revisions 1 and 2 by fresh
Claude authors, 2026-10-07. **Status:** revision 2 on the reviewed baseline `a7badd69`; nothing accepted,
and no review clearance or computed result is claimed.

**Deliverable:** specify the untruncated, conditional in-plane and normal momentum/support object
`ℬ_hold^live`, with its energy pairing, live profiles and every OPEN input retained, for independent
construction by SymPy and a blind Wolfram engine.

This is sub-step 3 of the O2 inventory. Its authority is the spec-authoring directive at `48647e55`.
The following source abbreviations are used throughout:

- **C §n:** `research/pde_ledger_v3/directives/O2_input_contract.md` at `217a92e9`, section n.
- **D:** `research/pde_ledger_v3/directives/O2_premise_decision_list.md` at `77d2c39a`.
- **I §n:** `research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md` at `59855a38`, section n.
- **v9:** `research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md` at `05b1a5d5`, with the source
  qualifications carried by C.

The equations, operands and domain qualifications needed here are written below. Source references
identify provenance; neither engine needs to read those files to obtain an input. Statuses belong to
the individual inputs. An equation described as **supplied** is an input the build cannot test, on
its stated domain. **Recorded**, **postulated**, **adopted premise** and **OPEN** retain their different
meanings; they are not relabelled supplied. Historical verification retains its recorded scope.

## 1. Object, setting and live quantities

**Supplied setting (C §1; I §1):** one isolated, spherically symmetric mass at rest, with its drain
flowing; far field, steady Eulerian material profiles, linear waves and lab time in the brane's
far-field rest frame. The associated optical object is leading eikonal. O2 concerns the material's
steady momentum and support accounting in that setting, rather than an optical observable or a wave
perturbation equation.

Keep `V`, `ρ_br`, `μ_⊥`, `ξ_w`, `h`, `δ`, `j_n`, the bulk state, stress and inertial responses,
reference/strain state, loading and boundary data live, including their spatial derivatives and
material-history dependence. Eulerian steadiness does not impose material constancy. No particular
radial function, constitutive family, relaxation kernel, stress measure, stress symmetry, derivative
order or finite list of internal variables is selected.

**Spherically symmetric profiles (a restriction used, from v9's setting).** The profiles are v9's
radial profiles about the mass. `ρ_br`, `μ_⊥`, `δ`, `ξ_w`, `h`, `j_n` and the bulk-density profile
`f` are general live functions of `r = |x|`, and `V` is the radial in-plane field `V^i = V_r(r) x^i/r`,
with `V_r` general and live. Component calculus may use the Cartesian `x^i` basis; the printed object
is on these profiles. The restriction is on the profiles only. It fixes no form of the OPEN operands
of §§3.2–3.3: no isotropy, parity, stress symmetry, absence of couple, chiral or rotational-reference
content (`𝒜_rot^live`), or constitutive family follows from it. Their component actions are taken on
these profiles. The transfer limits are in §7.

The component object has three in-plane directions and normal material content on the supplied
embedded brane. State the component basis and its relation to the far-field coordinates `x^i` and
bulk direction `w`. Distinguish coordinate projections from projections onto the embedded graph's
normal. Geometric component changes may be constructed from the supplied graph; a component basis
does not choose a sharp-sheet or finite-slab material reduction. Where a material or face quantity
cannot be projected without O6 or a live material response, its component action stays explicitly
unevaluated with that operand named.

**Declared measure (a convention).** Every density in this object is a density per coordinate volume
`d³x` of the far-field coordinates `x^i`, the measure on which the supplied mass law's divergence is
written (§3.1). This covers `ρ_br`, `j_n`, material momentum storage and transport, the carried and
other exchange occurrences, and the force, load, energy and power densities. In the optical ratio
`c_γ² ≡ μ_⊥/ρ_br`, `μ_⊥` is referred to the same measure as `ρ_br`. Induced-metric factors, such as
`g_ij`, `g^{ij}`, `det g_ij` and the graph normal, enter explicitly in the geometric actions; no
occurrence is moved to another measure. A quantity defined per native face area enters through its
native geometric factor and `𝒥_map` (§5). This fixes the volume measure of densities only. It
selects no stress measure, supplies no induced-measure mass balance, and does not remove the supplied
law's recorded qualification (§7).

`V` is the background in-plane velocity of the material displaced by light's `u`. It is distinct
from the perturbation velocity, outward face velocity `V_s`, bulk-normal drain `v_dr`, and a native
bulk velocity. No value for a normal material response or identification among these velocities
follows from those names (C §§2, 5–7). On the steady supplied graph `w = ξ_w` (§3.1), a brane
material point with in-plane coordinate velocity `V^i` has the bulk-direction coordinate velocity
that the graph and `V` determine. Construct that component from them; it stays live through `V` and
`ξ_w`, and there is no separate bulk-direction velocity operand. Normal velocity content that the
graph does not determine belongs to `𝒩_br^live` and `𝒥_map`. One example is different local
velocities at the faces of a finite-thickness realization, whose centre and thickness content `ξ_w`
does not fix. Such content stays an unevaluated action with those operands named. This material
velocity enters the material momentum entries (§4) and the force/power pairings (§6); the
bulk-direction carried exchange momentum is the OPEN reduction of §5.

The target is conditional on premises 1–4 below. It is a named balance with unresolved responses,
not a derived substrate law or a solution for the profiles. This specification supplies no assembled
momentum or energy balance, expected residual, sign, cancellation, profile or compatibility outcome.

## 2. Adopted premises and their accounting representation

Each row has **Status: adopted premise (user, 2026-10-06)** and the exact S21 label
**adopted substrate input to a conditional model (2026-10-06)** (C §1; D, user-selected premises).
None is a derived result.

| Premise | Content retained | Representation in this object |
| --- | --- | --- |
| **1 — material reference** | Elastic response in the optical shear regime and relaxation under steady load. The relaxation/reference evolution is a general unknown. Its power is carried explicitly, with sign undetermined and any net power naming its supplier. | The full live stress and material response keep the evolving reference/strain state and history. The energy object retains `𝒫_ref/relax^live`, `𝒮_E,net` and `𝒫_E,supply`. No relaxation law or steady-load modulus is chosen. |
| **2 — drive** | The committed dynamical order-conversion drain is the drive. No separate external body force. Its source and boundary/return data act through stress, traction and O3; the `GM` bridge is OPEN. | Use the permitted representation in which `F_drive` is **absent as a separate body-force entry**. Declare the drain-drive provenance of the material, boundary and exchange entries. Do not also add a drain-force aggregate beside them. |
| **3 — exchanged material** | Converted material carries the local brane material velocity `V`. This closes transported momentum only. Additional non-variational partners and their reaction system remain OPEN with S12. | Distinguish carried-material momentum, mechanical face/support loading and additional S12 momentum partners. Apply the outward exchange convention in §5: premise 3 closes the in-plane carried components at `V`, and the bulk-direction carried component stays the OPEN reduction stated there. Keep all live velocity and measure factors. |
| **4 — bulk traction** | Retain the postulated shear-free scalar bulk, with face-normal mechanical loading and no independent tangential bulk stress. O3 transport is separate. | The bulk part of the mechanical loading is `𝒯_bulk,n,s^live n̂_s`, with its amplitude OPEN. A tilted normal retains its in-plane projection. No tangential bulk stress is hidden in a support operand. |

The drive declaration, mechanical loading and O3 transfer are separately identifiable accounting
entries. Choosing the absent-body-force representation of premise 2 does not merge mechanical
traction with momentum carried by exchange. It also supplies neither an external support nor a
physical core holder (C §§1, 7, 10).

## 3. Live input register

### 3.1 Material identity, geometry, optical identifications and mass source

| Status and source in C | Equation or operand | Domain and use |
| --- | --- | --- |
| **Recorded material identity; supplied S9b identification**, C §2 | `u` displaces the material carrying `ρ_br`; `V` is its live background in-plane velocity. | Fixes the material whose momentum is accounted for. The homogeneous kinetic/continuity anchor is reported in §8.1; it supplies no flowing kinetic law. |
| **Supplied geometric identification**, C §5 | `g_ij = δ_ij + ∂_iξ_w ∂_jξ_w`, `g^{ij} = (g_ij)⁻¹`, `ξ_w = ℓh`. | Induced spatial metric on the recorded graph and retained L3 field identity. `ξ_w` is the brane's displacement into the bulk direction `w`, so the supplied graph is `w = ξ_w` over the far-field coordinates `x^i`. `ξ_w` has length, `h` is dimensionless, and `ℓ` is the fixed reduction scale, not a selected slab width. Geometry and slopes remain live. |
| **Supplied optical-regime live identifications**, C §6 | `c_γ(r)² ≡ μ_⊥(r)/ρ_br(r)`, `c_γ(r) ≡ c₀[1+δ(r)]`. | `μ_⊥` is premise 1's optical elastic stiffness. The ratio supplies no full stress, steady-load stiffness or inertial law. `c₀` is the supplied asymptotic light speed. |
| **Supplied steady mass balance**, C §6 | `∇·(ρ_br V) = −j_n`. | Use it as written on the declared measure (§1), retaining the full density factor and derivatives in the v9 source convention. This is the live mass input, not a proof of O6. A finite-slab or induced-measure replacement is not supplied by this equation. Its recorded relative-`O(ε)` qualification is a limit on claims, not a term in the object (§7). |
| **Supplied speed-profile anchoring**, C §6 | `Q_bg^L(x,t) = Q_bg(x)`, `Q_bg^M(x,t) = Q_bg(χ(x,t))`; v9 selects LAB_HELD for `c_γ`. | These are distinct physical anchorings on the recorded supplied background. `χ(x,t)` is the inverse material map, not `χ_B`. LAB_HELD does not impose a material-reference law or supply a holder. |
| **Supplied bulk-density inputs**, C §6 | `P = Kρ^n`, `c_s² = nKρ^(n−1)/m`, `f(r) ≡ ρ(r)/ρ₀ − 1`. | `ρ` is bulk number density, `m` particle mass and `n` a symbolic EOS exponent. The profile is unsolved; the recorded Part C response is first order in `f`. Neither `ρ_br(f)` nor a live normal traction follows. |

`c_s`, `c_γ` and the historical embedding-sector speed `c_E` remain distinct. The three Part C speed
responses are later compatibility inputs; none is selected here. The optical identifications and
LAB_HELD anchoring tag the O2 model point, rather than supplying mechanical support or a DC law
(C §§5–6).

### 3.2 OPEN material and constitutive inputs

Every entry below is a **general OPEN operand**. Notation supplies no closed argument list, locality,
instantaneous response, tensor realization, constitutive family or derivative cutoff. Live fields,
gradients and material history remain admissible dependences (C §1).

| Operand; status/source | Physical content and entry into the object | Owner retained |
| --- | --- | --- |
| `𝔅_A13`; **OPEN**, C §2 | Real/dissipative versus complex/inertial order-field branch, its action, degrees of freedom and conversion content. Branch dependence stays attached to the material and source responses; premise 1 does not decide A13. | S1's A13 gate, propagated through S5/S12. |
| `𝒯_br^live`; **OPEN form**, with character fixed by **adopted premise 1**, C §3 | Full live in-plane and normal stress, including any relaxation content. It enters the internal material-force accounting once, with the live reference/strain state. No conservative/dissipative split is adopted. | O2 carries it; S1.5/S8 antecedents, S22 nonlinear completion; relaxation ownership unassigned. |
| `𝒫_br^cons`; **OPEN conservative antecedent**, C §3 | Conservative material momentum input to the live momentum description. It is not an additional external force or a second momentum species. No live reduction into `ℐ_br^live` is supplied. | S1.5/S8. |
| `𝒯_br^cons`; **OPEN conservative antecedent**, C §3 | Conservative in-plane/normal-stress content, retained as an antecedent to the material-force description. It is not added beside the full live stress as another stress. Its relation to that stress remains unspecified. | S1.5/S8; nonlinear completion S22. |
| `ℐ_br^live`; **OPEN**, C §3 | Flowing/embedded inertial and kinetic response, entering momentum storage and transport. No identification with `ρ_br V`, a quadratic live kinetic energy, or a constant inertia is supplied. | S8, with S1.5 antecedents when available. |
| `𝒩_br^live`; **OPEN**, C §3 | Normal material response, keeping embedding, centre and thickness content distinct. It qualifies the normal material accounting jointly with the stress and inertia inputs; their overlap/identification remains unresolved, rather than being treated as three additive normal forces. | S5–S8 ingredients, Q1/Q2 static sector, S8/S22 live completion. |
| `𝒜_rot^live`; **OPEN**, C §3 | Whether internal angular momentum/couple stress is present and its physical rotational reference frame. Retain its effect on the admissible momentum/stress action and on the corresponding power accounting. No stress symmetry, vanishing couple content or chosen carrier is assumed. | S8 requirements; register assignments retain pass-2-review-pending status. |
| `ℛ_ref/strain^live`; **OPEN form**, **adopted premise 1**, C §4 | Reference/strain evolution and relaxation, including carrier, transport, formation/renewal through conversion/return, and work content. It enters through the material's evolving state/history and its explicit energy partner, not as an independently postulated body force. | Relaxation/reference ownership **unassigned**; O2 carries it. S8/S22 links are inferences; conversion/return functions remain S12's. |
| `ℳ_⊥` (O1); **OPEN**, C §6 | General stiffness response, including frequency regime and loading/material history, bulk/flow/embedding/thickness/projection dependence and gradients. It is a constitutive input to the linked optical/material description; no relation to `𝒯_br^live` is supplied. It is not a separate force. | S8, substrate reduction and S22 completion. |
| `ℛ_br` (O7); **OPEN**, C §6 | General brane-density response with live bulk, flow, embedding, thickness/projection dependence and gradients. It is input to the density in mass/momentum/energy accounting; it is not inferred from the bulk EOS or slab factorization. | S8, substrate reduction and S22 completion. |
| `ℰ_h^live` (O4); **OPEN**, C §5 | Live embedding/longitudinal relation with flow, exchange, variable coefficients and gradients. Keep it as a coupled input with `ξ_w=ℓh`. Its identity with or independence from O2's normal relation remains **unsettled**. Do not impose a second normal equation or use it as a duplicate normal force. | Q1/Q2 static ingredients; O4 live identification, S8/S22 completion. |

**Accounting convention, not a constitutive decomposition:** the material entries describe one
brane-material momentum/force accounting object. C's names supply no equation saying which response
contains or determines another. Keep the unresolved relations visible within that object. Where the
force action, momentum map, normal response or their compatibility is unsupplied, print the formal
component action of the named OPEN inputs. Do not manufacture an explicit action by choosing a
stress measure or by independently adding every named response. An auxiliary name for an unevaluated
action is notation for these existing OPEN inputs, not a new closed physical response.

### 3.3 Exchange, mechanical loads and boundary/support data

| Input; status/source | Content and use |
| --- | --- |
| `Π_n` (O3); transported-material content fixed by **adopted premise 3**, remaining partners **OPEN**, C §§1, 7, 10 | Momentum accounting for exchange. Its local carried-material convention is specified in §5. It is separate from mechanical traction; premise 3 does not close additional S12 partners. |
| `𝒥_map` (O6); **OPEN**, C §6 | Material/order weighting, live projection/window, measures, sheet/slab/face identifications and relation to bulk-normal `v_dr` and bulk/return data. It enters every native-to-v9 exchange or load reduction requiring those identifications, as detailed in §5. |
| `t_bulk,s^live = 𝒯_bulk,n,s^live n̂_s`; directional restriction from **adopted premise 4**, C §7 | Bulk mechanical traction on its native face. `𝒯_bulk,n,s^live` is a general signed **OPEN** normal-load amplitude; geometry and projection are live. No pressure/affinity/DC response law or sign of the amplitude is given. |
| `T_hold,s`; **OPEN**, C §7 | Complete mechanical face/support loading in O2. Its bulk part obeys the preceding restriction. The bulk amplitude qualifies that part of the full load; it is not a second load added beside a `T_hold,s` that already contains it. The support partition is unresolved. |
| An explicitly declared external support; **supplied held input only if declared**, C §7 | Would supply specified mechanical loads/work, with a conditional-on-that-hold domain. No such live support is selected or supplied here. Tangential external loading would have to be identified as external, not bulk shear. |
| `ℋ_core` (O5) and live mouth/core data; **OPEN**, C §§5, 7 | Physical core-holder response and selection of mouth displacement, flux/Robin or other boundary data, with any `GM` relation unresolved. They enter as boundary/support inputs; they are not solved by the local far-field object. Q2/S22 retain ownership. |
| S12 conversion/source/controller functions and additional momentum/energy partners; **OPEN**, C §§1–2, 6, 8, 10 | Native dynamical order-conversion and return functions, branch-dependent source balance, and the systems carrying reactions and supplying energy. They enter through the material, exchange and energy responses, without a new independent `F_drive` term. |
| S12 mouth/collar/return/IR and bulk-boundary data; **OPEN**, C §§1, 6–8, 10 | Native boundary/domain inputs, kept distinct from the local conversion source and its controllers. Map to the O2 face/exchange description only through the OPEN identifications. |

Names of OPEN boundary data do not fix a finite data set or choose a boundary law. A held mouth datum
in a historical static sector, or LAB_HELD speed anchoring, does not replace `ℋ_core`.

## 4. How material, force and transport content enter

The engines construct the material momentum/support relation from these physical roles, preserving
every unresolved material action. No algebraic sum for the assembled balance is provided here.

- **Material momentum and inertia:** use the live flowing/embedded response `ℐ_br^live`, with the
  conservative momentum antecedent `𝒫_br^cons` still named. Keep storage, spatial transport and
  material-history effects where the response requires them. Eulerian steadiness does not remove
  convective transport. The mass law identifies the mass-source convention; it does not close the
  momentum-to-velocity relation. The material velocity on which this response acts is that of §1,
  with in-plane components `V^i` and the bulk-direction component that the supplied graph and `V`
  determine. `ℐ_br^live` still supplies no map from that velocity to momentum, and normal content that
  the graph does not determine stays with `𝒩_br^live` and `𝒥_map`.
- **Internal material force:** use the full `𝒯_br^live` once. Retain the conservative antecedent
  `𝒯_br^cons`, evolving `ℛ_ref/strain^live` and admissibility/frame content `𝒜_rot^live` as its
  unresolved material inputs, without asserting a split or containment relation. Do not add the
  conservative antecedent or a separately invented relaxation stress to the full stress.
- **Normal material accounting:** keep `𝒩_br^live` and its unresolved relation to the stress/inertia
  actions explicit in the normal component object. Geometry can act on live material stress and
  momentum without identifying embedding, centre and width responses. Keep O4 separately named as
  a coupled relation, with equation identity/count unresolved.
- **Mechanical loading:** use `T_hold,s` as the face/support load entry, with its native bulk-normal
  restriction and OPEN support partition. Internal brane stress is not counted again as an external
  face load. Load integration/projection retains the native measures and any OPEN map. Tilted
  normal loading may have in-plane coordinate components; premise 4 is a directional restriction
  in the native face frame.
- **Exchange:** use the momentum carried with the local converted-material mass current, with the
  premise-3 velocity identification. Keep additional S12 momentum partners separately attributable
  and their reaction system named. Neither kind of transport is an extra mechanical traction.
- **Drive:** identify the dynamical drain's source/boundary provenance in these entries. The separate
  body-force representation is absent under premise 2. Do not invent a gravitational force profile
  from `GM`, density, flux or embedding amplitude.

These entry rules resolve bookkeeping roles, not the OPEN constitutive relations. If, for example,
the normal material input and a stress action are two descriptions of the same content, the
component object must retain that unresolved identification rather than count both as separate
forces. The engines may express the shared material contribution as an unevaluated joint response;
the original names and unresolved relation must remain visible. This is the conditional named
object allowed by C §1, not an engine-chosen closure.

## 5. O6, O3 and source/boundary separation

**O6 stays a general unknown (C §6; D, retained OPEN operands).** Neither a sharp-sheet reduction nor
a dynamical finite-slab region, face map, order weight or projection/window is selected. The supplied
v9 mass law can be used in its own convention while all identifications with native order loss,
face fluxes and bulk-normal flow remain `𝒥_map`-dependent.

For any representation requiring native-to-brane reduction, `𝒥_map` must carry the material/order
identification, live weights, projected/true measures and native geometry. The mass-current and
carried-momentum identifications must refer to the same exchanged material. Traction and work must
use compatible geometric/measure identifications, while remaining different physical contributions.
This consistency requirement supplies no formula or closed argument list for the map. Do not replace
it with the finite-slab equations reported in §8.4 or with the historical S11c face window.

**O3 convention (C §1, premise 3; C §6 mass-source convention):** `j_n` is the signed outward material
loss density in the supplied v9 relation `∇·(ρ_br V)=−j_n`, on the declared measure (§1). It need not
have a prescribed sign. Name the outward carried-momentum density in that same outward convention
`Π_n^carry`. Premise 3 identifies the carried velocity with the local brane material velocity. For
the in-plane `x^i` coordinate components it gives

```text
(Π_n^carry)^i = j_n V^i .
```

This is the carried-material input identification and orientation convention for those components,
not the assembled balance or a target for its residual. An engine using an inward-source notation
states the orientation change. Construct the exchange occurrence in the momentum accounting from that
convention.

On a native face description, transported momentum uses the native relative exchanged mass current
and the local brane material velocity of the converted material. A full normal/component reduction
retains `𝒥_map`, the live geometry and any unsupplied face-to-material velocity identification. Do
not choose an independent converted-material velocity, replace it with `v_dr`, or identify it with a
native bulk velocity. The premise-3 identification applies to the material after the specified
conversion/transfer; it does not fix how a native bulk current acquires that momentum. Any additional
non-variational conversion momentum partner and its reaction system remain OPEN with S12.

**Bulk-direction carried component.** Premise 3 applies to each native transfer at its own location.
It does not supply the reduced bulk-direction component `(Π_n^carry)^w` as `j_n` times one
bulk-direction velocity. In a finite-thickness realization, the transfers at different faces need
not share one current or one local bulk-direction velocity, and premise 3 supplies no such equality.
The reduced component therefore depends on the OPEN `𝒥_map` and the normal material response
`𝒩_br^live`, with the live geometry and any unsupplied face-to-material velocity identification. Print
it as an unevaluated OPEN action with those operands named, in the same outward convention.

The historical matter-stress current `Π_ij` includes the convective term `mρv_i v_j`, while
`σ^Q_ij` is its recorded quantum-stress content (C §7, stage 002 record). Its `v` is the native steady
radial bulk/reduced-lane inflow, with no supplied identification to `V`, `V_s` or `v_dr`. Convective
material transport belongs to the native momentum-current accounting; it is not an extra normal
mechanical pressure to add to `T_hold,s`. Its transfer into the O3 description requires the OPEN
material/map and S12 partner identifications. Do not use both a mapped convective current and the
same carried-material O3 current as separate transfers. Premise 4 excludes tangential mechanical
loading from `σ^Q_ij` as an adopted restriction, rather than a consequence of rest acoustics.

Maintain two distinct input inventories throughout:

| Inventory | Where it enters | What stays OPEN |
| --- | --- | --- |
| Local order-conversion/source and return-controller functions | Drive provenance, material/reference response, carried exchange and non-variational momentum/energy partners. | Branch-dependent source form, controllers, strengths/maps not supplied by C, reactions and energy supply. |
| Mouth/collar/return/IR, face and bulk-boundary/domain data | Boundary traction/work, native fluxes, O6 map and O5 data. | Boundary laws, live profiles, physical support/holder and the source-to-boundary connection. |

The committed drive is dynamical order conversion in conserved material. The historical frozen-wall
total-mass sink with remote return is not substituted for it. S14a's projected order-loss/far-field
flux bridge remains a different OPEN interface (C §§1, 6, 10; I §4).

## 6. Energy input and force/power pairing

**`ℬ_E^steady` is OPEN**, with the accounting requirement fixed by **adopted premise 1** and its S21
label in §2 (C §8). It denotes the steady energy relation to be constructed, not a precomputed
residual. Its required inputs are:

| OPEN operand (C §8) | Required content and entry | Owner |
| --- | --- | --- |
| `ℰ_br^live`, `𝒥_E^live` | Material energy storage and transport compatible with the live stress, inertia, normal response, reference evolution and rotational content. Retain live transport in the steady setting. | S1.5 antecedents, S8 material, S22 completion. |
| `𝒫_ref/relax^live` | Explicit power associated with `ℛ_ref/strain^live`, with no formula, sign or vanishing assumed. | O2 requirement; relaxation/reference ownership unassigned. |
| `𝒫_convert/exchange^live` | Order-conversion work, energy carried with exchanged material, and additional non-variational energy partners, with their reaction/supply systems. Premise 3 fixes transported momentum, not the carried total energy. | S12. |
| `𝒫_boundary^live` | Mechanical face/support work and other energy transfer through the live bulk/boundary data; any explicitly supplied hold is identified. | O2 accounting; S12 native data; Q2/S22 core response. |
| `𝒮_E,net`, `𝒫_E,supply` | Identity of the physical supplier of any net power and its stated power budget. Both stay general OPEN inputs. Naming the drain does not specify available energy or a budget. | O2 requirement; source/holder/supplier forms retain their owners. |

Pair each mechanical contribution with the velocity or generalized rate of the material point or
degree of freedom on which it acts, on the same measure (§1) and with the same geometric map as its
momentum occurrence. For a brane material point on the supplied graph, that velocity is the material
velocity of §1, the same velocity as in its momentum occurrence. In particular:

- The stress/normal-response contribution has its corresponding material stress work and energy
  transport. Any couple/frame content retains the matching rotational/generalized-rate work as
  OPEN `𝒜_rot^live` content. No ordinary elastic stored-energy functional is substituted for the
  full viscoelastic response.
- A face/support traction is paired with the actual velocity of its load application point.
  Where the live face-to-material identification is unavailable, retain it through the OPEN
  normal/material response and map. Do not use the light perturbation velocity or silently equate
  `V`, `V_s` and `v_dr`.
- Carried exchange energy and additional conversion/source power are accounted for separately from
  mechanical face work. Transported momentum at the local brane material velocity alone does not
  authorize a formula for total energy per converted mass, or an energy-free change of material
  reference.
- Keep `𝒫_ref/relax^live` explicit even when its work is represented within the material stress and
  internal-energy accounting. Identify its occurrence there instead of adding the same work again
  as an independent external power. The same rule applies when boundary or source power is also
  described within another energy operand. These names are accounting obligations, not an adopted
  additive decomposition into independent channels.

The supplier/budget operands must accompany any conditional relation requiring net supply. Their
unresolved identity or form remains visible in the output; it cannot be replaced by a declared
passive sign, reference reset, drain label or unsupported energy source. The historical order-work
and `C_ref` inputs in §8.6 retain their domains and do not close these energy operands.

No non-passive interface response is selected. If a later closure adopts one, it must name its
reservoir and state its power budget. This conditional obligation is separate from S12's
non-variational source partners and from premise 1's relaxation power. The historical S11b-C → S11c
handoff has no recorded live successor owner for it; the register's S12 target is an inference with
pass-2 review pending (C §8; I §4). Carry the obligation without inventing an owner or a response.

## 7. Counting, model point and transfer limits

**Supplied v9 counting (C §9):**

```text
ε(r) ≡ GM/(c₀²r) ,
δ = O(ε) ,                  (∂ξ_w)² = O(ε) ,
V/c₀ = O(ε^{1/2}) ,         (V/c₀)² = O(ε) ,
δ^a (V/c₀)^b ((∂ξ_w)²)^c ,  0≤a≤1 , 0≤b≤2 , 0≤c≤1
                              (nonnegative integer indices).
```

`GM` is the independent slow-test-matter orbital parameter, not a source, profile, mouth or drain
amplitude. Only `μ_⊥/ρ_br` inherits the speed-change grade. Keep the full density factor and
derivatives in the mass input. The recorded induced-metric qualification is a relative `O(ε)`
correction to `j_n`, not a live O6 law. O2 uses the law on the declared coordinate measure (§1) and
supplies no induced-metric mass balance. The qualification is carried as a limit on claims, not as a
term in the object: a claim that reads this `j_n` or `ρ_br` as a density per induced measure, or
compares it with an induced-metric mass balance, carries it. Fixed `ℓ` transfers the live `h` slope
to the geometric grade without supplying a mouth-amplitude grade. First order in bulk `f` is a
separate recorded response domain; no relation between `f` and `ε` is supplied.

**OPEN grades and derivative scales (C §9):** individual density and stiffness; inertia; stress and
normal response; reference/relaxation and power; source and exchange momentum; force and
traction/support; holder and mouth data; embedding-sector coefficients/source and longitudinal
field, wherever no grade is recorded. Retain these as named unknown grades/scales attached to their
inputs. Do not assign them a convenient higher order or separate density/modulus `O(ε)` variations.

O2 is **untruncated**. The optical monomial box records the model's optical counting; it is not an
order contract for removing material, exchange, boundary or power terms from O2.

| Restriction or freeze | Use and transfer limit in this specification |
| --- | --- |
| Far field, spherical isolated mass at rest, Eulerian steady profiles, lab time; linear optical waves and leading eikonal | Supplied model setting. No transfer to strong field, mouth/interior, moving or rotating mass, or a drain varying during a light crossing (C §1; I §1; v9 model point). |
| Spherically symmetric profiles: scalar profiles general functions of `r`, radial `V` (§1; v9's live radial profiles, v9:43, 53–54, 85–86, 91, 113–114, 144) | Used: the printed object is on these profiles. It is v9's supplied setting, not a new premise; the radial functions stay general and live. No transfer to profiles without that symmetry, including swirl or other angular dependence of `V` and angular dependence of a scalar profile. It fixes no form of the OPEN operands (§1): isotropy, parity, stress symmetry, and couple-stress, chiral or rotational-reference content (`𝒜_rot^live`) stay OPEN. |
| Isotropic optical speed, same shear-regime identification | Recorded optical scope. Direction-dependent stiffness, polarization-dependent propagation, coupled thickness/bulk branches, extra mixed `ωk` content and polarization/subprincipal transport remain later optical questions. This scope does not delete thickness/bulk dependences from OPEN mechanical responses (C §§4, 6; I §2). |
| S9 sharp zero-width sheet, continuum and vanishing wave amplitude | Historical optical limits. Continuum/linear optical scope is retained as recorded; no sharp-sheet material or exchange reduction is selected for O2 (C §§4, 6). |
| S9 no dissipation and frequency-independent moduli | **Revisited by premise 1**. Not imposed on the steady-load/reference relaxation response. Optical consequences remain deferred (C §4). |
| S9 zero in-plane background flow and zero background strain | v9 lifts flow and the isotropic strain freeze; `V` and the speed profile stay live. No direction-dependent optical extension is supplied (C §§1, 4). |
| LAB_HELD speed profile | Supplied anchoring in space; no material-reference constancy, physical support or physical holder follows (C §6). |
| S11b rest bulk, active normal drain excluded, uniform quadratic material inputs | Historical perturbation domain only. The uncarried normal-background-flow correction is not a live bulk law; its freeze is not applied to O2's `j_n` or OPEN `v_dr` map (C §§3, 7). |
| S11b breathing `k=0`, impermeable faces, no reciprocal traction | Separate historical thickness slice only; not a normal-flow or embedding closure (C §3). |
| S11c supplied profiles, first shape/background-jet order, frozen background current/exchange; uniform representative fold | Historical scope only, with full composition unresolved. No live coefficient substitution, density freeze, support balance or response is inherited (C §§3, 6–7). |
| L3 fixed reduction scale and postulated parent; L4 constant coefficients; L5 frozen sleeve/profile class; L6 static source-free constant-stiffness exterior and held mouth datum | Reported static/postulated sector only. Field identity is retained; governing equations and datum are not promoted to live flow or a physical holder (C §5). |
| Dynamical order-conversion drive, shear-free face-normal bulk restriction, optical-elastic/steady-relaxing material and local-velocity exchange | Adopted conditional-model premises 1–4, not inferred substrate properties. No transfer to a different branch, bulk shear law, reference law or exchange velocity without new physical input (C §1). |

The retained-reference/newly relaxed-reference comparison remains **EXPLORATORY / PAUSED**. Its
Claude-only verdict is literally **COHERENT CONDITIONAL COMPARISON**, for the qualified conditional
paper calculation. Its `B_*`, fixed-patch transport, formation data and proposed comparator laws are
not adopted inputs. The 2026-10-06 premise-1 choice supersedes the earlier compare-before-relaxation
ordering but selects neither comparator law and does not resume the comparison (C §4).

## 8. Recorded equations with restricted domains

These entries complete the standalone input packet. They are reported as equations or named
operands with C's statuses. They do not supply missing live material, map, embedding, exchange,
support or energy laws, and they are not extra terms to append to the live object. No historical
solved profile or prior-art result is an input.

### 8.1 Material identity and ontology (C §2)

**Recorded material identity; supplied S9b identification:** the homogeneous linear-displacement
anchor is

```text
T_u = ½ ρ_br |∂_t u|² ,       δρ_br = −ρ_br ∇·u .
```

This kinetic/continuity identification is restricted to homogeneous linear material displacement;
the second equation is not sourced steady continuity. S8's register target for material identity
and thickness inertia is an inference with pass-2 review pending, not a new live kinetic law.

**Postulated material-state ontology:**

```text
ρ = |ψ|² ,                   n_B = χ_B n_006 .
```

The historical ordered/disordered states are postulated states of one conserved material, without
deriving a wall. `n_006` is stage 006's conserved constituent number density, not v9's EOS exponent;
its velocity symbol `u` is distinct from the displacement above. The real-fraction split and its
order balance do not decide the OPEN A13 branch. Owners: S1, S5 and S12.

### 8.2 Quadratic material and wall antecedents (C §3)

**Supplied quadratic ingredients, reported only on the recorded uniform linear domain:**

```text
e_W ≡ δW/W₀ ,
U = ½ μ_R |∇×u|² + ½ B_ρ⁽³⁾ θ² + C W₀ θ e_W
    + ½ k_W W₀² e_W² + ½ κ_W W₀² |∇(δW)|² ,
T = ½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)² ,
μ_θ ≡ δU/δθ ,                p_W ≡ δU/δe_W ,
δ_vθ + δ_ve_W + ∇_x·δ_vu = 0 .
```

These are carried terms, not a complete basis or live law. Functional derivatives hold other
fields fixed. The final equation is the uniform linearization of an instantaneous no-transfer
**virtual** material-mass constraint, not physical evolution. Uniform inertia/thickness inputs and
the breathing slice have their own restrictions. `ζ_c` and `W` are independent centre/thickness
variables, not replacements for `ξ_w`.

S11c's supplied-profile first-background-jet operators have scoped repair support and unresolved
full composition. Their uniform representative fold does not extend to varying coefficients.
Neither those operators nor substituting live coefficients into `U` and `T` closes any live or
conservative material operand in §3.2.

**Recorded static wall tension:**

```text
σ_wall = ∫ dw κ_B (χ_B′)² .
```

The two-engine static one-dimensional single-kink checks were relative to postulated wall terms.
This is not a slab tension, width selection, G0 wall Hessian or flowing-brane stress. S6 owns tension;
S7 owns width. It supplies no live normal-response value.

### 8.3 Embedding reduction and static governing inputs (C §5)

**Supplied reduction, conditional on a postulated localized parent sector:**

```text
f₀(w) = 1/[ℓ cosh²(w/ℓ)] ,            N₀ = ∫ dw 2f₀² ,
h = P₀H ≡ N₀⁻¹ ∫ dw 2f₀H ,
M_h = N₀M₄ ,                         K_h = N₀K₄ = M_h c_E² .
```

The reduction was verified within its chosen postulated action. `{M₄,c_E}` are parent inputs,
not tension-derived normalization. The static/postulated Q1 domain supplies no live inertia.

**Recorded static-sector governing inputs, not adopted live equations:**

```text
A_eff = ρ_br + C_J²/κ_phase ,
S_Lh = ∫ dt d³x [½ A_eff (∂_t u_L)² + ½ M_h (∂_t h)²
                − ½ B_eff |∇u_L|² − ½ K_h |∇h|² − C_hu ∇u_L·∇h] ,
(δΩ/δh)_mouth = η_i(k_m h − g_χh s_i) ,
k_m = K_m/ℓ² ,          g_χh = J_m/ℓ ,          Q_χ[r_Σ,s_i] = s_i ,
d/dr(r² dh/dr) = 0 ,    h(a) = h_A ,            h → 0 at infinity ,
h_A ≡ ξ_w|_A/ℓ = P₀H|_A .
```

The action has constant coefficients and a postulated G0 parent. The mouth projection has a frozen
sleeve/profile class. `s_i` is puncture orientation, not a Part C response exponent; `J_m` is mouth
coupling, not `j_n`. The exterior is source-free, static, with constant generic `κ_ext>0`, distinct
from `μ_⊥` and the record's response ratio. No solved exterior profile/amplitude is supplied. There
is no map from `u_L` to steady `V` or from this mouth source to the mass's drain. Q1/Q2 retain their
domains; live extension is O4 and holder selection O5.

### 8.4 Density, slab and projection records (C §6)

**Recorded uniform transverse anchor:**

```text
ρ_br⁰ω² = μ_⊥k² .
```

Its basis qualifications do not define a varying-coefficient stress. The live optical ratio in
§3.1 is supplied separately.

**Supplied slab/projection kinematics, on the recorded finite-slab domain:**

```text
Σ_E = ρ_4D W ,      Σ_mat(X,t) = Σ_E(x(X,t),t) 𝒥_x(X,t) ,
δ_vΣ_mat = 0 ,
J_s^α = ρ_m(v_bulk,s − v_face,s^α)·n̂_s^α ,
a_s^α = sqrt(1+|∇_x h_s^α|²) ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋)                         (flat faces) ,
∂_tΣ^α + ∇_x·(Σ^α v) = −∑_s a_s^α J_s^α ,          v = ∂_t u .
```

`ρ_4D` and `ρ_m` are mass densities; `α` labels the recorded anchoring; `h_s^α` is a native face
graph, distinct from reduced `h`. `J_s^α` is outward relative mass flux per true face area. The
shape verification is confined to supplied profiles and first shape order, with background
current/exchange frozen. It does not verify a live O6 map. The virtual constraint is distinct from
physical sourced evolution; `v=∂_t u` does not identify steady `V` with light's perturbation velocity.

The historical `ρ_br=ρ_4D W` is slab kinematics. **RHO4-CONSTANT** fixes `ρ_4D,bg⁰`, whereas
**RHOBR-CONSTANT** fixes `ρ_br,bg⁰` while thickness varies. Neither representative is adopted, and
neither supplies `ℛ_br`. The historical S11c window is tied to its own face maps, not selected here.

**Postulated/PENDING shear projection**, with dimensional consistency asserted only:

```text
μ_R = ∫ dw χ_B μ_R⁽⁴⁾ .
```

This is neither an `ℳ_⊥` law nor a complete live map. S12 retains conversion and the separate native
source/boundary inventories; S14a retains the projected order-loss/far-field flux bridge.

### 8.5 Bulk, perturbation traction and held-support records (C §7)

**Postulated shear-free scalar bulk retained by adopted premise 4. Supplied rest-frame acoustic
model, on its rest-linearized domain:**

```text
v_bulk = ∇₄φ ,       δp = −ρ_m ∂_tφ ,       ∂_t²φ = c_s0² ∇₄²φ .
```

There is no bulk shear modulus; the record uses outgoing or decaying radiation conditions. Active
`v_dr` is excluded and the background-normal-flow correction is uncarried. This supplies no live
background pressure/load response.

**Supplied perturbation traction and face work, on their original linear/frozen-background domain:**

```text
𝒜_s = μ_θ/ρ_br⁰ − δp_s/ρ_m ,
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,       Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,
t_s = −(δp_s+Λ_X(ω)𝒜_s)n̂_s ,
δ_v𝒲_bulk^α = ∑_s a_s^α t_s^α·δ_vx_s^α .
```

The prescribed linear response constants/times are independent and real; the nonuniform extension
has first shape order with frozen background. These kernels are not a material-reference law, and
setting `ω=0` does not supply live normal loading or a background `j_n` law. The virtual work records
the native traction/geometry pairing, not a live energy closure.

**Recorded supplied support bundle on a held background:**

```text
𝒮_hold⁰ = {f_hold⁰,t_hold,s⁰} ,       V_s⁰ = J_s⁰ = 𝒜_s⁰ = 0 .
```

Its support-stabilized supplied-profile energy/geometry comparison supplies no O2 law with live
flow. The freeze is reported, not imposed here; neither it nor LAB_HELD selects a physical holder.

### 8.6 Order work and energy-reference records (C §8)

**Postulated historical real-fraction order-work contribution:**

```text
μ_χ = δF/δχ_B ,       P_order = ∫ d⁴X μ_χ D_tχ_B .
```

There is no extra number-density factor in this recorded contribution. It is restricted to the
postulated real-fraction free-energy/dynamics adjunct; it does not supply branch-dependent S12
completion, relaxation power or full steady energy accounting. The paused comparison does not
establish re-ordering as an available power source.

**Recorded energy-reference question, OPEN choice `C_ref`:** the historical action postulates
`U(ρ)=Kρ⁵/4`, while S1.5 records, specifically for its `n=5` EOS,

```text
P = ρU′ − U = Kρ⁵ ,       U(ρ) = Kρ⁵/4 + C_ref ρ .
```

`C_ref` is the unresolved chemical-potential/energy-reference choice, distinct from the quadratic
density–thickness coupling `C`. Do not select `C_ref=0` or extend this expression to v9's symbolic
`n`. S1.5 also owns the linked momentum-stress/quantum-energy/current improvement convention. The
reference choice remains named wherever conservative energy and source/exchange accounting depend
on it; choosing an EOS does not close that dependence.

## 9. Per-engine printed objects

Both engines independently construct the object from this packet. Their scripts print computed
objects and domain/input qualifications; interpretation belongs to the later record. There is no
expected-value acceptance criterion.

Print the following, with in-plane components and normal content separately identifiable:

| Printed object | Input trace and qualifications required |
| --- | --- |
| Component basis, native/reduced measures and geometric projection objects used | The declared coordinate measure for every density (§1); `g_ij`, `ξ_w=ℓh`, live geometry; distinguish a geometric change of basis from an unresolved material/face reduction through `𝒥_map`. |
| Material momentum storage/transport and inertial component actions | `ℐ_br^live`, `𝒫_br^cons`, live `V`/`ρ_br` where applicable, and the material velocity on the supplied graph (§1); state the OPEN momentum map and any unsupplied stress/normal-response identification, including normal content that the graph does not determine. Do not emit a chosen kinetic law as computed. |
| Internal material-force component actions, including normal-response and rotational content | `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live`; show unresolved relations without duplicate forces or an adopted constitutive split. |
| Mechanical face/support load components | `T_hold,s`, the OPEN `𝒯_bulk,n,s^live` bulk-normal restriction, geometry and map; preserve support/holder status and native measures. |
| Carried exchange-momentum components and any separately attributable source-partner actions | `Π_n`, premise 3's in-plane `j_n V^i`, the bulk-direction component as the OPEN reduction of §5 with `𝒥_map` and `𝒩_br^live` named, O6/native relative fluxes where needed; keep S12 partner and reaction systems OPEN. Mechanical load and native convective transport are not counted again. |
| Drive representation and source/boundary dependencies | Premise 2's absent separate body-force representation; drain action traced into stress, loading and O3, with local source/controller and boundary data separately named. No computed `F_drive(GM)` is implied. |
| The constructed `ℬ_hold^live` component objects | Each occurrence traced to the preceding physical entries, with its orientation, status, domain and remaining OPEN actions. No unresolved overlap is silently treated as an additive constitutive decomposition. |
| Mass input and any transformations actually used with the component objects | Supplied live `∇·(ρ_br V)=−j_n` on the declared coordinate measure (§1), full live density/derivatives and its recorded relative-`O(ε)` qualification (§7); any native-measure or native-flux identification remains O6-dependent. Print operands and transformed objects. |
| Force/power pairings and the constructed `ℬ_E^steady` object | All §6 energy operands, matched velocities/rates/measures (§6), explicit relaxation power, exchange energy, boundary work, and OPEN net supplier/budget. Preserve `C_ref` dependence where applicable and any conditional non-passive reservoir obligation. |
| Coupled-input and model-point qualifications | A13, O1, O4, O5, O6, O7, missing grades/scales, premises 1–4, the spherically symmetric profiles (§1), the restrictions in §7 and the historical-only domains in §8. O4 equation identity/count stays unsettled. |

Each term carries the section of this specification and the C section from which its input came.
OPEN actions are printed as OPEN actions with their physical role, rather than a substituted closed
formula or an absent term. An unperformed native reduction or unavailable constitutive action stays
named in the computed conditional object. If a historical relation is displayed alongside it,
print that relation's original domain separately; do not attach a historical static or frozen
equation to the live component object as an additional governing equation.

## 10. Interfaces and finite boundary

**Filing (D):** O2 is a bounded sub-step of its own. It declares these interfaces and owns none of
the downstream steps.

| Interface | Handoff and limits (C §10; I §4) |
| --- | --- |
| **S12** | Conversion/source/controller inventory, distinct boundary/domain inventory, O6 identifications, and additional non-variational momentum/energy partners with reaction/supply systems. O2's local-velocity transported momentum does not determine those partners, a return law, branch choice or energy supply. |
| **S14a / S14** | Projected order-loss source, drain label `J`, profile-dependent far-field flux map and controlled DC return bridge remain S14a obligations. S14 is **CONDITIONAL ON S14a**, with its recorded frozen DC-sink/zero-mode completion domain. Those completions and far-field profiles are not O2 premises. O2 supplies a conditional live material/load object, not that bridge or a source-to-`GM` identification. |
| **S16** | Response-side connection to the independent orbital `GM`. The inventory's designation of S16 as this matching interface is an **inference**. Its compact-body, supplied external-field/worldtube approximation retains calibrated potential, supplied mass/multipoles and discarded boundary-flux assumptions; its normalization is target-matched, not derived. A response-side match does not derive the drain's source or coupling. |
| **S21** | Requirements/provenance and every OPEN operand, with premises 1–4 labelled **adopted substrate input to a conditional model (2026-10-06)**. The later record carries register entries. S21 owns the revising-input/new-consequence sort; this spec performs no integration or substrate-completeness verdict. |
| **S1.5 / S8 / Q1–Q2 / S22** | Conservative/material antecedents, energy-reference/improvement convention, rotational obligations, static embedding qualifications, nonlinear material completion and physical core holder retain their stated owners. Unassigned relaxation/reference ownership remains unassigned. |

Prior art remains an oracle for later independently constructed objects in matched domains, never
a premise, constitutive choice or expected result (C §10; I §5).

**Deferred to the build:** executable representation of general response actions; component and
measure calculus retaining every live dependence; independent construction and literal emission
of the component and paired-energy objects; and executable script-control tests, including
variable-coefficient/FORM and independence checks. Methods, staging, resources, serialization,
transcript handling and guards belong to the later build directive. This spec selects no CAS route
or implementation. That directive cannot turn an OPEN physical input into a chosen response.

Constitutive elimination, additional premise selection, S11c repair/composition, S14a bridge work,
source/response `GM` matching, physical holder selection/solve, S21 integration and optical
compatibility remain at their named later steps. The named untruncated conditional balance is the
finite construction target. If a more explicit action requires a physical restriction not present
here, retain the named OPEN action where possible; otherwise report the missing restriction as a
question for the user, without choosing it. A second method failure, new unnamed sub-problem or
repair-to-repair triggers the repository stop rule.

**Authoring STOP:** this file is the sole deliverable. No assembled balance, computation, build,
review round, register edit or commit is part of this authoring task.
