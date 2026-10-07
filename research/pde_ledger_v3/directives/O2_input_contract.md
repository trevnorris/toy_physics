# O2 — input contract

**Author:** Codex, 2026-10-06. **Status:** authored; no review clearance claimed.

**Deliverable:** one input contract for O2's brane material, conservative momentum and normal-stress operands,
inertia and normal material response, reference evolution, geometry, density/stiffness/projection inputs,
bulk-traction premises and steady-state energy accounting.

This is sub-step 2 of `O2_steady_brane_balance_scoping.md:248` (`59855a38`), under
`O2_input_contract_directive.md` (`7bee773b`), applying `O2_premise_decision_list.md` (`77d2c39a`). Paths below
are relative to `research/pde_ledger_v3/`; `v2:` means `research/pde_ledger_v2/notes/`, and `docs/` means the
repository-root directory. **v9** means `directives/S9b_SHARED_PHYSICS.md` at `05b1a5d5`, with that version's
line numbers. Other source lines are from the files at `7bee773b2ea7dd0cdb2a601aed5449eb0bddb9cc`.

## 1. Domain and adopted premises

**Status: supplied setting.** One isolated, spherically symmetric mass at rest, with its drain flowing;
far field, steady material profiles, linear waves, and lab time in the brane's far-field rest frame
(v9:113–114). The optical object is leading eikonal (v9:66–69). Material velocity `V`, `ρ_br`, `μ_⊥`,
embedding `ξ_w`, speed change `δ`, normal exchange `j_n`, bulk state and their spatial derivatives stay
live. Steady Eulerian profiles do not supply material-reference constancy. The setting supplies no
constitutive response for the flowing state. **Owner:** bounded O2, with the interfaces recorded below.

Each of premises **1–4** has **Status: adopted premise (user, 2026-10-06)** and the decision list's exact
**Label: adopted substrate input to a conditional model (2026-10-06)**, for S21's later sort
(`directives/O2_premise_decision_list.md:12–14`). None is a derived result.

| Premise | Adopted content and source | Contract location / owner |
| --- | --- | --- |
| **1 — material reference** | Elastic response in the optical shear regime and relaxation under steady load; the relaxation response remains a general unknown. Steady energy accounting includes its power, with sign undetermined and any net power naming its supplier (`O2_premise_decision_list.md:16–27`). | Reference and energy inputs in §§4, 8. O2 carries the requirement; material content goes to S8, with nonlinear completion at S22 and source partners at S12. No specific relaxation law is supplied. |
| **2 — drive** | The committed order-conversion drain is the drive, with no separate external body force. Its action is through stress, traction and O3 operands; `F_drive` is identified with those terms or absent beside them. The `GM` bridge stays OPEN (`:28–31`). | Drive specification belongs to **sub-step 3**. S12 owns conversion/source and separate boundary data; S14a/S14 and S16 retain their source/response interfaces. This contract makes no choice between the two permitted representations of `F_drive`. |
| **3 — exchange momentum** | Converted material carries the local brane material velocity `V`. This closes transported momentum only; additional non-variational momentum partners and their reaction system stay OPEN with S12 (`:32–36`). | O3's momentum term and its convention belong to **sub-step 3**; no equation for `Π_n` is authored here. S12 owns the additional partners. |
| **4 — bulk traction** | Retain the postulated shear-free scalar bulk, with face-normal loading and no tangential bulk stress; O3 transport remains separate (`:37–39`). | Bulk premises and the OPEN normal-load response in §7. This does not fix the full `T_hold,s`, its physical support, or its power. |

Every **OPEN** object named below is a general operand. Naming it imposes no closed argument list,
locality, instantaneous response, finite set of internal variables, derivative order or constitutive
family beyond a cited record. Dependence on live fields, their gradients and material history is not
removed by notation. The steady setting does not license freezing those dependences.

The intended later claim is the **untruncated conditional O2 balance**, with `V` and `j_n` live and every
OPEN operand named (`O2_premise_decision_list.md:41–56`). The drive, O3 term and full in-plane/normal
face/support balance are sub-step 3's object, rather than additional equations supplied by this contract.

## 2. Material identity and branch

**Status: recorded material identity; supplied S9b identification.** `u` displaces the material carrying
`ρ_br`; `V` is that material's live background in-plane velocity. The recorded homogeneous kinetic and
continuity identification is

```text
T_u = ½ ρ_br |∂_t u|² ,       δρ_br = −ρ_br ∇·u .
```

**Source/domain:** `steps/S11_stray_longitudinal.md:32–38`; homogeneous linear material displacement.
The second relation is not the sourced steady continuity law. v9:15–16,30–45 supplies the material
identification on its live-profile object. `V` is distinct from light's perturbation velocity, outward
face velocity `V_s`, and bulk-normal `v_dr` (`directives/S11b_SHARED_PHYSICS.md:89–111`).
**Owner:** S8 for the brane material/kinetic content (`V3_STEP_PLAN.md:332–350`). The register's S8 target
for the material identity and thickness inertia is explicitly an inference, with pass-2 review pending
(`SUBSTRATE_REQUIREMENTS.md:3,386–400`); it does not supply a new live kinetic law.

**Status: postulated material-state ontology.** The historical bulk condensate and real-fraction split
record

```text
ρ = |ψ|² ,                    n_B = χ_B n_006 .
```

Here `n_006` is stage 006's conserved constituent **number** density `n`, not v9's EOS exponent; its
velocity `u` is also distinct from the displacement `u` above. **Source/domain:**
v2:`stages/ledger_stage004_gnls_action_dimensional_foundation.md:58–69` and
v2:`stages/ledger_stage006_two_phase_chiB_ontology.md:16–26,81–98`. The latter postulates ordered brane
and disordered bulk as states of one conserved material; it does not derive the wall from that material.
Its real-fraction range and order balance belong to its stated split, rather than selecting the live
branch. **Owner:** S1; branch-dependent wall content S5 and conversion content S12.

**Status: OPEN branch operand `𝔅_A13`.** Real/dissipative versus complex/inertial order-field content,
and the corresponding action, degree count and conversion balance, remain unresolved. Premise 1's
viscoelastic character does not decide A13. **Form/domain:** no branch is selected; the S5 expression is
explicitly a template. **Owner:** S1's A13 gate, propagated at S5 and S12
(`V3_STEP_PLAN.md:173–180,271–277,616–617`; `O2_premise_decision_list.md:43–44`).

## 3. Conservative momentum, stress, inertia and normal response

**Status: OPEN live material operands.** Retain separately named content:

| Operand | Physical input whose form is missing | Owner / boundary |
| --- | --- | --- |
| `𝒫_br^cons` | Conservative brane material momentum content. | S1.5 supplies conservative substrate balance antecedents; S8 supplies brane material content. No live reduction is supplied. |
| `𝒯_br^cons` | Conservative stress content, including the in-plane and normal-stress operands required by O2. | S8; nonlinear material completion S22. No stress measure, symmetry or full constitutive form is selected here. |
| `ℐ_br^live` | Inertial/kinetic response on the flowing, embedded material background. | S8, using S1.5 antecedents when available. The optical density identification alone does not supply this response. |
| `𝒩_br^live` | Normal material response, retaining distinctions among embedding, centre and thickness content until a model connects them. | S5–S8 for wall/width/compression ingredients; Q1/Q2 for their recorded static embedding sector; live identification with O4 remains OPEN. |
| `𝒜_rot^live` | Required internal angular-momentum/couple-stress content and the physical rotational reference frame, including whether such content is present. | S8's OPEN requirements; no carrier or frame is chosen. |

**Source/domain:** `O2_premise_decision_list.md:43–47` retains the live stress, inertia and normal response
as general unknowns. `V3_STEP_PLAN.md:185–204` assigns S1.5 conservative left-hand sides and keeps S12's
source partners OPEN; `:271–350,1207–1214` supplies the other material boundaries. The rotational
obligations remain OPEN in `SUBSTRATE_REQUIREMENTS.md:322–367` (pass-2 review pending at `:3`). This
contract carries those obligations without performing a material audit or adopting an active stress.

**Status: supplied quadratic ingredients, reported only on their recorded domain.** S11b carries

```text
e_W ≡ δW/W₀ ,
U = ½ μ_R |∇×u|² + ½ B_ρ⁽³⁾ θ² + C W₀ θ e_W
    + ½ k_W W₀² e_W² + ½ κ_W W₀² |∇(δW)|² ,
T = ½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)² ,
μ_θ ≡ δU/δθ ,                p_W ≡ δU/δe_W ,
δ_vθ + δ_ve_W + ∇_x·δ_vu = 0 .
```

**Source/domain:** `directives/S11b_SHARED_PHYSICS.md:255–299,318–352`. The displayed `U` is the set of
carried terms, not a complete basis or a live material law; additional invariants and representative
qualifications are recorded at `steps/S11b_interface_coupling_law.md:95–121`. The functional derivatives
hold the other fields fixed. The last equation is the uniform linearization of an instantaneous,
no-transfer **virtual** material-mass constraint, not a physical evolution equation. Inertia and the
thickness response are uniform linear model inputs; the breathing record has the separate slice
restrictions `k=0`, impermeable faces and no reciprocal traction
(`steps/S11bB_interface_assembly.md:80–85`). None supplies `𝒩_br^live` for `ξ_w`.

S11c's supplied-profile, first-background-jet operator content has scoped repair support and unresolved
full composition (`steps/S11c_b_variable_coefficient_operator.md:35–57`). Its uniform representative
fold does not extend to varying coefficients (`directives/S11c_b_SHARED_PHYSICS.md:165–174`). Neither
those operators nor a substitution of live coefficients into the displayed quadratic terms closes
`𝒫_br^cons`, `𝒯_br^cons`, `ℐ_br^live` or `𝒩_br^live`.

**Status: recorded wall tension, with static qualification.**

```text
σ_wall = ∫ dw κ_B (χ_B′)² .
```

**Source/domain:** v2:`stages/ledger_stage006_two_phase_chiB_ontology.md:95–97,133–154`; both engines
verified the static one-dimensional single-kink EL residual and tension integral relative to the
postulated wall terms. This is not a slab tension, width selection or flowing-brane stress. It is also
not the G0 wall Hessian (v2:`stages/ledger_stage030_electric_scalar_localized_h_closure.md:248`).
**Owner:** S6; width S7 (`V3_STEP_PLAN.md:285–288,315–330`). It supplies no value for a live normal operand.

## 4. Material-reference and strain evolution

**Status: adopted premise 1 (user, 2026-10-06); Label: adopted substrate input to a conditional model
(2026-10-06).** The optical shear response is elastic while the material relaxes under steady load.
**Status of its form: OPEN `ℛ_ref/strain^live`.** This single name denotes the unsupplied reference/strain
evolution and relaxation response, including its physical reference carrier, transport, formation or
renewal through conversion/return, and work content where applicable. It has no prescribed argument
list, tensor realization, relaxation kernel, rate, time scale or zero-frequency form. Premise 1 fixes
character, not an equation (`O2_premise_decision_list.md:16–25,44`).

**Recorded limits:** S9 took no dissipation and frequency-independent moduli, as well as a sharp sheet,
rest background, continuum and vanishing wave amplitude (`steps/S9_light_requires_shear.md:349–350`).
v9:93–106 lifts in-plane background flow and isotropic background-strain freezes while retaining its
other stated scope boundaries. Premise 1 expressly revisits the no-dissipation/frequency-independent
limits; those limits cannot be imposed on its steady-load relaxation response. The consequences for
optical propagation remain a later light-compatibility question, as required by the decision list
`:20–27`. LAB_HELD speed anchoring is not a reference-evolution law.

**Status: recorded exploratory comparison, not an adopted response.** The retained-reference/newly
relaxed-reference comparison is **EXPLORATORY / PAUSED** (`steps/S11c_PARTIAL_CLOSEOUT.md:35–37`;
`docs/light_em_investigation_handoff.md:7,11,19`). Its Claude-only verdict is literally
**COHERENT CONDITIONAL COMPARISON**, for the conditional paper calculation, with the qualifications in
`docs/elastic_reference_comparison_assessment.md:3–5,41–55`. Its proposed `B_*`, fixed-patch transport
and formation data are not inherited inputs (`docs/elastic_reference_limit_comparison.md:20–40`).
Neither comparator case supplies a relaxation rate or a complete physical stress. The earlier user
steer authorized comparing those limits before adding relaxation; the 2026-10-06 premise-1 choice
supersedes that order (`O2_premise_decision_list.md:16–19`). It selects neither comparator law and does
not resume the paused comparison.

**Owner:** O2 carries this input under premise 1; material/reference content belongs with S8, and
nonlinear completion with S22. The paused comparison remains an input for Q2/S22 with S12 connections;
S12 owns the actual conversion/return functions and partners, rather than a relaxation law supplied here.

## 5. Geometry, embedding and O4

**Status: supplied geometric identification.**

```text
g_ij = δ_ij + ∂_iξ_w ∂_jξ_w ,       g^{ij} = (g_ij)⁻¹ ,       ξ_w = ℓh .
```

**Source/domain:** v9:19–21,46–48,216–228; the induced spatial metric on its supplied graph and the
retained L3 field identity. Underlying identity:
v2:`stages/ledger_stage031_puncture_deflection_field_identity_source.md:60–76`, within the postulated
G0 sector. `ξ_w` is a length, `h` dimensionless and `ℓ` that reduction's fixed scale; `ℓ` is not a newly
selected slab width. `ξ_w`, `h` and their spatial derivatives remain live. `ζ_c` and `W` are independent
face-centre/thickness variables in S11b, rather than replacements for `ξ_w`
(`directives/S11b_SHARED_PHYSICS.md:89–96`). **Owner:** static identity Q1/Q2; live applicability O4.

**Status: supplied reduction, conditional on a postulated parent sector.**

```text
f₀(w) = 1/[ℓ cosh²(w/ℓ)] ,            N₀ = ∫ dw 2f₀² ,
h = P₀H ≡ N₀⁻¹ ∫ dw 2f₀H ,
M_h = N₀M₄ ,                         K_h = N₀K₄ = M_h c_E² .
```

**Source/domain:** v2:`stages/ledger_stage030_electric_scalar_localized_h_closure.md:46–51,63–66,88–112,127–130`.
The projection and reduction are verified within that chosen localized parent action; the parent is
postulated (`:15–22,228–233`). `{M₄,c_E}` are its inputs, not a tension-derived normalization. `c_E`,
`c_s` and `c_γ` have no supplied identification (`:259–261`). **Owner:** Q1, static/postulated domain
(`V3_STEP_PLAN.md:871–893`); this reduction is not a live inertial law.

**Status: recorded static-sector governing inputs, not adopted live equations.** The records give

```text
A_eff = ρ_br + C_J²/κ_phase ,
S_Lh = ∫ dt d³x [½ A_eff (∂_t u_L)² + ½ M_h (∂_t h)²
                − ½ B_eff |∇u_L|² − ½ K_h |∇h|² − C_hu ∇u_L·∇h] ,
(δΩ/δh)_mouth = η_i(k_m h − g_χh s_i) ,
k_m = K_m/ℓ² ,          g_χh = J_m/ℓ ,          Q_χ[r_Σ,s_i] = s_i ,
d/dr(r² dh/dr) = 0 ,    h(a) = h_A ,            h → 0 at infinity ,
h_A ≡ ξ_w|_A/ℓ = P₀H|_A .
```

**Source/domain:** action at stage 030 `:134–142`, with its constant coefficients and postulate boundary
`:228–233`; mouth inputs at stage 031 `:68–76,96–138,147–155`; exterior equation/data at `:157–174`.
The mouth projection assumes the recorded frozen sleeve/profile class (`:27–30,243–253`). The exterior
is source-free, static and has constant generic exterior stiffness `κ_ext>0`, distinct from `μ_⊥` and
from the record's response ratio. `s_i` is puncture orientation, not Part C's response exponent;
`J_m` is a mouth coupling, not `j_n`. No solved exterior profile or amplitude is supplied here.
**Owner:** Q1/Q2 with their static domains (`V3_STEP_PLAN.md:871–915`).

**Status: OPEN `ℰ_h^live` (O4).** It denotes the missing live embedding/longitudinal relation with flow,
exchange, variable coefficients and their gradients (v9:295–300). The recorded action supplies no map
from `u_L` to steady `V`, or from its mouth source to the mass's drain. **Relation to O2's normal content:
unsettled**, as fixed by `O2_premise_decision_list.md:46–47`. O4 remains separately named as a coupled
input; this contract chooses neither an identification with `𝒩_br^live` nor an independent equation
count. **Owner:** recorded static ingredients Q1/Q2; the live identification remains O4 for the later
O2 spec, with material completion at its S8/S22 owners.

## 6. Density, stiffness and projection

**Status: supplied live identifications and balance.**

```text
c_γ(r)² ≡ μ_⊥(r)/ρ_br(r) ,       c_γ(r) ≡ c₀[1+δ(r)] ,
∇·(ρ_br V) = −j_n .
```

**Source/domain:** v9:24–45,185–215, supplied on the live steady S9b object. The uniform anchor is
`ρ_br⁰ω² = μ_⊥k²`, recorded at `steps/S11b_interface_coupling_law.md:74–87`; its basis qualifications
do not define a varying-coefficient stress. `c₀` is the supplied asymptotic speed (v9:62–64).
`j_n` uses v9's normal-exchange source convention; the finite-slab identification is not thereby earned.
**Owner:** stiffness/inertia inputs S8; `j_n` and bulk profile gravity sector or S12 (v9:43–45,143–146).

**Status: supplied anchoring of the speed profile.** LAB_HELD keeps the steady `c_γ` profile at spatial
positions (v9:51–54,268–269). The source anchoring maps are

```text
Q_bg^L(x,t) = Q_bg(x) ,       Q_bg^M(x,t) = Q_bg(χ(x,t)) .
```

**Source/domain:** `directives/S11c_a_SHARED_PHYSICS.md:232–244`; these are distinct physical anchorings
on its supplied background. v9 selects the first for its speed. Its inverse material map `χ(x,t)` is
not the order field `χ_B`. This selection does not impose a material-constancy law or select a physical
holder. **Owner:** supplied S9b optical setting; material evolution remains §4.

**Status: supplied bulk-density inputs.**

```text
P = Kρ^n ,       c_s² = nKρ^(n−1)/m ,       f(r) ≡ ρ(r)/ρ₀ − 1 .
```

**Source/domain:** v9:140–146, with symbolic EOS exponent `n`, bulk **number** density `ρ` and particle
mass `m`; Part C's response is first order in `f`. Historical action content is postulated and its
foundation checks dimensional
(v2:`stages/ledger_stage004_gnls_action_dimensional_foundation.md:56–69,94–109`). This supplies neither
the bulk profile nor a brane-density response. **Owner:** bulk profile gravity sector or S12; substrate
EOS/action S1/S1.5. Part C's three speed-response choices (v9:164–172) are later compatibility inputs;
none is selected here as a constitutive closure.

**Status: OPEN `ℳ_⊥` (O1) and `ℛ_br` (O7).** These name the stiffness and brane-density responses,
respectively. Dependence on bulk state, flow, embedding, thickness and projection remains general,
as do their gradients. The speed ratio and bulk EOS determine neither response separately
(v9:283–286,309–312; `O2_premise_decision_list.md:45`). **Owner:** S8's brane inputs; unresolved material
reduction/completion stays with its substrate and S22 owners, rather than becoming a task here.

**Status: supplied slab/projection kinematics on the recorded finite-slab domain.**

```text
Σ_E = ρ_4D W ,      Σ_mat(X,t) = Σ_E(x(X,t),t) 𝒥_x(X,t) ,
δ_vΣ_mat = 0 ,
J_s^α = ρ_m(v_bulk,s − v_face,s^α)·n̂_s^α ,
a_s^α = sqrt(1+|∇_x h_s^α|²) ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋)                       (flat faces) ,
∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α ,         v = ∂_t u .
```

**Source/domain:** `directives/S11b_SHARED_PHYSICS.md:324–337` and
`directives/S11c_a_SHARED_PHYSICS.md:126–148,303–331,340–371`. Here `ρ_4D` and `ρ_m` are mass densities,
`α` labels the recorded anchoring, `h_s^α` is a face graph distinct from the reduced `h`, and `J_s^α`
is outward relative mass flux per true face area. The shape-derivative verification is confined to
supplied profiles and first shape order (`steps/S11c_a_interface_shape_derivatives.md:26–34,63–68`),
with background current/exchange frozen (`directives/S11c_a_SHARED_PHYSICS.md:369–390`). This domain
does not verify a live O6 map. The virtual constraint is separate from physical sourced evolution;
the displacement-model `v=∂_tu` is not an identification of steady `V` with light's wave velocity.

The historical factorization `ρ_br=ρ_4D W` is slab kinematics. RHO4-CONSTANT and RHOBR-CONSTANT select,
respectively, constant `ρ_4D,bg⁰` or constant `ρ_br,bg⁰` while thickness varies
(`directives/S11c_a_SHARED_PHYSICS.md:210–230`). Neither representative is adopted for O2 and neither
supplies `ℛ_br`.

**Status: OPEN `𝒥_map` (O6).** The material/order weighting, live projection/window and measure
identifications connecting sheet `j_n`, slab fluxes and bulk-normal `v_dr`, with bulk/return data, remain
general (v9:305–308; `O2_premise_decision_list.md:50`). The window in
`directives/S11c_a_SHARED_PHYSICS.md:394–402` is tied to its own face maps; it is not a newly selected
O2 projection. Stage 006's shear projection

```text
μ_R = ∫ dw χ_B μ_R⁽⁴⁾
```

is **postulated/PENDING**, with dimensional consistency asserted only
(v2:`stages/ledger_stage006_two_phase_chiB_ontology.md:98`); it is not a supplied `ℳ_⊥` law or a complete
live map. **Owner:** S12 for dynamical conversion and separate source/boundary inventories
(`V3_STEP_PLAN.md:579–617`; `steps/S11c_PARTIAL_CLOSEOUT.md:21,33`); S14a retains the distinct
projected order-loss/far-field flux bridge (`V3_STEP_PLAN.md:626–641`). Sub-step 3 must state its use of
O6 in the balance; this contract does not perform that reduction.

## 7. Bulk-traction premises and normal-load operand

**Status: postulated shear-free scalar bulk, retained by adopted premise 4 (user, 2026-10-06);
Label: adopted substrate input to a conditional model (2026-10-06).** Underlying status is explicitly
postulated at `steps/S9_light_requires_shear.md:180–187,331–337`. The supplied rest-frame acoustic model is

```text
v_bulk = ∇₄φ ,       δp = −ρ_m ∂_tφ ,       ∂_t²φ = c_s0² ∇₄²φ .
```

**Source/domain:** `directives/S11b_SHARED_PHYSICS.md:162–181`, with no bulk shear modulus and outgoing
or decaying radiation conditions at `:89–93`. Its operators are linearized about rest, with active
`v_dr` excluded (`:99–111`). The uncarried background-normal-flow correction remains a recorded scope
limit (`steps/S11b_interface_coupling_law.md:154–164`); this rest model supplies no live background
pressure response. **Owner:** inherited bulk premise from S9/S11b; live drain/boundary inputs S12.

Premise 4 supplies only the directional restriction on bulk mechanical traction:

```text
t_bulk,s^live = 𝒯_bulk,n,s^live n̂_s .
```

**Status of amplitude: OPEN `𝒯_bulk,n,s^live`.** The scalar is a general signed normal-load operand,
not an adopted pressure/affinity or DC response law. **Source:**
`O2_premise_decision_list.md:37–39,48`. Face geometry and projection stay live; loading normal to a
tilted face must not be replaced by zero in-plane projection on the brane. There is no independent
tangential bulk stress. **Owner:** O2 carries the OPEN load; S12 owns needed live bulk/drain data. No
new live traction law is assigned to S11c.

**Status: supplied perturbation traction and face work, reported on their original domain.**

```text
𝒜_s = μ_θ/ρ_br⁰ − δp_s/ρ_m ,
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,       Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,
t_s = −(δp_s+Λ_X(ω)𝒜_s)n̂_s ,
δ_v𝒲_bulk^α = Σ_s a_s^α t_s^α·δ_vx_s^α .
```

**Source/domain:** `directives/S11b_SHARED_PHYSICS.md:194–224,358–372` and
`directives/S11c_a_SHARED_PHYSICS.md:348–371`. These are prescribed linear responses with independent
real response constants/times; the affinity and virtual/mass accounting have independent reviewer
derivations recorded in `steps/S11bB_interface_assembly.md:126–139`. The nonuniform extension is first
shape order with its frozen background. These kernels are not premise 1's material-reference law;
putting `ω=0` does not supply `𝒯_bulk,n,s^live` or a background `j_n` law.

**Status: OPEN boundary/support inputs.** The complete `T_hold,s` remains OPEN. Its bulk component obeys
premise 4, whereas a declared external support would be a **supplied held input**, not a computed holder
(`O2_premise_decision_list.md:48–49`). The recorded support bundle
`𝒮_hold⁰={f_hold⁰,t_hold,s⁰}` is supplied with `V_s⁰=J_s⁰=𝒜_s⁰=0`
(`directives/S11c_a_SHARED_PHYSICS.md:246–279`); its energy/geometry comparison is not the live O2 law
(`directives/S11c_b_SHARED_PHYSICS.md:190–235`; v9:273–279).
**Owner:** full traction/support partition and force balance sub-step 3. Physical core holder and mouth
data `ℋ_core` (O5) remain with Q2/S22 (`V3_STEP_PLAN.md:896–939`; v9:301–304). Neither a held datum nor
LAB_HELD anchoring supplies that physical response.

## 8. Steady-state energy input

**Status: adopted premise 1's energy-accounting requirement (user, 2026-10-06);
Label: adopted substrate input to a conditional model (2026-10-06).**
**Form: OPEN `ℬ_E^steady`.** This denotes the steady state's energy balance, not a computed residual.
It explicitly carries the following general OPEN content:

| Named energy operand | Required physical content / owner |
| --- | --- |
| `ℰ_br^live`, `𝒥_E^live` | Material energy storage and transport compatible with the live stress, inertia, normal response and reference evolution. Conservative antecedents S1.5; brane material content S8, nonlinear completion S22. |
| `𝒫_ref/relax^live` | Power associated with the reference/strain evolution and relaxation response `ℛ_ref/strain^live`. O2 carries it under premise 1; no formula, sign or vanishing is supplied. |
| `𝒫_convert/exchange^live` | Order conversion, energy carried with exchanged material and any additional non-variational energy partners. S12 owns their forms and reaction/supply system; premise 3's momentum choice does not determine them. |
| `𝒫_boundary^live` | Mechanical work and energy transfer through the live bulk/face/boundary data, and any explicitly declared held support. O2 carries the accounting; its force/traction pairing is for sub-step 3, while native drain/return data belong to S12 and physical core response to Q2/S22. |
| `𝒮_E,net`, `𝒫_E,supply` | The physical supplier of any net power and its stated budget. Both identity and budget remain OPEN; naming the drain as drive supplies neither its available energy nor a numerical or functional power budget. |

These names identify accounting obligations, not an assumed additive constitutive decomposition or
independent channels. `ℬ_E^steady` must retain the relaxation power explicitly and identify the supplier
of any net power when a conditional closure is stated (`O2_premise_decision_list.md:22–25`). No supplier
mechanism, passive sign, energy-reference reset, or energy-free renewal is adopted. **Owner:** O2 input
requirement and sub-step 3's coupled accounting; underlying balances/partners retain the owners above.

**Status: postulated historical order-work contribution, qualified separately.** Stage 006 records

```text
μ_χ = δF/δχ_B ,       P_order = ∫ d⁴X μ_χ D_tχ_B .
```

**Source/domain:** v2:`stages/ledger_stage006_two_phase_chiB_ontology.md:69–75,100`, within its
postulated real-fraction free-energy/dynamics adjunct. This contribution uses no extra number-density
factor. It supplies neither the branch-dependent S12 completion nor `𝒫_ref/relax^live` nor the full
`ℬ_E^steady`. The paused comparison explicitly does not establish re-ordering as an available energy
source (`docs/elastic_reference_comparison_assessment.md:43–49`).

**Status: recorded energy-reference question; OPEN choice `C_ref`.** The historical action writes
`U(ρ)=Kρ⁵/4` as postulated content (stage 004 `:63–69`). S1.5 instead records, for that `n=5` EOS,

```text
P = ρU′ − U = Kρ⁵ ,       U(ρ) = Kρ⁵/4 + C_ref ρ .
```

**Source/domain:** `V3_STEP_PLAN.md:188–204`; `C_ref` renames the plan's chemical-potential/energy-reference
`C`, distinct from S11b's density–thickness coupling. The EOS leaves this reference choice unresolved.
**Owner:** S1.5, together with its momentum-stress/quantum-energy/current improvement convention.
This contract chooses neither `C_ref=0` nor an extension of this `n=5` expression to v9's symbolic `n`.

The existing non-passive-interface condition remains separate: if a later model adopts such a response,
it must name a reservoir and state a power budget (`steps/S11b_interface_coupling_law.md:57–63`;
`steps/S11bB_interface_assembly.md:195–197`). This does not select that response or settle relaxation
power. The inventory `§4` records that the historical S11b-C/S11c handoff has no live successor owner
for this condition; `SUBSTRATE_REQUIREMENTS.md:402–420` targets S12 only by register inference, with
review pending at `:3`. O2 must carry the condition if used, separately from S12's source partners.

## 9. Recorded counting and OPEN grades

**Status: supplied v9 counting.**

```text
ε(r) ≡ GM/(c₀²r) ,       δ = O(ε) ,       (∂ξ_w)² = O(ε) ,
V/c₀ = O(ε^{1/2}) ,      (V/c₀)² = O(ε) ,
δ^a (V/c₀)^b ((∂ξ_w)²)^c ,       0≤a≤1 , 0≤b≤2 , 0≤c≤1 .
```

**Source/domain:** v9:70–90,314–323; the box is the optical retained monomial set, with nonnegative
integer indices. `GM` is the independent slow-test-matter orbital parameter (v9:132–134), not any
profile, source or mouth amplitude. The ratio `μ_⊥/ρ_br` alone inherits the speed-change grade; the mass
balance retains its full density factor and derivatives. The induced-metric qualification is the
recorded relative `O(ε)` correction to `j_n`, not a new projection law. Fixed `ℓ` maps live `h` slopes
to the geometric grade without assigning a mouth-amplitude grade. First order in bulk `f` remains
separate; no relation between `f` and `ε` is supplied. **Owner:** gravity sector or S12 for v9's counting;
the orbital matching interface stays with the gravity work. Inventory §4 identifies S16's response-side
role as an inference, with its calibration qualifications intact.

**Status: OPEN missing grades and derivative scales.** Individual density, stiffness, inertia,
stress/normal-response, relaxation/power, source, exchange-momentum, force, traction/support, holder,
embedding-sector coefficient/source and longitudinal-field grades remain named unknowns wherever not
recorded. No derivative scale or further grade is chosen to remove an operand. This implements the
decision list `:51` and v9:309–323; the O2 deliverable remains untruncated. S11c's supplied-profile
bookkeepers and coefficient freezes are not an O2 order contract.

## 10. Remaining interfaces

S1.5/S8's missing conservative/material content stays with those owners or as the operands in §3;
reference evolution and the energy forms remain OPEN under the adopted character of premise 1.
O4's relation to O2's normal content is unsettled; `ℳ_⊥`, `ℛ_br`, `𝒥_map`, live normal loading,
`T_hold,s`, O5 data and missing grades remain OPEN. No additional physical choice is required merely to
carry them in the conditional named balance, and this contract adds no restrictive law.

Sub-step 3 owns drive representation, O3 momentum transport, the full face/support balance and their
energy pairing. S12 owns non-variational partners and separate source/boundary data; S14a owns the
dynamical-drain/far-field bridge and S14 remains conditional on it (`V3_STEP_PLAN.md:196–204,579–646`).
S16's orbital response interface retains its conditional/calibrated domain (`:670–698`); no drive/`GM`
coupling is supplied here. Q2/S22 retain holder and nonlinear material questions. S21 owns integration
and the revising-input/new-consequence sort (`:1156–1168`); premises 1–4 retain their adopted-input labels
at that handoff. Optical compatibility, S11c repair/composition and prior-art comparison supply no
premise or outcome in this contract.
