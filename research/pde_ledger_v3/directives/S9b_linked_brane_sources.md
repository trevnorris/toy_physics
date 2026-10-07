# S9b v9 — linked-brane sources and OPEN premises

**Deliverable:** the v9 addition to `S9b_SHARED_PHYSICS.md` supplies the recorded links and their domains for
Part D's Part B/Part C conditions and implied normal exchange, with `V` and `j_n` live.
**Status:** authored for review, 2026-10-06; no review or build clearance is claimed. The v8 baseline is
`c2f1cf2b`. This note supplies governing inputs, not solved profiles or light-observable results.

Paths below are relative to `research/pde_ledger_v3/`, except the two embedding records, which are under
`research/pde_ledger_v2/notes/stages/`. Line ranges refer to the source files read for this authoring.
References to `S9b_SHARED_PHYSICS.md` use the v8 baseline at `c2f1cf2b`, before the v9 addition changes its
line numbering. “Stage 030” and “stage 031” name the two embedding records listed in the directive.

## Supplied relations

| ID | Equation supplied | Source lines and applicability |
| --- | --- | --- |
| L1 | `c_γ² ≡ μ_⊥/ρ_br`, `c_γ ≡ c₀(1+δ)` | `directives/S9b_SHARED_PHYSICS.md:24–33,70–82` supplies the pointwise identification and counting. `steps/S11b_interface_coupling_law.md:74–87` supplies the uniform transverse anchor and basis qualification. This is not a constitutive density law. `directives/S11c_b_SHARED_PHYSICS.md:165–174` prohibits lifting the uniform coefficient fold into a variable-coefficient identity. |
| L2 | `∇·(ρ_br V)=−j_n` | `directives/S9b_SHARED_PHYSICS.md:37–45` already supplies the live steady balance. `directives/S11c_a_SHARED_PHYSICS.md:126–148` supplies `Σ=ρ_4D W`, `∂_tΣ+∇_x·(Σv)=−(J₊+J₋)` and `v=∂_tu`, and separates physical evolution from the virtual constraint. Its `:348–371` supplies `J_s^α=ρ_m(v_bulk,s−v_face,s^α)·n̂_s^α`, the true-face-area convention, and `∂_tΣ^α+∇_x·(Σ^αv)=−Σ_s a_s^αJ_s^α`. The record's linearized background has `J_s⁰=0`; that freeze is not inherited by the steady S9b balance. |
| L3 | `ξ_w=ℓh`, `h=P₀H=N₀⁻¹∫dw 2f₀H`, `f₀=1/[ℓcosh²(w/ℓ)]`, `N₀=∫dw 2f₀²`, `M_h=N₀M₄`, `K_h=N₀K₄=M_hc_E²` | Stage 030 `:46–51,63–66,88–112,127–130` supplies the parent, mode, projection and reduction. Stage 031 `:60–76` identifies the normal geometric displacement with this dimensionless reduced scalar. Stage 030 `:15–22,228–233` keeps the entire closure postulated; `:259–261` supplies no cone lock. `directives/S11b_SHARED_PHYSICS.md:95–96` distinguishes its face-centre perturbation `ζ_c` from this `h` field. |
| L4 | `A_eff=ρ_br+C_J²/κ_phase`; `S_Lh=∫dt d³x[½A_eff(∂_tu_L)²+½M_h(∂_th)²−½B_eff|∇u_L|²−½K_h|∇h|²−C_hu∇u_L·∇h]` | Stage 030 `:134–142` supplies the reduced coupled action. Its `:15–22,228–233` and `V3_STEP_PLAN.md:871–893` qualify it as the static electric sector, earned given a postulated G0 action, with C6 unresolved. It has constant coefficients and no live background flow or normal-exchange dynamics. It is supplied as a conditional sector input, not as the flowing-brane balance. |
| L5 | `(δΩ/δh)_mouth=η_i(k_mh−g_χhs_i)`, `k_m=K_m/ℓ²`, `g_χh=J_m/ℓ`, `Q_χ[r_Σ,s_i]=s_i` | Stage 031 `:68–76,96–138,147–155` supplies the coupling reduction, orientation projection and mouth Euler density. Its `:27–30,243–248` retains the postulated frozen sleeve/profile class. `s_i` is orientation, not Part C's exponent `s`; `J_m` is a coupling, not `j_n`. No identification with a mass-driven drain is supplied. |
| L6 | `d/dr(r²dh/dr)=0`, `h(a)=h_A`, `h_A=ξ_w|_A/ℓ=P₀H|_A`, `h→0` at infinity, for a source-free static exterior with constant generic `κ_ext>0` | Stage 031 `:157–174` supplies the radial governing equation, decaying boundary condition, held datum and unresolved core requirement. `κ_ext` is a rename of its generic exterior `κ` at `:164–165`, explicitly distinct from the response `κ=D/B_eff` at `:178–184`. `V3_STEP_PLAN.md:896–915` retains the static Q2 regime and holder debt. The solved exterior profile is not supplied to S9b. |

Only L1 and the already-supplied steady L2 are unconditional links of the S9b live-profile object. L3
identifies the embedding field within its postulated sector. L4–L6 supply that sector's governing inputs
with their static domains intact; their applicability to a mass with live flow is O4. No additional closed
steady force or embedding law with `V` and `j_n` live was supplied by the checked records.

## Reported equations and why they are not additional live-flow premises

- **Density representatives.** `ρ_br=ρ_4D W` is finite-slab kinematics
  (`directives/S11b_SHARED_PHYSICS.md:324–337`;
  `directives/S11c_a_SHARED_PHYSICS.md:132–148,210–230`). The latter supplies
  `ρ_4D,bg⁰=ρ_4D,ref⁰`, `ρ_br,bg⁰=ρ_4D,bg⁰W_bg` in RHO4-CONSTANT and
  `ρ_br,bg⁰=rho_br`, `ρ_4D,bg⁰=rho_br/W_bg` in RHOBR-CONSTANT. Those selected freezes are reported,
  not adopted; neither supplies `ρ_br(f)` for Part C.
- **Anchoring.** `Q_bg^L(x,t)=Q_bg(x)`, `Q_bg^M(x,t)=Q_bg(χ(x,t))`
  (`directives/S11c_a_SHARED_PHYSICS.md:232–244`) are distinct physical anchorings. S9b v8's
  `directives/S9b_SHARED_PHYSICS.md:51–54` selects LAB_HELD with live radial flow and varying speed.
  No material-constancy condition is added to that choice.
- **Face response.** `J_s=Λ_A(ω)𝒜_s+Λ_V(ω)V_s`,
  `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)`, `𝒜_s=μ_θ/ρ_br⁰−δp_s/ρ_m`,
  `μ_θ=(δU/δθ)|_{u,e_W,others fixed}`
  (`directives/S11b_SHARED_PHYSICS.md:194–224`;
  `directives/S11c_a_SHARED_PHYSICS.md:348–355`) are linear perturbation interface laws. The face normal
  velocity `V_s` is distinct from S9b's in-plane `V`. Substituting a DC frequency does not supply a
  background exchange law or the missing static chemical potential.
- **Wave momentum and exchange force.** `directives/S11b_SHARED_PHYSICS.md:353` names momentum balance
  for the wave displacement `u`; `:99–111` keeps `v_dr` out of its operators. Its `:356–363` supplies
  separate sourced evolution and `Q_J^direct=0` within that perturbation model. It supplies neither the
  live background balance nor its exchange momentum.
- **Frozen support balance.** `𝒮_hold⁰={f_hold⁰,t_hold,s⁰}` and
  `V_s⁰=J_s⁰=𝒜_s⁰=0`
  (`directives/S11c_a_SHARED_PHYSICS.md:246–279`;
  `directives/S11c_b_SHARED_PHYSICS.md:190–196`). The b admissibility operand is the background-order
  energy/geometry force-and-traction operand; its residual is operator operand minus declared support
  (`directives/S11c_b_SHARED_PHYSICS.md:215–235`). It is not a tested flowing-background balance.
- **Holder.** Stage 031 `:170–174` and `V3_STEP_PLAN.md:908–915` leave the physical holder unresolved.
  The plan's `:917–939` preserves competing candidates without selecting one. A held mouth datum does
  not close this debt.

## S11c step-record applicability check

The records were checked for an additional relation rather than importing an operator result as a
background law:

- `steps/S11c_a_interface_shape_derivatives.md:26–34,255–260`: geometric shape derivatives on supplied
  profiles; the normal-background-flow correction remains uncarried.
- `steps/S11c_b_variable_coefficient_operator.md:35–44,46–57`: displacement/constitutive/sourced-evolution
  operators on the supplied background; repaired equations have scoped review and unresolved composition.
- `steps/S11c_c1_curved_bulk_closure.md:132–155,185–196`: the density is live in the response; rest-frame
  and shape-order qualifications remain. The drain-tilt projection qualification is not a live-drain law.
- `steps/S11c_c2_self_energy_fold.md:192–195`: restore the spatial density dependence before the fold;
  no constitutive density response or background holder follows.
- `steps/S11c_d_profile_conditioned_scattering.md:25–29,60–72`: selected rest/supplied-profile results do
  not establish flowing-background behavior; live conversion goes to S12 and holder/response to Q2/S22.
- `steps/S11c_PARTIAL_CLOSEOUT.md:19–21,35–37`: current handoff retains those owners and treats exploratory
  support notes as inputs, not adopted laws or solved throats.
- `steps/S11c_SCOPE.md:39–43`: the normal background-flow correction remains uncarried; the frozen-width
  density mapping is an S11b modeling choice, not a flowing constitutive law.

## OPEN premises carried into Part D

| ID / live operand | What is missing | Source boundary |
| --- | --- | --- |
| O1 / `ℳ_⊥` | Constitutive dependence of `μ_⊥` on `ρ_br`, embedding, flow and bulk state; compatibility with Part C's imposed speed responses. | L1 is a ratio, not that law. `directives/S11c_a_SHARED_PHYSICS.md:171–208` supplies independent thickness/modulus profiles; `directives/S11c_b_SHARED_PHYSICS.md:165–174,199–206` supplies no variable-coefficient constitutive closure. |
| O2 / `ℬ_hold^live`, `F_drive`, `T_hold,s` | Steady in-plane and normal force/support balance with live flow, the force driving it, and the coupling to independent `GM`. | Frozen support sources above; `steps/S11c_PARTIAL_CLOSEOUT.md:19–21,35–37`. The wave momentum row is not the required relation. |
| O3 / `Π_n` | Momentum carried by normal exchange: vector direction and exchanged-material velocity, for the live steady balance. | `directives/S11b_SHARED_PHYSICS.md:356–363` supplies no direct generalized flux force in its linear perturbation scope; this cannot be promoted to a zero steady exchange-momentum premise. |
| O4 / `ℰ_h^live` | Applicability/extension of the postulated static coupled embedding closure to the mass with live `V` and `j_n`, variable coefficients and gradients; a map of static `u_L` to steady flow and of the mouth source to the drain. | Stage 030 `:15–22,134–142,228–233`; stage 031 `:40–58,157–174,243–253`; `V3_STEP_PLAN.md:871–915`. No flowing extension is recorded there. |
| O5 / `ℋ_core`, live boundary datum | Physical holder and boundary response selecting amplitude or mouth flux/Robin data, and their relation to `GM`. | Stage 031 `:159–174`; `V3_STEP_PLAN.md:908–939` (R1/R61, competing mechanisms remain unsettled). |
| O6 / `𝒥_map` | Sheet/finite-slab flux and measure identification, relation of `j_n` to bulk-normal `v_dr`, and bulk/drain/return boundary data. | `directives/S11b_SHARED_PHYSICS.md:99–111,194–224`; `directives/S11c_a_SHARED_PHYSICS.md:348–390`; `steps/S11c_d_profile_conditioned_scattering.md:64–72`. L2 fixes the supplied conservation law, not this response. |
| O7 / `ℛ_br` and symbolic grades | Brane-density response to bulk `f`, any thickness/projection relation, and separate density/stiffness/source/support/holder orders and derivative scales. No relation of `f` to `ε` is supplied. | `directives/S11c_a_SHARED_PHYSICS.md:210–230` only supplies selected density representatives; `directives/S9b_SHARED_PHYSICS.md:70–89,139–145` supplies the speed/slope/velocity counting and separate first order in `f`. |

**Order placement.** L1 enters the `δ=O(ε)` speed-ratio grade, without separate density/modulus grades.
L2 contains one live velocity factor, at grade `b=1`, with the full density derivative; the sheet's
induced-metric correction retains v8's relative `O(ε)` qualification. L3 maps slopes into the geometric
grade `(∂ξ_w)²=O(ε)`, without assigning an amplitude order. L4–L6 remain qualified static inputs; their
extension, source/support/holder orders and density response remain O1–O7 operands. None justifies
dropping a term, selecting a profile, changing v8's retained box, or claiming a Part B/Part C outcome.
