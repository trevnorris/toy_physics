# Measurements: S15 decision list, setup lookups

Generator: `_scratch/s15/gen/s15_setup_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval by line, literal-match counts.

## L0. The list

```
$ sha256sum research/pde_ledger_v3/directives/S15_decision_list.md && git log -1 --format='%h %s'
d83e69e5046ec71be88acd240b4167971b466a2d07ed5f52936f4bf6d4239a01  research/pde_ledger_v3/directives/S15_decision_list.md
d6e4f74a Preserve the S15 decision list as reviewed in pass 1 (not accepted)
```

## L1. The plan: S15, the S14a conditional, the S22 rows, S19's heading

```
$ sed -n '652,669p' research/pde_ledger_v3/V3_STEP_PLAN.md
### S15 · Two drains attract ⚠ RUN S16 FIRST for regime context — ⛔ NOT as S15's premise
⛔ **CONDITIONAL ON S14a.** Not a carried result until the dynamical-`Γ_B` bridge succeeds.
⚠ **Corrected twice.** Round 1: run S16 first, since the corpus establishes the worldtube closure before
the particle law. ⛔ **Round 2 corrects that correction: S16 does NOT license S15.** S16 is
*response-side* (how a body moves in a field); S15 is *source-side* (what force two drains exert) and
carries **its own** Noether / control-surface assumptions. ⇒ Run S16 first for regime context, but
⛔ derive S15's premises independently — do not inherit them.

Noether stress at infinity, drain strengths `Q_i` **given**: `F_12 = −I_F,12 · (m N_∞,3 Q1Q2 / 4πr²) r̂`  ⛔ **`I_F,12` is not 1** — it is listed as debt in S22; ⛔ do not display the coefficient-one form.

⛔ **Two corrections to the original wording:**
- **The source map is missing.** S13 introduces `J`; this step uses `Q_i` — the corpus relation
  `Q_i = Θ_Qi J_i / N_∞,3` must be an explicit step, not an unremarked substitution. Add
  `{J, Θ_Q, N_∞,3, Q}` to S22.
- **The sign claim is overstated.** This is the **leading matter-stress** contribution only; the source
  says the full sign retains quantum, confinement, Maxwell and profile residuals, and the normalization
  is **calibrated**. ⇒ Bank it as *"leading matter-stress sign"*, ⛔ not *"the force sign, earned"*.

```

```
$ sed -n '639p;1305,1307p' research/pde_ledger_v3/V3_STEP_PLAN.md
⛔ **Until this succeeds, S14, S15 and S20 are CONDITIONAL, not carried results.** ⚠ Mark them so in
| `M` and the supplied multipoles | inputs of the S16 reduction, not outputs |
| `G` / `N_inf,3` normalization, `I_F,12` | the force magnitude is calibrated through them |
| `Theta_Q` and the `J -> Q` map | S15 must not substitute silently |
```

```
$ sed -n '766p' research/pde_ledger_v3/V3_STEP_PLAN.md
### S19 · Gravitomagnetism ⛔ NOT a cite-and-record step — it has a known open problem
```

## L2. The flow framing, R-S21-01 and R-S14a-01

```
$ sed -n '759p;765,766p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  record M01/M18): each force states the flow it needs, S21 compares them, and a mismatch
- **Flow comparison decision, user 2026-10-19** — Each force states the flow it needs; S21
  compares those requirements, and a mismatch reopens S9b. S12 states the source compatibility
```

```
$ sed -n '714p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
### R-S21-01 — conditional deflection and radar-slope profile requirements
```

```
$ sed -n '974,986p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
### R-S14a-01 — the dynamically selected slow limit and far behaviour

- **source** S14a part 2 (`steps/S14a_2_slow_limit.md`, §§1–5/7; retrievals M01–M13/M16/M17),
  for the conditional consumers S14, S15 and S20.
- **target** jointly S12 (source and flux forms), S8/S22 (the slab's material law), and
  **no owner named** for the face and interior laws, by the user's V2 decision at `12aa65d5`.
- **status** OPEN
- **requirement** — the slow limit and far behaviour selected by the slab's own dynamics,
  for S14a's still-underived controlled return and profile-dependent `J→Q` obligations.
  The selected steady-throughput and exchange-live decisions keep spreading, exchange and
  storage live. No average, solution family, density gate, order of limits, law, material,
  profile, premise or candidate is selected here.
- **on failure** — S14, S15 and S20 remain **CONDITIONAL, not carried results**.
```

## L3. P5 and P6: the S9b record rows and the S9b spec

```
$ sed -n '667,668p' research/pde_ledger_v3/steps/S9b_light_bending.md
| P5, D2, 2026-10-07; carriage 2026-10-08 | Momentum density ρbrV and carried current ρbrVⁱVʲ. Part D only; no A–C object. Additional stress-carried momentum/retained grade remains OPEN, S8. |
| P6, D2, revised 2026-10-08 | Isotropic pressure plus linear viscous Cauchy stress; live pbr(ρbr), ccomp²=dpbr/dρbr, η, ζ. Part D only; no A–C object. Not derived from P1; power-law creep excluded from that adopted form, relaxation and supplier/power identification OPEN. |
```

```
$ sed -n '340,378p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md

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
```

## L4. The projected domain, S16's general balance and the core

```
$ sed -n '62,64p;76,86p;89,90p' research/pde_ledger_v3/directives/S16_decision_list.md
## Supplied equations (unfalsifiable in this build)
On S14a part 1's projected domain: in-plane coordinates `x^i` (three) and time `t`, with the slab integrated in `w`
between its live faces.
- **The general projected momentum balance** (the user's "General balance"):

  ```
  ∂_t G^i + ∂_j Π^{ij} = f^i
  ```

  - `G^i` is the projected momentum density, `Π^{ij}` its flux and `f^i` the force density. All three are general
    live fields of `x` and `t`.
  - `f^i` includes every exchange of momentum with the faces and the bulk. The projected inputs do not divide it by
    face or by the faces' motion, so it stays undivided (E2).
  - `G^i` is not related to `F^i` or `Σ` by any supplied law, and `Π^{ij}` is not assumed symmetric.
- **The throat's centre** is a live point `X_c(t)`. The excluded core lies inside an inner surface that moves with
  it. The core's content stays with S22.
```

## L5. O2: the drain premise, the historical matter stress, the transfer limit

```
$ sed -n '94p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-400
| **2 — drive** | The committed dynamical order-conversion drain is the drive. No separate external body force. Its source and boundary/return data act through stress, traction and O3; the `GM` bridge is OPEN. | Use the permitted representation in which `F_drive` is **absent as a separate body-force entry**. Declare the drain-drive provenance of the material, boundary and exchange entries. Do no
```

```
$ sed -n '250,257p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
The historical matter-stress current `Π_ij` includes the convective term `mρv_i v_j`, while
`σ^Q_ij` is its recorded quantum-stress content (C §7, stage 002 record). Its `v` is the native steady
radial bulk/reduced-lane inflow, with no supplied identification to `V`, `V_s` or `v_dr`. Convective
material transport belongs to the native momentum-current accounting; it is not an extra normal
mechanical pressure to add to `T_hold,s`. Its transfer into the O3 description requires the OPEN
material/map and S12 partner identifications. Do not use both a mapped convective current and the
same carried-material O3 current as separate transfers. Premise 4 excludes tangential mechanical
loading from `σ^Q_ij` as an adopted restriction, rather than a consequence of rest acoustics.
```

```
$ sed -n '351p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-400
| Far field, spherical isolated mass at rest, Eulerian steady profiles, lab time; linear optical waves and leading eikonal | Supplied model setting. No transfer to strong field, mouth/interior, moving or rotating mass, or a drain varying during a light crossing (C §1; I §1; v9 model point). |
```

## L6. Prior art: the gravity section's rows G7, G9 and G13

```
$ sed -n '417p;419p;423p' docs/medium_requirements_and_prior_art.md | cut -c1-700
| `G7` the force a drain feels from the flow around it. ⛔ Not S16's object: S16 computes no passive factor, and files as requirements a named field of the surroundings, a rule separating it from the body's own field, and its normalization. This prior result is what the record will check that open requirement against | **KNOWN** | ⭐⭐ **Lagally**, *ZAMM* **2** 409 (1922): for **steady** flow, the force and moment on a closed surface depend only on the sources and vortices and on the velocity at their location, and on the boundary's shape only through that velocity. Unsteady flow and moving bodies: **Landweber & Yih**, *JFM* **1** 319 (1956). Multipoles and deformable bodies: **Landweber 
| `G9` two drains attracting. ⛔ S16 gives no force law; this is an oracle for a later step | **KNOWN** | ⭐ **Lamb**, *Hydrodynamics* (1906) Art. 144, p. 192: *"two sources of like sign attract, and two of unlike sign repel, with forces varying inversely as the square of the distance"*. Also **Wang**, Corollary 4: the force on a sink of strength `Q₂` moving at `v₂` is the inverse-square interaction term `−ρQ₂Q₁/(4πr²) r̂₂₁` **plus a velocity term `ρQ₂v₂`** (`v₂` *"is the velocity of the source with strength Q2"*), so the inverse square alone holds for a sink at rest. Pulsating bodies: **V. Bjerknes**, *Fields of Force* (1906) pp. 32–33, after **C. A. Bjerknes**.
| `G13` objections to the sink ether | **KNOWN** | A moving drain feels a drag: Wang's formula (32) gives *"a universal damping force"*. Energy disappears without explanation, some drag must arise, and no process that creates matter is known (Wikipedia, secondary only). For pulsation, every pulsation in the universe would have to be in phase (same source) |
```

## L7. The earlier construction of the force (earlier material, not a premise)

```
$ sed -n '44,53p' software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md
Convention: n_hat outward from U2; F_12 is force on defect 2 by defect 1; stationary F_12=-int_boundary Pi_ij n_j dS+int_U2 f_i^body dV.

Control surface: reduced-3D sphere around defect 2 with core scale << a << r12; v2=-Q2*n_hat/(4*pi*a^2), v1 is constant over the sphere to leading order.

Reduced-3D compact lane integral results:

- Convective cross flux: `int Pi_conv,cross.n dS=-(4/3)*m_GNLS*N_infty,3*Q2*v1`.
- Bernoulli/EOS pressure cross flux: `delta P_cross=-m_GNLS*N_infty,3*(v1.v2), so int delta_ij P_cross n_j dS=+(1/3)*m_GNLS*N_infty,3*Q2*v1`.
- Matter cross flux: `int Pi_matter,cross.n dS=-m_GNLS*N_infty,3*Q2*v1`.
- Force result: `F_12^matter=m_GNLS*N_infty,3*Q2*v1=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*rhat_12`.
```

```
$ sed -n '34p' software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md
Force form: `F_12=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*I_F,12*rhat_12 in the compact reduced-3D lane, with Qi=Theta_Qi*Ji/N_infty,3.`
```

```
$ sed -n '141,144p;201,202p' research/pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md
Take a control sphere of radius `a` (core scale `<< a <<` separation) around
defect 2, with defect 2's own inflow `v_2 = -Q_2 n_hat/(Omega_{d-1} a^{d-1})`
and defect 1's inflow `v_1` approximately constant across it. Let `N` be the
ambient number density. The traction is `Pi_ij n_j`; integrate over the sphere.
reduced-3D (Omega_2 = 4 pi):   F_12 = - m N_3   Q_1 Q_2 / (4 pi   r_12^2) r_hat_12,
bulk-4D    (Omega_3 = 2 pi^2): F_12 = - m rho_4 Q_1 Q_2 / (2 pi^2 R_12^3) R_hat_12.
```

## L8. The charter's conditional, C3, the slab material owner, R-S16-03

```
$ sed -n '238,241p' research/pde_ledger_v3/CHARTER.md
⛔ **CONDITIONAL — neither carried nor discarded, pending S14a:** {#conditional-s14a}
- the `1/r²` law, the attractive sign, and the stage009 `RETURN_RESIDUAL_PREDICTION` falsifier.
  Correcting the drain to the dynamical `Γ_B` law severed the chain to these — they were derived under a
  Gauss **number-flux** drain. ⛔ Do not cite them as earned until the S14a bridge succeeds.
```

```
$ sed -n '171p' research/pde_ledger_v3/DEFECT_REGISTER.md | cut -c1-300
| **C3** | **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING` / `BLOCKS_MASS_EMERGENCE`. Only a dimensional bridge `α_J ħJ/c_γ²` exists; pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`, *"No action-level, boundary-source, Noether, or Hamiltonian equation"* | `stage004:181`; `parameter
```

```
$ sed -n '292p' research/pde_ledger_v3/steps/S12_order_source.md | cut -c1-300
| Interior-order law/history/strain; faces' law; source/flux forms | No separate owner named for the first two; source/flux forms remain S12's. Slab material completion remains separately S8/S22. None of the forms is delivered here (M01, M11, M18). |
```

```
$ /usr/bin/grep -n -E '^### R-S16-03' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md && sed -n "$(/usr/bin/grep -n -E '^### R-S16-03' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -d: -f1),+8p" research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | /usr/bin/grep -E 'target|status'
1054:### R-S16-03 — momentum exchange at each slab face
- **target** no owner named
- **status** OPEN
```

## L9. Literal-match counts in the list (leak probe; each hit is read by the orchestrator)

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/S15_decision_list.md)"
14:- **The plan's step.** S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
count: 1
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/S15_decision_list.md)"
146:    printed pieces, a stated dependence on the separation, including which give an inverse-square one, as the user's
count: 1
```

```
$ /usr/bin/grep -n -F -- '1/r' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- '1/r' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r²' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'r²' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r^2' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'r^2' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'drag' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'drag' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'sign' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'sign' research/pde_ledger_v3/directives/S15_decision_list.md)"
21:  - limits the sign claim to "the **leading matter-stress** contribution only" (`:666–668`).
150:  - the sign of a printed piece only on its stated conditions, and only for that piece (`:666–668`).
173:- **No owner assigned here:** the forms of `f^i` and `Π_o^{ij}`. The slab's material law goes to S8/S22
179:- **The pieces of the full sign that P5 and P6 do not supply:** the plan's quantum, confinement and Maxwell residuals
count: 4
```

```
$ /usr/bin/grep -n -F -- 'Q1' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Q1' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q_1' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Q_1' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q₁' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Q₁' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '4π' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- '4π' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'repel' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'repel' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Bjerknes' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Bjerknes' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lamb' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Lamb' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Newton' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'Newton' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'damping' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'damping' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'positive' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'positive' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'negative' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F -- 'negative' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```
