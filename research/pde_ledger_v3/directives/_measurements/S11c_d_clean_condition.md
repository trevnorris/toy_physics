# S11c_d_clean_condition.md (v4) — grounding commands (rule 2 / E1)

Mechanical lookups, run from the repo root on 2026-09-30, HEAD `c1e96e76`. Regenerated from the commands below;
nothing transcribed. Fences are longer than fences in quoted source text.

````
$ sed -n 14,17p research/pde_ledger_v3/CHARTER.md
- ⭐⭐ **v3 then TOOK a method change (2026-08-01, user decision): REQUIREMENTS-FIRST.** Each force
  sector states what it needs to survive. **Brane and bulk are defined LAST, at the knit**, by asking
  whether one medium can satisfy every requirement at once. **A no-go between requirements IS the
  falsification.**
````

````
$ sed -n 3p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
The current pilot answers a strict rest-bulk, LAB_HELD/RHO4_CONSTANT development-input question. It does not yet answer leakage in the calibrated, draining medium. These limitations were explicit in the governing scope, but should have been foregrounded in the pilot decision. No result can remove them through numerical convergence alone.
````

````
$ sed -n 11,13p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
The development input sets mu_R=1, rho_br=1 and c_s0=10 in the L_ref/T_ref unit frame. Under the bare S9/R4 definition c_gamma^2=mu_R/rho_br, this gives c_gamma/c_s=0.1, not 1. The step plan explicitly distinguishes the derived ratio definition from the calibrated, uncommitted equality lambda_gamma=1; the equality has not been imposed on these inputs.

There is an additional distinction: that bare coefficient formula must not be substituted for the actual full S11c-d transverse dispersion. The saved omega=3 LEFT incoming branch has normal momentum 2.439262183530097 and tangential norm squared 1/20. Its fixed-point phase speed omega/sqrt(k_normal^2+1/20) is 1.2247448713915874, or 0.12247448713915873 of c_s0. At the full-contrast RIGHT end the corresponding ratio is 0.12186666955535794; its zero-contrast ratio returns to the LEFT value. These are arithmetic readouts of saved full-pencil modes, not a new root solve or a demonstrated identification with the calibrated light cone. The extra retained elastic structure and coefficient conventions cannot be silently discarded.
````

````
$ sed -n 5p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
On the user's request to assess before continuing and stop for go/no-go, the existing central-balance-v2 scientific child was suspended in memory, not terminated or restarted. The pause receipt records PID 4097233, its process start identity, cgroup and exact command. It had completed two matrix groups after twenty restorations and saved 12,288 batches of the next group; no finite solve had been reached. The guard/supervisor remain active, no deadline was added, and pinned sources are unchanged. Resume requires the user's decision and rechecking that process identity. Do not launch another job alongside it.
````

````
$ sed -n 1107,1111p research/pde_ledger_v3/V3_STEP_PLAN.md
⛔ **Two different objects — round 2. Do not merge them:** the **ratio** `λγ = c_γ/c_s` is a *derived
definition* (registry `R3`); the **equality `λγ = 1`** is a separate *calibrated / uncommitted* cone
lock. Classifying the ratio says nothing about the lock.

⇒ ⛔ **Both classifications are fixed, not open:** the **ratio** `λγ = c_γ/c_s` is **`derived`** (it is `R3`); the **equality `λγ = 1`** is **`calibrated`/uncommitted**. ⛔ Do not re-open the ratio's class. Introduce `λγ` with provenance, confront
````

````
$ sed -n 1116,1126p research/pde_ledger_v3/V3_STEP_PLAN.md
⭐⭐ **`λγ = 1` is not a free landing — observation already constrains it.** In this model `c_s` is the
**gravity-change/phonon speed** and `c_γ` is the **light-cone speed**
(`research/pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md:75-77`):

> *"`c_s` is the phonon/gravity-change speed; `c_gamma` is the light-cone speed; `v_b` is the condensate
> flow."*

Gravitational-wave and electromagnetic arrival-time observations (**GW170817 / GRB 170817A**) constrain
those two speeds to agree to roughly **1 part in 10¹⁵**. ⚠ **The observational value IS known** — `λγ = 1`
to that precision is not in doubt as a *target*.

````

````
$ sed -n 35,37p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
density representatives ρ4D/ρbr, two anchorings LAB_HELD/MATERIAL_ADVECTED), S11c-b computes: (1) the **§3a energy
basis** — the O(3)-Kronecker field-bilinear invariant family, corrected to **40 = 10 uniform + 15 ∂W_bg-spurion + 15
∂μ_R,bg-spurion**; (2) the **variable-coefficient slab OPERATOR ROWS** — the equations of motion for U/θ/e_W with the
````

````
$ sed -n 90,91p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
`v_bulk_normal_0` (the bulk normal drain, `S11b_SHARED_PHYSICS.md:104`) is a scope-limit parameter, not an
active DOF, and appears in no derived operator (§0).
````

````
$ sed -n 95,97p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
`∂_t²φ=c_s0²∇₄²φ`; the current and conservation law `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; and the
dynamic, anchored slab window `Ω` supplied in S11c-a §3. S11c-b performs no curved-bulk response solve (§0).
````

````
$ sed -n 145,148p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
````

````
$ sed -n 343,354p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
```text
δ_vx_s^α ≡ δ_vR_s^α|_X ,              v_face,s^α ≡ ∂_tR_s^α|_X ,
V_s^α ≡ V_{n,s}^α ≡ n̂_s^α·v_face,s^α .
```

Use the same `V_s^α` in every object below. With all bulk quantities traced at `R_s^α`, define once

```text
J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α ,
n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m ,
𝒜_s^α ≡ μ_s^α − δp_s^α/ρ_m ,
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α)n̂_s^α ,
````

````
$ sed -n 365,366p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
δ_v𝒲_bulk^α ≡ Σ_s a_s^α t_s^α·δ_vx_s^α ,
∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α ,       v = ∂_t u .
````

````
$ sed -n 26p docs/toy_model_ontology_summary.md
The model is a classical \(4+1\)-dimensional, one-medium analog in which our three-dimensional space is an ordered, shear-supporting brane; particles are finite, oriented throats through it; a trapped transverse brane-shear standing mode helps hold each throat open; throat transport, order conversion, stress, return, and reservoir response are the proposed material mechanism underlying a calibrated gravity observable; and freely propagating light, electric coupling, magnetism, radiation, passive response, inertia, and possible cosmic expansion are target behaviors of the same complete material configuration.
````

````
$ sed -n 100p docs/toy_model_ontology_summary.md
The bulk is therefore an internal reservoir, not an external supply. In a closed-loop realization, localized throat drainage transfers material, momentum, and energy from ordered brane degrees of freedom into de-structured bulk degrees of freedom, while distributed return transfers material back into the ordered state. The bulk may store the corresponding response as compression, flow, internal energy, entropy, or other unresolved excitations of the same medium. Any sustained cycle must include the evolution of that reservoir rather than treating it as an inexhaustible battery.
````

````
$ sed -n 315p docs/toy_model_ontology_summary.md
The fields \(h_+\) and \(h_-\) describe the interfaces only while each remains a single-valued graph over the brane coordinates. They are useful background, far-field, and separated-interface collective coordinates, not a complete nonlinear throat topology. Overhangs, necks, tubes, folded sheets, reconnection, pinch-off, and topology change must be represented by the parent fields \(\chi_B(\mathbf x,w),n(\mathbf x,w),\theta(\mathbf x,w)\), the support fields, and core level sets. The earlier single coordinate \(h\) should be read only as schematic notation for a derived normal eigencombination.
````

````
$ sed -n 362p docs/toy_model_ontology_summary.md
- **Structural support:** energy \(E_{\rm support}\) in a trapped transverse brane-shear standing mode—the model's light or photon mode—helps hold the aperture open.
````

````
$ sed -n 957p docs/toy_model_ontology_summary.md
Freely propagating light and the throat-support mode must be distinguished. Free light is a background-brane \(\mathbf u_T\) excitation. The support mode is a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator; it is not a second photon substance and is not assumed to propagate as a shear wave through fully de-structured bulk material. Any \(\Omega_{\rm support}\) is a diagnostic energy/stress-localization region extracted from that solved mode, not a hard PDE domain imposed by \(\mu_R>0\).
````

````
$ sed -n 1366p docs/toy_model_ontology_summary.md
- \(\Gamma_{\rm drain}\) converts ordered brane material into de-structured bulk material at throats;
````

````
$ sed -n 113,116p docs/native_light_em_and_vortex_throat_interpretation.md

- Compare candidate carriers: circulating intake, brane-tangent vortex flow,
  mixed \(a\)-\(w\) circulation, trapped chiral shear, and an independent
  microrotation of the ordered substructure.
````

````
$ sed -n 1343p docs/native_light_em_and_vortex_throat_interpretation.md
The assumption that the de-structured bulk carries no comparable shear channel helps confine a transverse disturbance, but it is not sufficient by itself. Interface motion, compressional conversion, evanescent tails, and nonlinear radiation can still drain energy. A geon-like state requires a normalizable or sufficiently long-lived bound configuration after all of those channels are included.
````

````
$ sed -n 2377,2381p docs/native_light_em_and_vortex_throat_interpretation.md
### 9.8 A trapped standing wave is not automatically spinning

A real linearly polarized standing wave can have zero time-averaged angular
momentum. Two degenerate modes with a relative phase can form a circularly
polarized bound pattern that carries angular momentum. Orbital winding can
````

## Leg evidence cited in P2/P4 — literal stdout of the round-1 leg scripts (committed at `7b38e9dc`)

````
$ grep -n -E 'RESIDUAL|D_DVZ|D_DUZ' research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/codex/reflection_operator_audit.stdout.txt
2:DOT_VECTOR_RESIDUAL= 0
3:TRACE_GRADIENT_RESIDUAL= 0
4:CONSTRAINT_FOLD_SCALAR_RESIDUAL= 0
5:MATERIAL_ANCHOR_MINUS_U_DOT_G_RESIDUAL= 0
6:TILTED_NORMAL_COVARIANCE_RESIDUAL= Matrix([[0], [0], [0], [0]])
7:FACE_NORMAL_VELOCITY_SCALAR_RESIDUAL= 0
8:Z_INVARIANT_FACE_SCALAR_D_DVZ= 0
9:Z_INVARIANT_ANCHOR_D_DUZ= 0
16:K2_ACTIVE_MIXED_DUZ_DQX= g_y
17:K2_ACTIVE_MIXED_DUZ_DQY= -g_x
28:ODD_SQUARED_REFLECTION_RESIDUAL= 0
````

````
$ grep -n -E 'COEFF_u_z|coeff of d_t u_z|n_plus\[2\]|n_minus\[2\]|has dv_uz|d_z dv_uz' research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/grok/02_face_normal_constraint.out
14:COEFF_u_z_in_n_plus_dot_u  = 0
15:COEFF_u_z_in_n_minus_dot_u = 0
20:V_plus  coeff of d_t u_z = 0
21:V_minus coeff of d_t u_z = 0
29:n_plus[2] exact = 0
30:n_minus[2] exact = 0
34:has dv_uz = False
35:div_P has d_z dv_uz = False
39:COEFF_u_z in material pullback of Q(y) = 0
````

## Round-2 leg evidence cited in the disposition (R2-7, R2-8)

````
$ grep -n -E 'LINEAR_ODD_SCALAR_HESSIAN_AT_BACKGROUND|NONLINEAR_SCALAR_SOURCE' research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/r2_codex/operator_symmetry_audit.stdout.txt
31:NONLINEAR_SCALAR_SOURCE lambda*odd_amp**2
32:LINEAR_ODD_SCALAR_HESSIAN_AT_BACKGROUND 0
````

````
$ cat research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/r2_codex/smooth_profile_fourier_audit.stdout.txt
PROFILE exp(-1/(1-x^2)) for |x|<1, else 0 (C-infinity, nonanalytic)
K 20 ABS_F 0.000561902949952957520530988 NEG_LOG_OVER_K 0.374209070494956553 NEG_LOG_OVER_SQRT_K 1.67351383884746745
K 40 ABS_F 0.000128744611871061807549287 NEG_LOG_OVER_K 0.223941996721032406 NEG_LOG_OVER_SQRT_K 1.41633354680884247
K 80 ABS_F 0.00000839213464027346394662635 NEG_LOG_OVER_K 0.146102695538938832 NEG_LOG_OVER_SQRT_K 1.30678223568409
K 160 ABS_F 1.34109966830722522083066e-8 NEG_LOG_OVER_K 0.113294942615915806 NEG_LOG_OVER_SQRT_K 1.43308026417747616
K 320 ABS_F 4.13712643579755546390872e-10 NEG_LOG_OVER_K 0.0675182796272936264 NEG_LOG_OVER_SQRT_K 1.20780370376374171
````

## v4 author evidence — scripts and literal stdout

````
$ python3 research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v4_author_scripts/symmetry_domain_audit.py
MIRROR_AXIS_1 [[1, 0, 0], [0, 1, 0], [0, 0, -1]]
O2_RANK2_POLAR_FORM [[{'A': 1}, 0, 0], [0, {'B': 1}, 0], [0, 0, {'B': 1}]]
O2_RANK2_POLAR_MIRROR_IMAGE [[{'A': 1}, {}, {}], [{}, {'B': 1}, {}], [{}, {}, {'B': 1}]]
AXIAL_23_SO2_CANDIDATE [[0, 0, 0], [0, 0, {'C': 1}], [0, {'C': -1}, 0]]
AXIAL_23_MIRROR_IMAGE [[{}, {}, {}], [{}, {}, {'C': -1}], [{}, {'C': 1}, {}]]
AXIAL_23_FIXED_BY_MIRROR False
AXIS_1_QUARTER_TURN [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
O2_RANK2_POLAR_QUARTER_TURN_IMAGE [[{'A': 1}, {}, {}], [{}, {'B': 1}, {}], [{}, {}, {'B': 1}]]
PARITY_BLOCK v_face_3 FROM u_3 SAME_SECTOR True
PARITY_BLOCK delta_v_x_3 FROM u_3 SAME_SECTOR True
PARITY_BLOCK V FROM u_3 SAME_SECTOR False
PARITY_BLOCK J FROM u_3 SAME_SECTOR False
PARITY_BLOCK affinity FROM u_3 SAME_SECTOR False
W_COORDINATE_MIRROR_SIGN 1
NORMAL_DERIVATIVE_PRESERVES_INPLANE_PARITY True
ROUND_SECTOR 0 SCALAR_PARITY 1 SPHEROIDAL_PARITY 1 TOROIDAL_PARITY ZERO_SECTOR TOROIDAL_EXISTS False
ROUND_SECTOR 1 SCALAR_PARITY -1 SPHEROIDAL_PARITY -1 TOROIDAL_PARITY 1 TOROIDAL_EXISTS True
ROUND_SECTOR 2 SCALAR_PARITY 1 SPHEROIDAL_PARITY 1 TOROIDAL_PARITY -1 TOROIDAL_EXISTS True
ROUND_SECTOR 3 SCALAR_PARITY -1 SPHEROIDAL_PARITY -1 TOROIDAL_PARITY 1 TOROIDAL_EXISTS True
ROUND_SECTOR 4 SCALAR_PARITY 1 SPHEROIDAL_PARITY 1 TOROIDAL_PARITY -1 TOROIDAL_EXISTS True
CONTROL_ORDER SPECIALIZE_BASELINE_THEN_APPEND_CONTROL
K1_FIXED_DATUM e=(sin(beta),0,cos(beta))
K1_ON_P_U3_THETA_CARRIERS ['a_K1*cos(beta)*d1(u3)*d1(theta)', 'a_K1*cos(beta)*d2(u3)*d2(theta)']
K2_ON_R1_P_CARRIER a_K2*theta*g1*d2(u3)
````

````
$ python3 research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v4_author_scripts/potential_trace_pullback_audit.py
FORMAL_TRACE_COORDINATE 1 delta_p_plus
FORMAL_TRACE_COORDINATE 2 d_w_delta_p_plus
FORMAL_TRACE_COORDINATE 3 delta_v_bulk_plus_1
FORMAL_TRACE_COORDINATE 4 delta_v_bulk_plus_2
FORMAL_TRACE_COORDINATE 5 delta_v_bulk_plus_3
FORMAL_TRACE_COORDINATE 6 delta_v_bulk_plus_4
FORMAL_TRACE_COORDINATE 7 d_w_delta_v_bulk_plus_1
FORMAL_TRACE_COORDINATE 8 d_w_delta_v_bulk_plus_2
FORMAL_TRACE_COORDINATE 9 d_w_delta_v_bulk_plus_3
FORMAL_TRACE_COORDINATE 10 d_w_delta_v_bulk_plus_4
FORMAL_TRACE_COORDINATE 11 delta_p_minus
FORMAL_TRACE_COORDINATE 12 d_w_delta_p_minus
FORMAL_TRACE_COORDINATE 13 delta_v_bulk_minus_1
FORMAL_TRACE_COORDINATE 14 delta_v_bulk_minus_2
FORMAL_TRACE_COORDINATE 15 delta_v_bulk_minus_3
FORMAL_TRACE_COORDINATE 16 delta_v_bulk_minus_4
FORMAL_TRACE_COORDINATE 17 d_w_delta_v_bulk_minus_1
FORMAL_TRACE_COORDINATE 18 d_w_delta_v_bulk_minus_2
FORMAL_TRACE_COORDINATE 19 d_w_delta_v_bulk_minus_3
FORMAL_TRACE_COORDINATE 20 d_w_delta_v_bulk_minus_4
SUPPLIED_CLASS_P_ANSATZ phi=Phi(x1,w)*exp(i*(k2*x2-omega*t))
SUPPLIED_TRACE_PRESSURE delta_p=i*rho_m*omega*Phi
SUPPLIED_TRACE_VELOCITY (d1(Phi), i*k2*Phi, 0, Psi) Psi=d_w(Phi)
SUPPLIED_NORMAL_JET_PRESSURE d_w(delta_p)=i*rho_m*omega*Psi
SUPPLIED_NORMAL_JET_VELOCITY (d1(Psi), i*k2*Psi, 0, d_w(Psi))
SUPPLIED_WAVE_EQUATION_PULLBACK d_w(Psi)=(k2^2-omega^2/c_s0^2)*Phi-d1^2(Phi)
CLASS_P_ODD_BULK_TRACE_COMPONENT 0
CLASS_P_ODD_BULK_NORMAL_JET_COMPONENT 0
VIRTUAL_TEST_COORDINATE 1 delta_v_u_1
VIRTUAL_TEST_COORDINATE 2 delta_v_u_2
VIRTUAL_TEST_COORDINATE 3 delta_v_u_3
VIRTUAL_TEST_COORDINATE 4 delta_v_e_W
VIRTUAL_TEST_COORDINATE 5 delta_v_zeta_c
VIRTUAL_OUTPUT_OBJECT FaceSource.virtual_displacement
D_VIRTUAL_X3_D_TEST_U3 1
D_VIRTUAL_X3_D_PHYSICAL_U3 0
D_FACE_VELOCITY3_D_PHYSICAL_U3 d_t
VIRTUAL_MAP_IS_PHYSICAL_FRECHET_ROW False
INTERFACE_SOURCE_LINE 150 # Traced bulk perturbations and their normal jets at the flat reference faces.
INTERFACE_SOURCE_LINE 154 delta_v_bulk = {
INTERFACE_SOURCE_LINE 158 dw_delta_v_bulk = {
INTERFACE_SOURCE_LINE 163 1: symbol("d_w_delta_p_plus", "COORDINATE", "upper pressure-perturbation normal jet at the flat reference face"),
INTERFACE_SOURCE_LINE 164 -1: symbol("d_w_delta_p_minus", "COORDINATE", "lower pressure-perturbation normal jet at the flat reference face"),
INTERFACE_SOURCE_LINE 638 pressure_perturbation = affine_bulk_perturbation(
INTERFACE_SOURCE_LINE 641 velocity_perturbation = tuple(
GOVERNING_SPEC_LINE 95 Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
````

````
$ python3 research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v4_author_scripts/flow_horizon_source_scope_audit.py
ONTOLOGY_LINE 100 The bulk is therefore an internal reservoir, not an external supply. In a closed-loop realization, localized throat drainage transfers material, momentum, and energy from ordered brane degrees of freedom into de-structured bulk degrees of freedom, while distributed return transfers material back into the ordered state. The bulk may store the corresponding response as compression, flow, internal energy, entropy, or other unresolved excitations of the same medium. Any sustained cycle must include the evolution of that reservoir rather than treating it as an inexhaustible battery.
ONTOLOGY_LINE_HAS_WORD 100 material True
ONTOLOGY_LINE_HAS_WORD 100 energy True
ONTOLOGY_LINE_HAS_WORD 100 particle False
ONTOLOGY_LINE_HAS_WORD 100 coherent False
ONTOLOGY_LINE 1366 - \(\Gamma_{\rm drain}\) converts ordered brane material into de-structured bulk material at throats;
ONTOLOGY_LINE_HAS_WORD 1366 material True
ONTOLOGY_LINE_HAS_WORD 1366 energy False
ONTOLOGY_LINE_HAS_WORD 1366 particle False
ONTOLOGY_LINE_HAS_WORD 1366 coherent False
````

## §P · Prior-art existence checks (WebSearch, 2026-09-30)

Existence and abstract only. This is a verbatim retrieval of what each search returned; ⛔ no full text was read,
and no number or formula quoted by the external search AI was verified.

| query (verbatim) | returned record |
|---|---|
| `Molz Beamish "Leaky plate modes" radiation into a solid medium JASA 1996` | E. B. Molz & J. R. Beamish, "Leaky plate modes: Radiation into a solid medium", JASA 99, 1894–1900 (1996). Abstract: a 60 μm alumina membrane in helium; when the helium freezes, the attenuation of SH(0) and L(1) jumps, attributed to radiation into the solid; SH(0) generates only shear waves in the solid. |
| `Demma Cawley Lowe "Scattering of the fundamental shear horizontal mode from steps and notches in plates" JASA 2003` | A. Demma, P. Cawley, M. Lowe, JASA 113(4), 1880–1891 (2003): FE + modal decomposition, SH0 reflection/transmission at thickness steps and notches. |
| `"Mode conversion of the fundamental shear horizontal wave at a defect" NDT&E International 2025` | C. Peyton, S. Dixon, B. Dutton, W. Vesga, R. S. Edwards, NDT & E Int. 158 (March 2026), doi:10.1016/j.ndteint.2025.103534: SH0 incident on a defect gives a mode-converted S0 reflection depending on defect width/length/depth; FE verified experimentally. |
| `Kubrusly von der Weid Dixon first four SH guided wave modes symmetric non-symmetric discontinuities plates NDT&E 2019` | A. C. Kubrusly, J. P. von der Weid, S. M. Dixon, NDT & E Int. 108 (Dec 2019): symmetric discontinuities create only modes sharing the incident mode's symmetry; asymmetric ones can convert to modes of different symmetry. |
| `Gu Fuller "subsonic wave scattering from discontinuities on fluid-loaded plates" JASA` | Y. Gu & C. R. Fuller, "Active control of sound radiation due to subsonic wave scattering from discontinuities on fluid-loaded plates. I: Far-field pressure", JASA 90, 2020–2026 (1991). |
| `Friedland Giannotti "Astrophysical bounds on photons escaping into extra dimensions" PRL 2008` | A. Friedland & M. Giannotti, PRL 100, 031602 (2008), doi:10.1103/PhysRevLett.100.031602: warped single-brane photon localization; photon metastability in plasma; stellar-cooling bounds. |
| (title appeared in the Demma search listing; content not retrieved) | "The scattering of the fundamental torsional mode from axi-symmetric defects with varying depth profile in pipes", JASA 127(6), 3440. |
