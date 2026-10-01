# S11c_d_clean_condition.md (v5) — grounding commands

Mechanical lookups, rerun from the repo root on 2026-10-01 against the v4 baseline `5693e861` plus the v5
working-tree edits. Each command and its literal output are shown. Fences are longer than fences in quoted text.

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
$ sed -n 5p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
On the user's request to assess before continuing and stop for go/no-go, the existing central-balance-v2 scientific child was suspended in memory, not terminated or restarted. The pause receipt records PID 4097233, its process start identity, cgroup and exact command. It had completed two matrix groups after twenty restorations and saved 12,288 batches of the next group; no finite solve had been reached. The guard/supervisor remain active, no deadline was added, and pinned sources are unchanged. Resume requires the user's decision and rechecking that process identity. Do not launch another job alongside it.
````

````
$ sed -n 11,13p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
The development input sets mu_R=1, rho_br=1 and c_s0=10 in the L_ref/T_ref unit frame. Under the bare S9/R4 definition c_gamma^2=mu_R/rho_br, this gives c_gamma/c_s=0.1, not 1. The step plan explicitly distinguishes the derived ratio definition from the calibrated, uncommitted equality lambda_gamma=1; the equality has not been imposed on these inputs.

There is an additional distinction: that bare coefficient formula must not be substituted for the actual full S11c-d transverse dispersion. The saved omega=3 LEFT incoming branch has normal momentum 2.439262183530097 and tangential norm squared 1/20. Its fixed-point phase speed omega/sqrt(k_normal^2+1/20) is 1.2247448713915874, or 0.12247448713915873 of c_s0. At the full-contrast RIGHT end the corresponding ratio is 0.12186666955535794; its zero-contrast ratio returns to the LEFT value. These are arithmetic readouts of saved full-pencil modes, not a new root solve or a demonstrated identification with the calibrated light cone. The extra retained elastic structure and coefficient conventions cannot be silently discarded.
````

````
$ sed -n 47p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_equal_speed_feasibility.md
The same PID4097233 was resumed with SIGCONT after checking start ticks85113744, exact command/cgroup, 2GiB/zero-swap/task containment, systemd infinity/Restart=no and the existing waiting completion hook. No source/gate/helper was edited and no new scientific process was launched. See `S11c_d_numerical_radiating_flow_calibration_resume.json`.
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
$ sed -n 90,97p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
`v_bulk_normal_0` (the bulk normal drain, `S11b_SHARED_PHYSICS.md:104`) is a scope-limit parameter, not an
active DOF, and appears in no derived operator (§0).

### 1b · Bulk acoustics and the projection law

Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
`∂_t²φ=c_s0²∇₄²φ`; the current and conservation law `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; and the
dynamic, anchored slab window `Ω` supplied in S11c-a §3. S11c-b performs no curved-bulk response solve (§0).
````

````
$ sed -n 115,121p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
The uniform-basis inputs are the fields and their first gradients `{u,∇u,θ,∇θ,e_W,∇e_W}`, with the S11b
symmetry group in full (`S11b_SHARED_PHYSICS.md:280–288`): **in-plane translation invariance** (so `u` enters
only through its gradients, never undifferentiated), **in-plane `O(3)` isotropy and parity**, **reflection
`w→−w`**, **equivalence modulo total in-plane divergences**, and **no time-reversal, positivity, or
boundedness** assumption. Independence of invariants is judged as field bilinears with B1's constraint **not**
applied; every independent invariant is carried with its own free symbolic coefficient. These are the same
construction rules S11b task B0-energy used; §3a inherits them with variable coefficients.
````

````
$ sed -n 145,154p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
δ_vΣ_mat = 0 ,   (uniform linearisation)  δ_vθ + δ_ve_W + ∇_x·δ_vu = 0 .
```

The three `τ_I` are independent. Equations of motion are obtained by S11b's method — balance laws, the
binding virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed
external virtual work — **not** by putting an irreversible response kernel in an ordinary action.
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
$ sed -n 112,120p docs/native_light_em_and_vortex_throat_interpretation.md
### 9. Spin, chirality, helicity, and magnetic moment — *drafted*

- Compare candidate carriers: circulating intake, brane-tangent vortex flow,
  mixed \(a\)-\(w\) circulation, trapped chiral shear, and an independent
  microrotation of the ordered substructure.
- Begin with the angular-momentum tensor \(J_{IJ}\), not a borrowed
  three-dimensional spin vector.
- Separate circulation sign, throat orientation, chirality, helicity,
  intrinsic angular momentum, and magnetic moment.
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

## v5 author computations — commands and literal stdout

````
$ python3 research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v5_author_scripts/bulk_pullback_pairing_audit.py
PRESSURE_TRACE I*omega*rho_m*Phi(x1, w)
VELOCITY_TRACE [Derivative(Phi(x1, w), x1), I*k2*Phi(x1, w), 0, Derivative(Phi(x1, w), w)]
PRESSURE_NORMAL_JET I*omega*rho_m*Derivative(Phi(x1, w), w)
VELOCITY_NORMAL_JET [Derivative(Phi(x1, w), w, x1), I*k2*Derivative(Phi(x1, w), w), 0, k2**2*Phi(x1, w) - Derivative(Phi(x1, w), (x1, 2)) - omega**2*Phi(x1, w)/c_s0**2]
DENSITY_TRACE I*omega*rho_m*Phi(x1, w)/c_s0**2
DENSITY_MINUS_PRESSURE_OVER_C2 0
DENSITY_NORMAL_JET I*omega*rho_m*Derivative(Phi(x1, w), w)/c_s0**2
DENSITY_TIME omega**2*rho_m*Phi(x1, w)/c_s0**2
CURRENT_TRACE [rho_m*Derivative(Phi(x1, w), x1), I*k2*rho_m*Phi(x1, w), 0, rho_m*Derivative(Phi(x1, w), w)]
CURRENT_NORMAL_JET [rho_m*Derivative(Phi(x1, w), w, x1), I*k2*rho_m*Derivative(Phi(x1, w), w), 0, k2**2*rho_m*Phi(x1, w) - rho_m*Derivative(Phi(x1, w), (x1, 2)) - omega**2*rho_m*Phi(x1, w)/c_s0**2]
TRACE_EXACT_MINUS_REFERENCE_AFFINE 0
TRACE_EXACT_MINUS_BACKGROUND_AFFINE -W0*eta*s*w1(x1)*Derivative(Phi(x1, w), w)/2
PAIRING_SAME_PHASE 0
PAIRING_CONJUGATE_PHASE 2*I*pi*test_amplitude*trial_amplitude
````

````
$ python3 research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v5_author_scripts/round_spin_carrier_audit.py
HELICITY_TOROIDAL 0
HELICITY_TOROIDAL_QUADRATURE 0
HELICITY_TOROIDAL_POLOIDAL 5*sqrt(2)*pi**(3/2)*mixing/4
PARITY_TOROIDAL 1
PARITY_POLOIDAL -1
ANGULAR_MOMENTUM_TOROIDAL_QUADRATURE [0, 0, sqrt(2)*pi**(3/2)*omega/16]
CORIOLIS_SCALAR_L0_PROJECTION -sqrt(2)*pi**(3/2)*Omega/4
CORIOLIS_SCALAR_L2_PROJECTION -sqrt(2)*pi**(3/2)*Omega/8
O3_FIXED_VECTOR_SOLUTIONS [{Lx: 0, Ly: 0, Lz: 0}]
````

## Earlier computed term-structure evidence retained by P4

````
$ grep -n -E 'CONSTRAINT_FOLD_SCALAR_RESIDUAL|MATERIAL_ANCHOR_MINUS_U_DOT_G_RESIDUAL' research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/codex/reflection_operator_audit.stdout.txt
4:CONSTRAINT_FOLD_SCALAR_RESIDUAL= 0
5:MATERIAL_ANCHOR_MINUS_U_DOT_G_RESIDUAL= 0
````

````
$ grep -n -E 'COEFF_u_z_in_n_(plus|minus)_dot_u|COEFF_u_z in material pullback' research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/grok/02_face_normal_constraint.out
14:COEFF_u_z_in_n_plus_dot_u  = 0
15:COEFF_u_z_in_n_minus_dot_u = 0
39:COEFF_u_z in material pullback of Q(y) = 0
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
