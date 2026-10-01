# S11c_d_clean_condition.md (v2) — grounding commands (rule 2 / E1)

Mechanical lookups, run from the repo root on 2026-09-30, HEAD `7b38e9dc`. Regenerated from the commands below; nothing transcribed.

```
$ sed -n 3p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
The current pilot answers a strict rest-bulk, LAB_HELD/RHO4_CONSTANT development-input question. It does not yet answer leakage in the calibrated, draining medium. These limitations were explicit in the governing scope, but should have been foregrounded in the pilot decision. No result can remove them through numerical convergence alone.
```

```
$ sed -n 11,15p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
The development input sets mu_R=1, rho_br=1 and c_s0=10 in the L_ref/T_ref unit frame. Under the bare S9/R4 definition c_gamma^2=mu_R/rho_br, this gives c_gamma/c_s=0.1, not 1. The step plan explicitly distinguishes the derived ratio definition from the calibrated, uncommitted equality lambda_gamma=1; the equality has not been imposed on these inputs.

There is an additional distinction: that bare coefficient formula must not be substituted for the actual full S11c-d transverse dispersion. The saved omega=3 LEFT incoming branch has normal momentum 2.439262183530097 and tangential norm squared 1/20. Its fixed-point phase speed omega/sqrt(k_normal^2+1/20) is 1.2247448713915874, or 0.12247448713915873 of c_s0. At the full-contrast RIGHT end the corresponding ratio is 0.12186666955535794; its zero-contrast ratio returns to the LEFT value. These are arithmetic readouts of saved full-pencil modes, not a new root solve or a demonstrated identification with the calibrated light cone. The extra retained elastic structure and coefficient conventions cannot be silently discarded.

Thus neither the bare input ratio nor the actual selected modal ratio realizes equal light/bulk speeds. Identifying this development slice with the calibrated analog-light band remains OPEN. Changing that calibration would change physical inputs and the branch/phase-matching problem, not just relabel a plot.
```

```
$ sed -n 5p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md
On the user's request to assess before continuing and stop for go/no-go, the existing central-balance-v2 scientific child was suspended in memory, not terminated or restarted. The pause receipt records PID 4097233, its process start identity, cgroup and exact command. It had completed two matrix groups after twenty restorations and saved 12,288 batches of the next group; no finite solve had been reached. The guard/supervisor remain active, no deadline was added, and pinned sources are unchanged. Resume requires the user's decision and rechecking that process identity. Do not launch another job alongside it.
```

```
$ sed -n 1104,1111p research/pde_ledger_v3/V3_STEP_PLAN.md
`CONE_LOCK_CALIBRATED` with `λγ = 1` as a **calibration, not a derivation**; and `R3` already encodes the
ratio in the registry seed.

⛔ **Two different objects — round 2. Do not merge them:** the **ratio** `λγ = c_γ/c_s` is a *derived
definition* (registry `R3`); the **equality `λγ = 1`** is a separate *calibrated / uncommitted* cone
lock. Classifying the ratio says nothing about the lock.

⇒ ⛔ **Both classifications are fixed, not open:** the **ratio** `λγ = c_γ/c_s` is **`derived`** (it is `R3`); the **equality `λγ = 1`** is **`calibrated`/uncommitted**. ⛔ Do not re-open the ratio's class. Introduce `λγ` with provenance, confront
```

```
$ sed -n 1116,1126p research/pde_ledger_v3/V3_STEP_PLAN.md
⭐⭐ **`λγ = 1` is not a free landing — observation already constrains it.** In this model `c_s` is the
**gravity-change/phonon speed** and `c_γ` is the **light-cone speed**
(`research/pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md:75-77`):

> *"`c_s` is the phonon/gravity-change speed; `c_gamma` is the light-cone speed; `v_b` is the condensate
> flow."*

Gravitational-wave and electromagnetic arrival-time observations (**GW170817 / GRB 170817A**) constrain
those two speeds to agree to roughly **1 part in 10¹⁵**. ⚠ **The observational value IS known** — `λγ = 1`
to that precision is not in doubt as a *target*.

```

```
$ sed -n 35p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
density representatives ρ4D/ρbr, two anchorings LAB_HELD/MATERIAL_ADVECTED), S11c-b computes: (1) the **§3a energy
```

```
$ sed -n 90,91p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
`v_bulk_normal_0` (the bulk normal drain, `S11b_SHARED_PHYSICS.md:104`) is a scope-limit parameter, not an
active DOF, and appears in no derived operator (§0).
```

```
$ sed -n 95,97p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
`∂_t²φ=c_s0²∇₄²φ`; the current and conservation law `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; and the
dynamic, anchored slab window `Ω` supplied in S11c-a §3. S11c-b performs no curved-bulk response solve (§0).
```

```
$ sed -n 145,148p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
```

```
$ sed -n 350,354p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
```text
J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α ,
n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m ,
𝒜_s^α ≡ μ_s^α − δp_s^α/ρ_m ,
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α)n̂_s^α ,
```

```
$ sed -n 26p docs/toy_model_ontology_summary.md
The model is a classical \(4+1\)-dimensional, one-medium analog in which our three-dimensional space is an ordered, shear-supporting brane; particles are finite, oriented throats through it; a trapped transverse brane-shear standing mode helps hold each throat open; throat transport, order conversion, stress, return, and reservoir response are the proposed material mechanism underlying a calibrated gravity observable; and freely propagating light, electric coupling, magnetism, radiation, passive response, inertia, and possible cosmic expansion are target behaviors of the same complete material configuration.
```

```
$ sed -n 100p docs/toy_model_ontology_summary.md
The bulk is therefore an internal reservoir, not an external supply. In a closed-loop realization, localized throat drainage transfers material, momentum, and energy from ordered brane degrees of freedom into de-structured bulk degrees of freedom, while distributed return transfers material back into the ordered state. The bulk may store the corresponding response as compression, flow, internal energy, entropy, or other unresolved excitations of the same medium. Any sustained cycle must include the evolution of that reservoir rather than treating it as an inexhaustible battery.
```

```
$ sed -n 315p docs/toy_model_ontology_summary.md
The fields \(h_+\) and \(h_-\) describe the interfaces only while each remains a single-valued graph over the brane coordinates. They are useful background, far-field, and separated-interface collective coordinates, not a complete nonlinear throat topology. Overhangs, necks, tubes, folded sheets, reconnection, pinch-off, and topology change must be represented by the parent fields \(\chi_B(\mathbf x,w),n(\mathbf x,w),\theta(\mathbf x,w)\), the support fields, and core level sets. The earlier single coordinate \(h\) should be read only as schematic notation for a derived normal eigencombination.
```

```
$ sed -n 362p docs/toy_model_ontology_summary.md
- **Structural support:** energy \(E_{\rm support}\) in a trapped transverse brane-shear standing mode—the model's light or photon mode—helps hold the aperture open.
```

```
$ sed -n 957p docs/toy_model_ontology_summary.md
Freely propagating light and the throat-support mode must be distinguished. Free light is a background-brane \(\mathbf u_T\) excitation. The support mode is a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator; it is not a second photon substance and is not assumed to propagate as a shear wave through fully de-structured bulk material. Any \(\Omega_{\rm support}\) is a diagnostic energy/stress-localization region extracted from that solved mode, not a hard PDE domain imposed by \(\mu_R>0\).
```

```
$ sed -n 1366p docs/toy_model_ontology_summary.md
- \(\Gamma_{\rm drain}\) converts ordered brane material into de-structured bulk material at throats;
```

```
$ sed -n 113,116p docs/native_light_em_and_vortex_throat_interpretation.md

- Compare candidate carriers: circulating intake, brane-tangent vortex flow,
  mixed \(a\)-\(w\) circulation, trapped chiral shear, and an independent
  microrotation of the ordered substructure.
```

```
$ sed -n 1343p docs/native_light_em_and_vortex_throat_interpretation.md
The assumption that the de-structured bulk carries no comparable shear channel helps confine a transverse disturbance, but it is not sufficient by itself. Interface motion, compressional conversion, evanescent tails, and nonlinear radiation can still drain energy. A geon-like state requires a normalizable or sufficiently long-lived bound configuration after all of those channels are included.
```

```
$ sed -n 2377,2380p docs/native_light_em_and_vortex_throat_interpretation.md
### 9.8 A trapped standing wave is not automatically spinning

A real linearly polarized standing wave can have zero time-averaged angular
momentum. Two degenerate modes with a relative phase can form a circularly
```

```
$ ls research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/codex research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/grok
research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/codex:
actual_source_structure_audit.sh
actual_source_structure_audit.stdout.txt
citation_grounding_audit.py
citation_grounding_audit.stdout.txt
reflection_operator_audit.py
reflection_operator_audit.stdout.txt
rotation_and_spin_audit.py
rotation_and_spin_audit.stdout.txt
round_toroidal_audit.py
round_toroidal_audit.stdout.txt
speed_gap_scattering_audit.py
speed_gap_scattering_audit.stdout.txt

research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_scripts/grok:
01_parity_and_invariants.out
01_parity_and_invariants.py
02_face_normal_constraint.out
02_face_normal_constraint.py
03_k1_k2_form_controls.out
03_k1_k2_form_controls.py
04_round_parity_and_rotation.out
04_round_parity_and_rotation.py
05_citation_verify.out
05_citation_verify.py
06_levi_civita_presence.out
```

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
