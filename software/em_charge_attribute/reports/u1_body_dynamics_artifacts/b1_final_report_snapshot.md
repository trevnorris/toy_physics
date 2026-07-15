# U1 throat-body dynamics — Decision-16 OPERATIVE Phase-A baseline

Input SHA-256: `eef5cdc2178ec76f10f3e5e361cdd1461eb2eee3980ffa99f814df02475866ee`. Phase A Axis 1 is **`COMPUTATION_VALID`** and Axis 2 is **`U1_BASE_OK`**.
Baseline status is **`OPERATIVE`** under `decision_16_retire_brane_polar_field`. This report halts after spec 7.0/7.1; Phase B and Phase C are `NOT_RUN(upstream)` by process.

## Frozen setup and honest scope

The medium-rest lab frame, co-moving steady family, co-moving `Omega_c`, ambient-subtracted exterior-ball IR scheme, fixed `mdot`, and future `C_mdot` ownership of acceleration-like outer flux were declared before residual evaluation. The family is an action-derived exterior solution joined at `r=a` to field-level core traces through the typed throat surface functional. This is a collective-coordinate/effective-action computation, not a nonlinear throat simulation.

No production input contains a tail exponent, eigenvalue sign, class boolean, or verdict. `Mh` and `cE` remain symbolic positive OPEN action coefficients; no collective-coordinate or multipole functional is an OPEN leaf. The independent brane polar field, its traces, and its action knobs are absent from the operative input.

Retained finding (historical, not part of the operative action): commit `a5c079eb` established the long-wavelength helical instability of the now-retired massless-P plus `lambdaPu` baseline. Decision 16 retains that result as evidence for retirement.

## 7.0 — declared-input ledger

| root | status | type | domain | [L,T,M] | arguments | symmetry | source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sleeve_core_trace | OPEN_INPUT | PRIMITIVE_OPEN | scalar_field_trace_on_Sigma | (0, 0, 0) | field_point, outward_normal, wall_side, a, L | BODY_CONJUGATE_SCALAR | nonlinear_throat_core_open |
| geon_core_bundle | OPEN_INPUT | PRIMITIVE_OPEN | field_trace_bundle_on_Sigma | (0, 0, 0) | density_trace, h_trace, phase_flux_trace, wall_side | BODY_CONJUGATE_BUNDLE | spec_7_0_geon_profile_open |
| throat_surface_functional | OPEN_INPUT | ACTION | field_functional_on_Sigma | (2, -2, 1) | field_traces, normal_derivatives, induced_metric | BODY_CONJUGATE_SCALAR | ledger_stage006_f_throat |
| outer_surface_functional | OPEN_INPUT | ACTION | field_functional_on_partial_Omega_c | (2, -2, 1) | field_traces, outward_normal, ambient_branch | AMBIENT_BRANCH_COVARIANT | return_surface_action_open |
| native_continuity | WHITELIST_SOURCED | BALANCE | fields_on_time_x_Omega_c | (-4, -1, 0) | number_density_field, velocity_field | R_W_COVARIANT | ledger_stage004_stage006 |
| native_momentum | WHITELIST_SOURCED | BALANCE | fields_on_time_x_Omega_c | (-3, -2, 1) | momentum_density_field, native_stress_tensor, action_force_density | R_W_COVARIANT | GNLS_wall_shear_h_action |
| E4_shear_lock | OPEN_INPUT | CONSTRAINT | velocity_trace_functional_on_collar | (1, -1, 0) | bulk_velocity_trace, shear_velocity_trace, collar_kernel | BODY_CONJUGATE_VECTOR | U2_seam |
| E5_rayleigh | OPEN_INPUT | RAYLEIGH | tangential_velocity_trace_functional_on_Sigma | (2, -3, 1) | bulk_velocity_trace, sleeve_velocity_trace, gammaSigma | BODY_CONJUGATE_SCALAR | U2_seam |
| return_closure | OPEN_INPUT | RETURN | structured_field_flux_map_on_partial_Omega_c | (0, -1, 1) | mass_flux_density, momentum_flux_density, destination_field | AMBIENT_BRANCH_COVARIANT | pathA_29_open_return |
| hbar | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (2, -1, 1) | none | R_W_EVEN | ledger_stage004 |
| m_GNLS | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (0, 0, 1) | none | R_W_EVEN | ledger_stage004 |
| rho_inf | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-4, 0, 0) | none | R_W_EVEN | ledger_stage004 |
| K_EOS | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (18, -2, 1) | none | R_W_EVEN | ledger_stage004_stage005 |
| aB | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-2, -2, 1) | none | R_W_EVEN | ledger_stage006 |
| kappaB | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (0, -2, 1) | none | R_W_EVEN | ledger_stage006 |
| muR4 | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-2, -2, 1) | none | R_W_EVEN | ledger_stage006 |
| rhoBr | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-3, 0, 1) | none | R_W_EVEN | ledger_stage007 |
| muR | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-1, -2, 1) | none | R_W_EVEN | ledger_stage007 |
| OmegaW | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (0, -1, 0) | none | R_W_EVEN | ledger_stage007 |
| ellg | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (1, 0, 0) | none | R_W_EVEN | ledger_stage007 |
| Mh | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, -2, 1) | none | R_W_EVEN | pathA_38_dynamic_normalization_open |
| cE | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (1, -1, 0) | none | R_W_EVEN | pathA_38_dynamic_speed_open |
| mdot | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, -1, 1) | none | R_W_EVEN | spec_7_0_fixed_intake |
| gammaSigma | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, 0, 0) | none | R_W_EVEN | U2_E5_field_Rayleigh |
| tangentDtN | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, 0, 0) | none | R_W_EVEN | U2_E5_field_surface_DtN |

### Action source completeness

Both engines loaded `software/stage1_solver/decisions/16_retire_brane_polar_field.md`, verified its SHA-256 `7683c9e8b747b49b80974b58967e455581f02035587db52ed6bf8673c9663e87`, parsed all three required Decision-16 fragments, and confirmed status `OPERATIVE`. They also loaded every cited stage-note file and parsed every executable action expression.

The exact P-retired cover has 15 terms: GNLS Berry/flow/EOS, Madelung quantum pressure, wall well/gradient/shear and typed throat/mix functionals, MacCullagh shear plus `u_w`, and the OPEN `h` kinetic/gradient normalization. The legacy three-term T0 polar block and `L_Pu` are parsed only as retired evidence; they are not assembled.

| assembled term | loaded source | parsed source fragment |
| --- | --- | --- |
| bulk_berry | research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md | `i hbar psi* d_t psi` |
| bulk_flow_kinetic | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `(1/2) m_GNLS n \|u\|^2` |
| quantum_pressure | software/stage1_solver/reports/pathA_35_G0_freeze.md | `QP := (hbar^2/(8 m rho))` |
| bulk_EOS | research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md | `U(rho) = (K/4) rho^5` |
| wall_double_well | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_B(chi_B) = a_B chi_B^2 (1-chi_B)^2` |
| wall_gradient | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `(kappa_B/2)\|grad_4 chi_B\|^2` |
| wall_shear_gate | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `chi_B f_shear` |
| throat_source | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_throat` |
| wall_mix | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_mix` |
| brane_shear_kinetic | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br (partial_t u^a)(partial_t u^a)` |
| brane_shear_gradient | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) mu_R Omega_u^a Omega_u^a` |
| uw_kinetic | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br (partial_t u_w)^2` |
| uw_gap | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br Omega_w^2 u_w^2` |
| h_kinetic | software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md | `dynamic` |
| h_gradient | software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md | `radial branch` |

SOURCE_COMPLETENESS compares the assembled IDs to the exact Decision-16 P-retired cover. It fails on any missing live term, any unexpected term, and any P-sector expression. A reintroduced P term without an explicit Decision-16 citation fails at the citation check; even a cited P term remains forbidden from the operative assembly. The production removal probes independently show that quantum pressure and the MacCullagh shear gradient reach their derived operator entries.

### Inline dimensions, IR, and G9 declaration

| constructed expression | computed [L,T,M] |
| --- | --- |
| bulk_berry | (-2, -2, 1) |
| bulk_flow_kinetic | (-2, -2, 1) |
| quantum_pressure | (-2, -2, 1) |
| bulk_EOS | (-2, -2, 1) |
| wall_double_well | (-2, -2, 1) |
| wall_gradient | (-2, -2, 1) |
| wall_shear_gate | (-2, -2, 1) |
| brane_shear_kinetic | (-2, -2, 1) |
| brane_shear_gradient | (-2, -2, 1) |
| uw_kinetic | (-2, -2, 1) |
| uw_gap | (-2, -2, 1) |
| h_kinetic | (-2, -2, 1) |
| h_gradient | (-2, -2, 1) |
| momentum_rate_density | (-3, -2, 1) |
| control_surface_force | (1, -2, 1) |

Predeclared G9 force scale: `mdot*sqrt(5*K_EOS*rho_inf^4/m_GNLS)`; relative bound: `(a/R)^2+epsilon_rigid^2+epsilon_quad` with `epsilon_rigid=1/100` and `epsilon_quad=1/10000`. It was not fitted to a measured residual.

### Computed co-moving laws

```text
partial_t n+div_4(n v)=0
partial_t|y n+div_4[n(v-V)]=0
partial_t g_i+partial_j Pi_ij=f_i[action]
partial_t|y g_i+partial_j(Pi_ij-V_j g_i)=f_i[action]
integral_Omega f_i d4y-integral_partialOmega(Pi_ij-V_j g_i)N_j dSigma3
```

The composite field `f(-V*t + x)` is differentiated before either balance law is reduced. The scalar chain-rule, continuity, momentum, and declared-flux residuals are all zero. The surface force is a native control-volume balance; no particle `F=dP/dt` form is imported.

## 7.1 — far-field solve that decides G1

For a bulk channel of radial dimension `d` and harmonic number `ell`, both engines construct and solve

```text
A_c [f''+(d-1)f'/r-ell(ell+d-2)f/r^2] - B_c f = 0.
```

`B_density=e''(rho_inf)` is differentiated in-engine from `K rho^5/4`, giving `5 K rho_inf^3`; `B_chi` is likewise differentiated from the wall potential. The translated norm uses the full measure: `S_(d-1) integral r^(d-1) |f'|^2 dr`. For the brane shear channel the four-dimensional measure factorizes as `d^3x g_ell(w)^2 dw`, so its reported radial dimension is three while the `w` integral is finite.

| channel | d | solved radial equation | tail class | nu or gap | translated norm | value at R=a=1 | normalizable |
| --- | --- | --- | --- | --- | --- | --- | --- |
| density_EOS | 4 | Eq(-5*f(r) + Derivative(f(r), (r, 2))/4 + 3*Derivative(f(r), r)/(4*r), 0) | EXPONENTIAL_GAP | 2*sqrt(5) | Integral[1,inf] r^3 (d/dr(r^-1.0 K_1(gap*r)))^2 dr | 0.000204528635461 | True |
| wall_chiB | 4 | Eq(-3*f(r) + 2*Derivative(f(r), (r, 2)) + 6*Derivative(f(r), r)/r, 0) | EXPONENTIAL_GAP | sqrt(6)/2 | Integral[1,inf] r^3 (d/dr(r^-1.0 K_1(gap*r)))^2 dr | 0.409684761529 | True |
| bound_phase | 4 | Eq(Derivative(f(r), (r, 2)) + 3*Derivative(f(r), r)/r, 0) | POWER_LAW | 2 | 2/R**2 | 2 | True |
| brane_shear | 3 | Eq(5*Derivative(f(r), (r, 2))/4 + 5*Derivative(f(r), r)/(2*r) - 5*f(r)/(2*r**2), 0) | POWER_LAW | 2 | 4/(3*R**3) | 1.33333333333 | True |
| uw | 3 | 196/75*f(r)=0 | ALGEBRAIC_GAP | algebraic | 0 | 0 | True |
| h | 4 | Eq(Mh*(Derivative(f(r), (r, 2)) + 3*Derivative(f(r), r)/r), 0) | POWER_LAW | 2 | 2/R**2 | 2 | True |

### Surviving coupled operator result

The drain background produces a density–phase mixed Hessian. Its radial degree differs from the phase diagonal by `-2` in `d=4`; both engines compute `SUBLEADING_FOR_D_GT_2`, so it does not change the scalar tail classes or the six-channel verdict. There is no P–u block in the operative operator.

The bound-flow phase is selected from the solved continuity family by the `mdot` surface normalization; its computed normalization residual is `0`. No exponent is supplied in YAML.

### Force-balance operator and quotient

```text
D_E=(-A_n Laplacian+U''(n0)) on delta_n, plus translated throat-source block; analogous Hessian/balance blocks for chiB,u,h and endpoint traces
Z_A=-partial_A Phi_0 including the translated throat-source block
D_E Z residual = 0
Gram = z1**2 + z2**2
Q Z = ['0', '0']
```

The base source is constructed from the declared profile and the native force-balance equation; changing the field after constructing the source makes `BASE_BALANCE` fail. Translation is then checked by Cartesian substitution into `D_E`, not by a generator-scale identity.

### Endpoint response solves

| endpoint | field condition | class | normal response | tangent response | shear response |
| --- | --- | --- | --- | --- | --- |
| E1 | v.normal=V.normal and v.tangent=V.tangent on Sigma | holonomic_field_BC | 1 | 1 | 0 |
| E2 | v.normal=V.normal and tangential_traction=0 on Sigma | holonomic_field_BC | 1 | 0 | 0 |
| E3 | permeable phase texture; no velocity trace tied to V | bulk_action | 0 | 0 | 0 |
| E4 | g_A=V_A-C_A[uT_dot at collar]=0 | nonholonomic_constraint | 0 | 0 | 1 |
| E5 | v.normal=V.normal; tangentDtN*v.tangent+gammaSigma*(v.tangent-V.tangent)=0 | Rayleigh | 1 | 4/9 | 0 |

E1–E5 are five distinct exact rank-3 boundary/constraint solves rather than identity systems: no-slip, free-slip, permeable texture, nonholonomic shear-lock, and Robin/Rayleigh respectively. E4 additionally constructs `delta W=lambda_A(delta V_A-C_A[delta udot_T])`, whose allowed-variation residual is zero.

### Evaluated reduced moments and L_eff

All angular factors and radial integrals are executed. The following production values show the computed endpoint dependence; the machine result retains every symbolic field-trace expression and dependency set.

| moment | E1 | E2 | E3 | E4 | E5 |
| --- | --- | --- | --- | --- | --- |
| B_X | 0 | 0 | 0 | 0 | 0 |
| B_Xp | 0 | 0 | 0 | 0 | 0 |
| B_p | 0 | 0 | 0 | 0 | 0 |
| H_PP | 0.068538919452 | 0.068538919452 | 0.068538919452 | 0.068538919452 | 0.068538919452 |
| H_XP | 0 | 0 | 0 | 0 | 0 |
| H_XX | 0.201420497981 | 0.201420497981 | 0.201420497981 | 0.201420497981 | 0.201420497981 |
| I_density_grad | 0.756796365754 | 0.756796365754 | 0.756796365754 | 0.756796365754 | 0.756796365754 |
| I_density_l2 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 |
| I_h_grad | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 |
| I_shear_grad | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 |
| I_wall_grad | 1.28479139854 | 1.28479139854 | 1.28479139854 | 1.28479139854 | 1.28479139854 |
| I_wall_l2 | 0.204279992179 | 0.204279992179 | 0.204279992179 | 0.204279992179 | 0.204279992179 |
| N_UU | 10.3295924476 | 2.5823981119 | 0 | 0 | 4.11270810413 |
| N_UW | 1.14773249418 | 0.286933123544 | 0 | 0 | 0.669510621603 |
| N_WW | 0.127525832686 | 0.127525832686 | 0.127525832686 | 0.127525832686 | 0.127525832686 |
| U_PP | 0.034618100866 | 0.034618100866 | 0.034618100866 | 0.034618100866 | 0.034618100866 |
| U_XP | 0.0634665182543 | 0.0634665182543 | 0.0634665182543 | 0.44426562778 | 0.0634665182543 |
| U_XX | 0.209439510239 | 0.209439510239 | 0.209439510239 | 10.2625360017 | 0.209439510239 |

For each endpoint the directly embedded action reduces to

```text
L_eff = A_X V + A_p pdot + C_Xp p V + G_VV V^2/2 + G_Vp V pdot + G_pp pdot^2/2 - K_pp p^2/2.
```

- `E1`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 2.2999402325143*density_delta*m_GNLS + 9.8696044010894*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=2.2999402325143*density_delta*m_GNLS*tilt_phase + 9.8696044010894*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 7.5398223686155*muR*tilt_shear**2`.
- `E2`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 0.57498505812856*density_delta*m_GNLS + 2.4674011002723*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=0.57498505812856*density_delta*m_GNLS*tilt_phase + 2.4674011002723*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 7.5398223686155*muR*tilt_shear**2`.
- `E3`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 7.5398223686155*muR*tilt_shear**2`.
- `E4`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 7.5398223686155*rhoBr*shear_transverse**2 + 15.079644737231*rhoBr*shear_transverse + 7.5398223686155*rhoBr`; `G_Vp=4.1887902047864*rhoBr*shear_transverse*tilt_shear + 4.1887902047864*rhoBr*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 7.5398223686155*muR*tilt_shear**2`.
- `E5`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 0.91571694442697*density_delta*m_GNLS + 3.9295647152485*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=1.3416318023*density_delta*m_GNLS*tilt_phase + 5.7572692339688*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 7.5398223686155*muR*tilt_shear**2`.



<!-- PHASE_A_MOMENT_CORRECTION_NOTE_START -->
**Phase-A amendment note (surfaced by Phase B1):** the `G_VV` rows immediately above display the frozen pre-amendment brane-shear unit coefficient `12π/5` (`7.539822368616`). `PHASE_A_MOMENT_CORRECTION(brane_shear)` corrects the governed `U_XX → G_VV` chain to `8π/3` (`8.377580409573`). The separately undetermined `U_XP`, `U_PP`, and `I_shear_grad` tilt-profile rows remain frozen/`UNRESOLVED(tilt_profile)`. Payload digest: `a32c25f4325671d280b54df6c51abd9b25008ef5e6008b98972bac1ed81f7e69` → `b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d`.
<!-- PHASE_A_MOMENT_CORRECTION_NOTE_END -->
Direct differentiation computes both canonical momenta and `Q_p`. Reconstruction independently iterates over the parsed action expressions, substitutes the rigid field embedding, performs the moment reductions, and compares that sum with the claimed coefficient form; the residual is zero for E1–E5. `K_pp` contains the surviving action-derived density, wall, shear, and `h` stiffness terms and no retired P contribution. E4 multiplier reaction and E5 Rayleigh loss remain outside `L_eff` in their typed channels.

## Gates and able-to-fail evidence

| gate | production result |
| --- | --- |
| G1 | PASS |
| G10 | NOT_RUN(phase_C) |
| G11 | NOT_RUN(phase_C) |
| G2 | PASS |
| G3 | PASS |
| G4 | NOT_RUN(phase_B) |
| G5 | PASS |
| G6 | PASS_ENDPOINT_MAP_REPORTED |
| G7 | NOT_RUN(phase_B) |
| G8 | NOT_RUN(phase_C) |
| G9 | NOT_RUN(phase_B;tolerance_predeclared) |

G1 passes because all six translated tails are normalizable, while G2 passes because the displayed norms and reduced coefficients are finite in the predeclared ambient-subtracted exterior scheme. G3 is `PASS` from the computed Axis-2 verdict. G5 passes with body-only, symmetric-postulate, and one-sided-asymmetry tags kept distinct.

### Per-tooth object ablations

| tooth | SymPy | Mathematica |
| --- | --- | --- |
| INPUT_LEDGER | EXIT_1@ASSERT_FAIL:INPUT_LEDGER | EXIT_1@ASSERT_FAIL:INPUT_LEDGER |
| SOURCE_COMPLETENESS | EXIT_1@ASSERT_FAIL:SOURCE_COMPLETENESS | EXIT_1@ASSERT_FAIL:SOURCE_COMPLETENESS |
| DIMENSIONAL | EXIT_1@ASSERT_FAIL:DIMENSIONAL | EXIT_1@ASSERT_FAIL:DIMENSIONAL |
| COMOVING_CONTINUITY | EXIT_1@ASSERT_FAIL:COMOVING_CONTINUITY | EXIT_1@ASSERT_FAIL:COMOVING_CONTINUITY |
| COMOVING_MOMENTUM | EXIT_1@ASSERT_FAIL:COMOVING_MOMENTUM | EXIT_1@ASSERT_FAIL:COMOVING_MOMENTUM |
| BASE_BALANCE | EXIT_1@ASSERT_FAIL:BASE_BALANCE | EXIT_1@ASSERT_FAIL:BASE_BALANCE |
| TAIL_ODE | EXIT_1@ASSERT_FAIL:TAIL_ODE | EXIT_1@ASSERT_FAIL:TAIL_ODE |
| ZERO_MODE | EXIT_1@ASSERT_FAIL:ZERO_MODE | EXIT_1@ASSERT_FAIL:ZERO_MODE |
| PROJECTOR | EXIT_1@ASSERT_FAIL:PROJECTOR | EXIT_1@ASSERT_FAIL:PROJECTOR |
| ENDPOINT_RESPONSE | EXIT_1@ASSERT_FAIL:ENDPOINT_RESPONSE | EXIT_1@ASSERT_FAIL:ENDPOINT_RESPONSE |
| MOMENT_INTEGRALS | EXIT_1@ASSERT_FAIL:MOMENT_INTEGRALS | EXIT_1@ASSERT_FAIL:MOMENT_INTEGRALS |
| RECONSTRUCTION | EXIT_1@ASSERT_FAIL:RECONSTRUCTION | EXIT_1@ASSERT_FAIL:RECONSTRUCTION |
| CANONICAL_VARIATION | EXIT_1@ASSERT_FAIL:CANONICAL_VARIATION | EXIT_1@ASSERT_FAIL:CANONICAL_VARIATION |
| CHANNEL_UNIQUENESS | EXIT_1@ASSERT_FAIL:CHANNEL_UNIQUENESS | EXIT_1@ASSERT_FAIL:CHANNEL_UNIQUENESS |
| TYPED_DATAFLOW | EXIT_1@ASSERT_FAIL:TYPED_DATAFLOW | EXIT_1@ASSERT_FAIL:TYPED_DATAFLOW |
| PROVENANCE_FORBIDDEN | EXIT_1@ASSERT_FAIL:PROVENANCE_FORBIDDEN | EXIT_1@ASSERT_FAIL:PROVENANCE_FORBIDDEN |
| ANCESTRY | EXIT_1@ASSERT_FAIL:ANCESTRY | EXIT_1@ASSERT_FAIL:ANCESTRY |
| NATIVE_PADDING | EXIT_1@ASSERT_FAIL:NATIVE_PADDING | EXIT_1@ASSERT_FAIL:NATIVE_PADDING |
| PARITY | EXIT_1@ASSERT_FAIL:PARITY | EXIT_1@ASSERT_FAIL:PARITY |
| OUTCOME_REACHABILITY | EXIT_1@ASSERT_FAIL:OUTCOME_REACHABILITY | EXIT_1@ASSERT_FAIL:OUTCOME_REACHABILITY |

Comparator mutations independently perturb one canonical coefficient and inject a symbol into one canonical monomial; they fail `ENGINE_CANONICAL` and `ENGINE_DEPENDENCIES` respectively.

### Outcome reachability through field data

| fixture | computed outcome |
| --- | --- |
| fat_tail | U1_BASE_ILL_POSED(NO_NORMALIZABLE_ZERO_MODE:bound_phase,h) |
| gapped_and_localized | U1_BASE_OK |
| open_h_sign | U1_BASE_UNRESOLVED(Mh,cE) |
| unstable_eos | U1_BASE_UNSTABLE(density_EOS) |

These controls use the same radial solver/classifier. The fat-tail fixture changes the spatial field domain, the unstable fixture changes the EOS action coefficient, and the unresolved fixture removes the positivity stratum of the OPEN `h` coefficient; none supplies an outcome token or restores a P field.

## Dual-engine agreement

- `axis verdicts derived identically`
- `Decision-16 citation, stage-note manifests, and exact P-retired action coverage`
- `parsed action expressions and two active remove/re-derive source probes`
- `15 inline dimensional constructions`
- `continuity and momentum co-moving reductions`
- `6 solved radial channels, residuals, and translated norms`
- `action-derived operator Hessians and density-phase indicial/degree verdict`
- `force-balance operator, translated source block, Gram/projector algebra`
- `five solved endpoint responses and all evaluated reduced moments`
- `canonical L_eff coefficient monomials and expression-derived dependency sets`
- `substitution reconstruction and no-double-count channel partition`
- `66 per-structure ancestry/native-padding ablations`
- `three-way parity tags and four data-driven outcome classes`

Agreement is on source-derived action coverage, solved ODE residuals, norm values, endpoint BVP coefficients, evaluated moments, canonical additive monomials, and dependency sets extracted from those monomials—not on copied verdict strings.

## Provenance, ancestry, and parity

```text
ACTION + primitive field traces -> balanced exterior family -> D_E -> translated zero mode / quotient
ACTION + primitive field traces + endpoint field BVP -> evaluated moments -> L_eff -> momenta / Q_var
BALANCE + CONSTRAINT(E4) + RAYLEIGH(E5) + RETURN -> total generalized force (not L_eff)
```

The machine graph contains `29` typed nodes and `29` edges, all produced by traversal of the parsed quadratic and rigid-embedding expressions. Forbidden import injection is detected by node-type traversal. `66` nonzero additive structures are actually re-derived with their ancestor removed: certified native structures change from `NONZERO_NATIVE_STRUCTURE` to `ABSENT`, while OPEN `h` remainders stay labeled and receive no native certification.

Parity tags: embedding `BODY_CONJUGATION_ONLY`, symmetric background `BODY_PLUS_AMBIENT_POSTULATE`, one-sided control `ONE_SIDED_ASYMMETRY_MAP`.

## Verdict and Phase-A halt

- Axis 1: **`COMPUTATION_VALID`**.
- Axis 2: **`U1_BASE_OK`** in all E1–E5 Phase-A cells.
- Decision-16 baseline status: **`OPERATIVE`**; six channels are present and the independent brane P sector is absent.
- The result is conditional on the declared positive coefficient stratum and the field-level core/surface family; the controls show the classifier can still return unstable, unresolved, and ill-posed outcomes.
- Phase B: `NOT_RUN(upstream)`.
- Phase C: `NOT_RUN(upstream)`.

## Proposed parameter-register rows/edges (not applied)

Proposed OPEN rows: `Mh`, `cE`, `mdot`, `gammaSigma`, `tangentDtN`, `sleeve_core_trace`, `geon_core_bundle`, `throat_surface_functional`, `outer_surface_functional`, `E4_shear_lock`, `E5_rayleigh`, and `return_closure`, with the domains/dimensions/arguments/symmetries in the ledger above.

Retired-row proposal:

| row/knob | proposed status | disposition |
| --- | --- | --- |
| lambdaPu | RETIRED_DECISION_16 | remove the twist–shear coupling row |
| alphaAniso | RETIRED_DECISION_16 | remove the easy-plane wall row |
| T0_polar_kinetic | RETIRED_WITH_P_SECTOR | remove inherited polar kinetic normalization |
| T0_polar_gradient | RETIRED_WITH_P_SECTOR | remove inherited polar gradient normalization |
| T0_polar_quartic | RETIRED_WITH_P_SECTOR | remove inherited polar quartic normalization |

Proposed live edges: `(hbar,m,rho_inf,K_EOS)->density D_E`; `(aB,kappaB)->chiB D_E`; `(rhoBr,muR,ellg)->u D_E`; `(Mh,cE)->h D_E`; `core traces+surface functional->balanced exterior family`; `endpoint trace systems->N_**,U_**`; `ACTION moments->L_eff`; `E4_shear_lock->F_constraint`; `E5_rayleigh->F_Rayleigh`; `return_closure->F_flux`. There is no live P edge.

**HALT: Phase A is complete; no Phase-B computation is included.**
---

## Phase B1 — indexed mechanics remediation 3

The independently recomputed Phase-A payload is `a32c25f4325671d280b54df6c51abd9b25008ef5e6008b98972bac1ed81f7e69` and remains protected. [claim `phase_a_digest` → `source_contract.phase_a_payload_sha256`; recompute `sha256(normalized_phase_a_payload)`]

### Phase-A amendment carried into B1

`PHASE_A_MOMENT_CORRECTION(brane_shear)` changes the brane-shear unit gradient `12π/5 → 8π/3` (`7.539822368616` → `8.377580409573`). The normalized Phase-A payload digest moves `a32c25f4325671d280b54df6c51abd9b25008ef5e6008b98972bac1ed81f7e69` → `b23993cca80dc3e6a790abcf68c1af63aa804fc47b06b153b9f224ccf27f899d`.

The correction closure contains exactly 10 paths: `evaluated_moments.E1.U_XX`, `evaluated_moments.E2.U_XX`, `evaluated_moments.E3.U_XX`, `evaluated_moments.E4.U_XX`, `evaluated_moments.E5.U_XX`, `endpoint_effective_actions.E1.coefficients.GVV`, `endpoint_effective_actions.E2.coefficients.GVV`, `endpoint_effective_actions.E3.coefficients.GVV`, `endpoint_effective_actions.E4.coefficients.GVV`, `endpoint_effective_actions.E5.coefficients.GVV`.

Tilt-profile rows are frozen as `UNCHANGED;UNRESOLVED(tilt_profile)`; byte semantics unchanged = `true`.

### Derived positive content and honest exits

The manifest contains 7 typed field routes and 2 surface variations. It distinguishes rigid-advection from endpoint-response tangents and carries each field-local radial dimension and measure. [claim `field_manifest` → `field_manifest.fields`; recompute `join(indexed_routes,phase_a.tail_channels)`]

Failed indexed-tangent lookups emit 8 leaves: `indexed_density_tilt_profile, indexed_flow_tilt_response, indexed_h_tilt_profile, indexed_phase_tilt_profile, indexed_shear_tilt_profile, indexed_sleeve_surface_normal_profile, indexed_sleeve_tilt_profile, indexed_uw_tilt_profile`. [claim `emitted_leaves` → `indexed_profile_missing_leaves`; recompute `failed_phase_a_indexed_tangent_lookups`]

All 10 endpoint/ambient cells emit a derived isotropic `M_XX(p=0)`, zero reconstruction residual, a binding scalar regression, and executable conditional native-slice residuals for `M_Xp`/`M_pp`. [claim `scalar_regression` → `scalar_regression`; recompute `eV.T*M_XX*eV-GVV`] [claim `native_slice_constraints` → `indexed_cells`; recompute `native_MXp/Mpp_projection-minus-GVP/GPP`]

The authoritative component aggregation is `{'UNRESOLVED': 10}`; no missing indexed or OPEN remainder was promoted to a value. [claim `mechanics_map` → `mechanics_map`; recompute `aggregate(cells.component_findings)`] [claim `cell_count` → `cells`; recompute `active_endpoint-ambient-stratum_product`]

| endpoint | computed source | headline | intrinsic Ω |
| --- | --- | --- | --- |
| E1 | E1 | UNRESOLVED | ZERO_COMPUTED |
| E2 | E2 | UNRESOLVED | ZERO_COMPUTED |
| E3 | E3 | UNRESOLVED | ZERO_COMPUTED |
| E4 | E3 | UNRESOLVED | ZERO_COMPUTED |
| E5 | E2 | UNRESOLVED | ZERO_COMPUTED |

### Traversal, dimensions, and congruence

Expression/root traversal produces 42 OPEN root/block records. The shared finite generator yields an empty control (`STRUCTURAL_ZERO`) and a witness control (`NONZERO_WITNESS`). [claim `open_reachability` → `open_root_reachability`; recompute `typed_ledger-union-traversal`] [claim `finite_controls` → `reachability_analysis.finite_bound_controls`; recompute `same_generator_parity_filter`]

Units were restored termwise on the emitted expressions; computed dimension sets are `{'L_translation': [[2, -2, 1]], 'L_native_Xp_slice': [[2, -2, 1]], 'L_native_pp_slice': [[2, -2, 1]], 'M_XX': [[0, 0, 1]], 'M_Xp_native_slice': [[1, 0, 1]], 'M_pp_native_slice': [[2, 0, 1]], 'Omega_XX_control': [[-2, -1, 1]]}`. [claim `dimensions` → `dimensions.records`; recompute `termwise_units_restoration`]

The produced blocks reduce by exact transverse/longitudinal congruence to 6 conditional pivots; the full signature remains unresolved because their coefficients contain emitted remainders. [claim `derived_congruence` → `derived_congruence`; recompute `produced_blocks-to-covariant-coefficients`]

### Berry/G4 and endpoint machinery

The two traversal-derived Berry DAGs cover 10 production cells and agree after the executed zero-mode quotient. Production winding is `0`. [claim `berry_coverage` → `berry.production_pullback_coverage`; recompute `production_cells-cross-ambient_branches`]

The control contour gives `k=1` and downstream `sigma=-1`; the consumed sheet area leaves a zero total-to-per-area residual. [claim `g4_winding` → `g4_control.computed_winding`; recompute `contour_integral/(2*pi)`] [claim `g4_sigma` → `g4_control.computed_sigma`; recompute `Omega_xy/(rho_mass*Gamma*epsilon_xy)`] [claim `sheet_area` → `g4_control.total_to_per_area_residual`; recompute `Omega_total/sheet_cell_area-Omega_per_area`]

E4 emits its pre-constraint action before differentiating `M_aug`; every stored Hessian residual is zero and `J` comes from the constraint/moduli equations. [claim `e4_action_hessian` → `E4.M_aug_hessian_residual`; recompute `hessian(preconstraint_extended_action)-M_aug`]

E5 consumes and deletes root `E5_rayleigh`, then re-solves to the computed E2 conservative response. [claim `e5_root` → `E5.root_deleted_conservative_solution`; recompute `delete(rayleigh_root)-and-resolve`]

### Records, gates, and halt

Return closure is absent by graph reachability, not by declaration. [claim `closure_absence` → `provenance_graph.global_return_closure_absence`; recompute `reachability(return_closure,mechanics_targets)`]

| gate | computed status |
| --- | --- |
| G1 | INHERITED_PHASE_A_REPRODUCED |
| G2 | KNOWN_COEFFICIENTS_FINITE;FULL_REMAINDERS_UNRESOLVED |
| G3 | DERIVED_BLOCK_CONGRUENCE;FULL_SIGNATURE_UNRESOLVED |
| G4 | PASS |
| G5 | KNOWN_M_AND_OMEGA_COVARIANT;FULL_REMAINDERS_UNRESOLVED |
| G6 | ENDPOINT_MAP_COMPUTED |
| G7 | NOT_RUN(phase_B2) |
| G8 | NOT_RUN(phase_C) |
| G9 | NOT_RUN(phase_B2) |
| G10 | NOT_RUN(phase_C) |
| G11 | NOT_RUN(phase_C) |

Gate rows and candidate ownership are artifact-backed. [claim `gate_statuses` → `gate_evidence`; recompute `engine_check-derived-gates`] [claim `partition` → `partition_ledger`; recompute `computed-candidate-ownership`]

The external unchanged-executable gauntlet passed 35 focused mutation cases, 2 metamorphic controls, and 346 per-key liveness cases. [typed mutation artifact `mechanics_mutations`]

The engines agree on 9 independently recomputed groups: per-field d_f angular Gram contractions, isotropy, scalar regression, and native GVP/GPP residuals; endpoint variation solves plus E4 action-first Hessian/J lift and E5 root-deletion solve; 42 traversal-derived OPEN-root dispositions and two exercised finite controls; route-separated Berry DAGs with 10 production pullback cells and a real quotient; G4 sheet-area reduction, k=1 signed control, sigma=-1, and production k=0; units restored on real expressions and derived-coefficient symmetry-block congruence; component records, closure absence, provenance, and all-UNRESOLVED map reaggregated; typed input roles/sinks cover all 346 declared scalar leaves; independent SymPy/Wolfram Phase-A protection and nonshared representations. [typed agreement artifact `mechanics_dual_engine`]

### Known residual constructs

These residual constructs remain explicit and are not accepted silently or presented as independently eliminated:

| construct | bound artifact field | comparator-side backstop |
| --- | --- | --- |
| sheet_area_X_minus_X_residual | symplectic_mechanics.g4_control.total_to_per_area_residual | B1_C_ENGINE_MATH recomputes/nonzero-checks sheet_cell_area and requires the emitted residual matrix to vanish |
| hand_typed_provenance_edge_atoms | mechanics_provenance_graph.edges | B1_C_RECORD_MAP rebuilds the exact edge union from contraction, reachability, Berry, E4, and E5 records |
| oracle_ban_stamps | indexed_mechanics.cells.*.field_contraction_integrals.*.{oracle_ancestry_forbidden,oracle_paths_consumed} | B1_C_ENGINE_MATH requires both engines' stamps and separately compares every forward-derived contraction tensor and quadrature tooth |

### Typed report-claim bindings

| claim | schema path | type | recomputation |
| --- | --- | --- | --- |
| phase_a_digest | source_contract.phase_a_payload_sha256 | sha256 | sha256(normalized_phase_a_payload) |
| field_manifest | field_manifest.fields | per_field_manifest | join(indexed_routes,phase_a.tail_channels) |
| emitted_leaves | indexed_profile_missing_leaves | derived_string_set | failed_phase_a_indexed_tangent_lookups |
| scalar_regression | scalar_regression | residual_mapping | eV.T*M_XX*eV-GVV |
| native_slice_constraints | indexed_cells | conditional_residual_mapping | native_MXp/Mpp_projection-minus-GVP/GPP |
| g4_winding | g4_control.computed_winding | sympy_expression | contour_integral/(2*pi) |
| g4_sigma | g4_control.computed_sigma | sympy_expression | Omega_xy/(rho_mass*Gamma*epsilon_xy) |
| sheet_area | g4_control.total_to_per_area_residual | canonical_matrix | Omega_total/sheet_cell_area-Omega_per_area |
| berry_coverage | berry.production_pullback_coverage | cell_selection | production_cells-cross-ambient_branches |
| e4_action_hessian | E4.M_aug_hessian_residual | canonical_matrix | hessian(preconstraint_extended_action)-M_aug |
| open_reachability | open_root_reachability | generator_records | typed_ledger-union-traversal |
| finite_controls | reachability_analysis.finite_bound_controls | emptiness_and_witness_records | same_generator_parity_filter |
| cell_count | cells | mapping_length | active_endpoint-ambient-stratum_product |
| dimensions | dimensions.records | computed_dimension_records | termwise_units_restoration |
| derived_congruence | derived_congruence | block_congruence | produced_blocks-to-covariant-coefficients |
| mechanics_map | mechanics_map | component_aggregate_map | aggregate(cells.component_findings) |
| closure_absence | provenance_graph.global_return_closure_absence | graph_predicate | reachability(return_closure,mechanics_targets) |
| e5_root | E5.root_deleted_conservative_solution | root_ablation_solve | delete(rayleigh_root)-and-resolve |
| gate_statuses | gate_evidence | gate_mapping | engine_check-derived-gates |
| partition | partition_ledger | ownership_ledger | computed-candidate-ownership |

Axis 1 is `COMPUTATION_VALID`; every mechanics cell remains honestly `UNRESOLVED`. B2 and Phase C remain `NOT_RUN(phase)`. [claim `mechanics_map` → `mechanics_map`; recompute `aggregate(cells.component_findings)`] [claim `partition` → `partition_ledger`; recompute `computed-candidate-ownership`]
