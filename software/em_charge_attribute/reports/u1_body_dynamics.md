# U1 throat-body dynamics — Phase A remediation

Input SHA-256: `c55acef7564261eb881ddf841328ffb23ca6720620ef68c72571840e50723bce`. Phase A Axis 1 is **`COMPUTATION_VALID`** and Axis 2 is **`U1_BASE_UNSTABLE(Pu_coupled_helical)`**.
This report halts after spec 7.0/7.1. Phase B and Phase C are `NOT_RUN(upstream)` by process. The R1 coupled analysis changes the former scalar-only physics verdict, and that changed verdict is reported without repair toward the earlier target.

## Frozen setup and honest scope

The medium-rest lab frame, co-moving steady family, co-moving `Omega_c`, ambient-subtracted exterior-ball IR scheme, fixed `mdot`, and future `C_mdot` ownership of acceleration-like outer flux were declared before residual evaluation. The family is an action-derived exterior solution joined at `r=a` to field-level core traces through the typed throat surface functional. This is a collective-coordinate/effective-action computation, not a nonlinear throat simulation.

No production input contains a tail exponent, eigenvalue sign, class boolean, or verdict. `Mh` and `cE` remain symbolic positive OPEN action coefficients; no collective-coordinate or multipole functional is an OPEN leaf.

## 7.0 — declared-input ledger

| root | status | type | domain | [L,T,M] | arguments | symmetry | source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sleeve_core_trace | OPEN_INPUT | PRIMITIVE_OPEN | scalar_field_trace_on_Sigma | (0, 0, 0) | field_point, outward_normal, wall_side, a, L | BODY_CONJUGATE_SCALAR | nonlinear_throat_core_open |
| geon_core_bundle | OPEN_INPUT | PRIMITIVE_OPEN | field_trace_bundle_on_Sigma | (0, 0, 0) | density_trace, polar_trace, h_trace, phase_flux_trace, wall_side | BODY_CONJUGATE_BUNDLE | spec_7_0_geon_profile_open |
| throat_surface_functional | OPEN_INPUT | ACTION | field_functional_on_Sigma | (2, -2, 1) | field_traces, normal_derivatives, induced_metric | BODY_CONJUGATE_SCALAR | ledger_stage006_f_throat |
| outer_surface_functional | OPEN_INPUT | ACTION | field_functional_on_partial_Omega_c | (2, -2, 1) | field_traces, outward_normal, ambient_branch | AMBIENT_BRANCH_COVARIANT | return_surface_action_open |
| native_continuity | WHITELIST_SOURCED | BALANCE | fields_on_time_x_Omega_c | (-4, -1, 0) | number_density_field, velocity_field | R_W_COVARIANT | ledger_stage004_stage006 |
| native_momentum | WHITELIST_SOURCED | BALANCE | fields_on_time_x_Omega_c | (-3, -2, 1) | momentum_density_field, native_stress_tensor, action_force_density | R_W_COVARIANT | GNLS_wall_polar_shear_h_action |
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
| alphaAniso | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-2, -2, 1) | none | R_W_EVEN | ledger_stage006 |
| rhoBr | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-3, 0, 1) | none | R_W_EVEN | ledger_stage007 |
| muR | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-1, -2, 1) | none | R_W_EVEN | ledger_stage007 |
| lambdaPu | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (-1, -2, 1) | none | R_W_EVEN | ledger_stage007 |
| OmegaW | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (0, -1, 0) | none | R_W_EVEN | ledger_stage007 |
| ellg | WHITELIST_SOURCED | ACTION_COEFFICIENT | coefficient | (1, 0, 0) | none | R_W_EVEN | ledger_stage007 |
| Mh | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, -2, 1) | none | R_W_EVEN | pathA_38_dynamic_normalization_open |
| cE | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (1, -1, 0) | none | R_W_EVEN | pathA_38_dynamic_speed_open |
| mdot | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, -1, 1) | none | R_W_EVEN | spec_7_0_fixed_intake |
| gammaSigma | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, 0, 0) | none | R_W_EVEN | U2_E5_field_Rayleigh |
| tangentDtN | OPEN_INPUT | ACTION_COEFFICIENT | coefficient | (0, 0, 0) | none | R_W_EVEN | U2_E5_field_surface_DtN |

### Action source completeness

Both engines loaded the cited stage-note files, parsed every executable action expression, and parsed the T0/G0 action blocks. The assembled source cover contains the GNLS Berry/flow/EOS terms, Madelung quantum pressure, wall well/gradient/shear/anisotropy and typed throat/mix functionals, all three T0 polar terms, all MacCullagh/P–u/`u_w` terms, and the OPEN `h` kinetic/gradient normalization.

| assembled term | loaded source | parsed source fragment |
| --- | --- | --- |
| bulk_berry | research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md | `i hbar psi* d_t psi` |
| bulk_flow_kinetic | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `(1/2) m_GNLS n \|u\|^2` |
| quantum_pressure | software/stage1_solver/reports/pathA_35_G0_freeze.md | `QP := (hbar^2/(8 m rho))` |
| bulk_EOS | research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md | `U(rho) = (K/4) rho^5` |
| wall_double_well | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_B(chi_B) = a_B chi_B^2 (1-chi_B)^2` |
| wall_gradient | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `(kappa_B/2)\|grad_4 chi_B\|^2` |
| wall_shear_gate | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `chi_B f_shear` |
| wall_anisotropy | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `alpha_aniso chi_B (P·w_hat)^2` |
| throat_source | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_throat` |
| wall_mix | research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md | `f_mix` |
| polar_kinetic | software/stage1_solver/reports/pathA_24_T0_freeze.md | `(1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)` |
| polar_gradient | software/stage1_solver/reports/pathA_24_T0_freeze.md | `(1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)` |
| polar_quartic | software/stage1_solver/reports/pathA_24_T0_freeze.md | `(1/4) m rho c_s^2(rho) (P^i P^i - 1)^2` |
| brane_shear_kinetic | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br (partial_t u^a)(partial_t u^a)` |
| brane_shear_gradient | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) mu_R Omega_u^a Omega_u^a` |
| Pu_coupling | software/stage1_solver/reports/pathA_35_G0_freeze.md | `L_Pu := - lambda_Pu varpi_a Omega_u^a` |
| uw_kinetic | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br (partial_t u_w)^2` |
| uw_gap | software/stage1_solver/reports/pathA_35_G0_freeze.md | `(1/2) rho_br Omega_w^2 u_w^2` |
| h_kinetic | software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md | `dynamic` |
| h_gradient | software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md | `radial branch` |

The mandatory pieces omitted by the first build are present. The immutable source manifest independently discovers quantum pressure, the MacCullagh shear gradient, and `L_Pu`; deleting any of the three changes its derived operator entry, and their joint SOURCE_COMPLETENESS ablation fails before classification.

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
| wall_anisotropy | (-2, -2, 1) |
| polar_kinetic | (-2, -2, 1) |
| polar_gradient | (-2, -2, 1) |
| polar_quartic | (-2, -2, 1) |
| brane_shear_kinetic | (-2, -2, 1) |
| brane_shear_gradient | (-2, -2, 1) |
| Pu_coupling | (-2, -2, 1) |
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

`B_density=e''(rho_inf)` is differentiated in-engine from `K rho^5/4`, giving `5 K rho_inf^3`; `B_chi` and the radial-polar Hessian are likewise differentiated from their potentials. The translated norm uses the full measure: `S_(d-1) integral r^(d-1) |f'|^2 dr`. For the brane shear channel the four-dimensional measure factorizes as `d^3x g_ell(w)^2 dw`, so its reported radial dimension is three while the `w` integral is finite.

| channel | d | solved radial equation | tail class | nu or gap | translated norm | value at R=a=1 | normalizable |
| --- | --- | --- | --- | --- | --- | --- | --- |
| density_EOS | 4 | Eq(-5*f(r) + Derivative(f(r), (r, 2))/4 + 3*Derivative(f(r), r)/(4*r), 0) | EXPONENTIAL_GAP | 2*sqrt(5) | Integral[1,inf] r^3 (d/dr(r^-1.0 K_1(gap*r)))^2 dr | 0.000204528635461 | True |
| wall_chiB | 4 | Eq(-3*f(r) + 2*Derivative(f(r), (r, 2)) + 6*Derivative(f(r), r)/r, 0) | EXPONENTIAL_GAP | sqrt(6)/2 | Integral[1,inf] r^3 (d/dr(r^-1.0 K_1(gap*r)))^2 dr | 0.409684761529 | True |
| bound_phase | 4 | Eq(Derivative(f(r), (r, 2)) + 3*Derivative(f(r), r)/r, 0) | POWER_LAW | 2 | 2/R**2 | 2 | True |
| polar_tangent | 4 | Eq(5*Derivative(f(r), (r, 2)) + 15*Derivative(f(r), r)/r, 0) | POWER_LAW | 2 | 2/R**2 | 2 | True |
| polar_radial | 4 | Eq(-10*f(r) + 5*Derivative(f(r), (r, 2)) + 15*Derivative(f(r), r)/r, 0) | EXPONENTIAL_GAP | sqrt(2) | Integral[1,inf] r^3 (d/dr(r^-1.0 K_1(gap*r)))^2 dr | 0.238971420511 | True |
| brane_shear | 3 | Eq(5*Derivative(f(r), (r, 2))/4 + 5*Derivative(f(r), r)/(2*r) - 5*f(r)/(2*r**2), 0) | POWER_LAW | 2 | 4/(3*R**3) | 1.33333333333 | True |
| uw | 3 | 196/75*f(r)=0 | ALGEBRAIC_GAP | algebraic | 0 | 0 | True |
| h | 4 | Eq(Mh*(Derivative(f(r), (r, 2)) + 3*Derivative(f(r), r)/r), 0) | POWER_LAW | 2 | 2/R**2 | 2 | True |

### R1 coupled operator result

The displayed scalar channel rows remain the diagonal solves, but the verdict is taken from the full action-derived block. The drain background produces a density–phase mixed Hessian; its radial degree is two powers below the phase diagonal in `d=4`, so the computed indicial degree test leaves the scalar decay exponents unchanged.

For the localized P–u term, each engine independently solves the bulk polar half-space profile and substitutes it back into the quadratic action. The resulting surface Hessian is

```text
[['10*k', 'k/10'], ['k/10', '5*k**2/4']]
det = k**2*(1250*k - 1)/100
witness k = 1/2500; det(witness) = -1/1250000000
```

Its computed class is **`UNSTABLE_HELICAL`**. Thus the normalizable translational scalar tails still exist, but the declared homogeneous base has a negative long-wavelength coupled mode and Axis 2 changes honestly to **`U1_BASE_UNSTABLE(Pu_coupled_helical)`**.

The bound-flow phase is selected from the solved continuity family by the `mdot` surface normalization; its computed normalization residual is `0`. No exponent is supplied in YAML.

### Force-balance operator and quotient

```text
D_E=(-A_n Laplacian+U''(n0)) on delta_n, plus translated throat-source block; analogous Hessian/balance blocks for chiB,P,u,h and endpoint traces
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
| I_Pu_cross | 0.0475998886908 | 0.0475998886908 | 0.0475998886908 | 0.0475998886908 | 0.0475998886908 |
| I_density_grad | 0.756796365754 | 0.756796365754 | 0.756796365754 | 0.756796365754 | 0.756796365754 |
| I_density_l2 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 | 0.0217853543438 |
| I_h_grad | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 | 0.0771062843835 |
| I_polar_grad | 0.173489139863 | 0.173489139863 | 0.173489139863 | 0.173489139863 | 0.173489139863 |
| I_polar_quartic | 0.154212568767 | 0.154212568767 | 0.154212568767 | 0.154212568767 | 0.154212568767 |
| I_shear_grad | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 | 0.0623125815588 |
| I_wall_grad | 1.28479139854 | 1.28479139854 | 1.28479139854 | 1.28479139854 | 1.28479139854 |
| I_wall_l2 | 0.204279992179 | 0.204279992179 | 0.204279992179 | 0.204279992179 | 0.204279992179 |
| N_UU | 10.3295924476 | 2.5823981119 | 0 | 0 | 4.11270810413 |
| N_UW | 1.14773249418 | 0.286933123544 | 0 | 0 | 0.669510621603 |
| N_WW | 0.127525832686 | 0.127525832686 | 0.127525832686 | 0.127525832686 | 0.127525832686 |
| P_PP | 0.154212568767 | 0.154212568767 | 0.154212568767 | 0.154212568767 | 0.154212568767 |
| P_XP | 0 | 0 | 0 | 0 | 0 |
| P_XX | 0.736306821818 | 0.736306821818 | 0.736306821818 | 0.736306821818 | 0.736306821818 |
| U_PP | 0.034618100866 | 0.034618100866 | 0.034618100866 | 0.034618100866 | 0.034618100866 |
| U_XP | 0.0634665182543 | 0.0634665182543 | 0.0634665182543 | 0.44426562778 | 0.0634665182543 |
| U_XX | 0.209439510239 | 0.209439510239 | 0.209439510239 | 10.2625360017 | 0.209439510239 |

For each endpoint the directly embedded action reduces to

```text
L_eff = A_X V + A_p pdot + C_Xp p V + G_VV V^2/2 + G_Vp V pdot + G_pp pdot^2/2 - K_pp p^2/2.
```

- `E1`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 11.945654675034*a**2*m_GNLS*polar_radial**2 + 9.8696044010894*a**2*m_GNLS*polar_tangent**2 + 2.2999402325143*density_delta*m_GNLS + 9.8696044010894*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=2.2999402325143*density_delta*m_GNLS*tilt_phase + 9.8696044010894*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 9.8696044010894*a**2*m_GNLS*tilt_polar**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=55.516524756128*K_EOS*a**2*rho_inf**5*tilt_polar**2 + 2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 8.3775804095728*lambdaPu*tilt_polar*tilt_shear + 7.5398223686155*muR*tilt_shear**2`.
- `E2`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 11.945654675034*a**2*m_GNLS*polar_radial**2 + 9.8696044010894*a**2*m_GNLS*polar_tangent**2 + 0.57498505812856*density_delta*m_GNLS + 2.4674011002723*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=0.57498505812856*density_delta*m_GNLS*tilt_phase + 2.4674011002723*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 9.8696044010894*a**2*m_GNLS*tilt_polar**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=55.516524756128*K_EOS*a**2*rho_inf**5*tilt_polar**2 + 2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 8.3775804095728*lambdaPu*tilt_polar*tilt_shear + 7.5398223686155*muR*tilt_shear**2`.
- `E3`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 11.945654675034*a**2*m_GNLS*polar_radial**2 + 9.8696044010894*a**2*m_GNLS*polar_tangent**2 + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 9.8696044010894*a**2*m_GNLS*tilt_polar**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=55.516524756128*K_EOS*a**2*rho_inf**5*tilt_polar**2 + 2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 8.3775804095728*lambdaPu*tilt_polar*tilt_shear + 7.5398223686155*muR*tilt_shear**2`.
- `E4`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 11.945654675034*a**2*m_GNLS*polar_radial**2 + 9.8696044010894*a**2*m_GNLS*polar_tangent**2 + 7.5398223686155*rhoBr*shear_transverse**2 + 15.079644737231*rhoBr*shear_transverse + 7.5398223686155*rhoBr`; `G_Vp=4.1887902047864*rhoBr*shear_transverse*tilt_shear + 4.1887902047864*rhoBr*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 9.8696044010894*a**2*m_GNLS*tilt_polar**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=55.516524756128*K_EOS*a**2*rho_inf**5*tilt_polar**2 + 2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 8.3775804095728*lambdaPu*tilt_polar*tilt_shear + 7.5398223686155*muR*tilt_shear**2`.
- `E5`: `G_VV=9.8696044010894*Mh*h_scalar**2/cE**2 + 11.945654675034*a**2*m_GNLS*polar_radial**2 + 9.8696044010894*a**2*m_GNLS*polar_tangent**2 + 0.91571694442697*density_delta*m_GNLS + 3.9295647152485*m_GNLS*rho_inf + 7.5398223686155*rhoBr*shear_transverse**2`; `G_Vp=1.3416318023*density_delta*m_GNLS*tilt_phase + 5.7572692339688*m_GNLS*rho_inf*tilt_phase + 4.1887902047864*rhoBr*shear_transverse*tilt_shear`; `G_pp=9.8696044010894*Mh*tilt_h**2/cE**2 + 9.8696044010894*a**2*m_GNLS*tilt_polar**2 + 2.2999402325143*density_delta*m_GNLS*tilt_phase**2 + 9.8696044010894*m_GNLS*rho_inf*tilt_phase**2 + 4.1887902047864*rhoBr*tilt_shear**2`; `K_pp=55.516524756128*K_EOS*a**2*rho_inf**5*tilt_polar**2 + 2.7231692929794*K_EOS*density_delta**2*rho_inf**3 + 11.103304951226*Mh*tilt_h**2 + 3.6770398592144*aB*wall_delta**2 + 4.7299772859649*density_delta**2*hbar**2/(m_GNLS*rho_inf) + 11.563122586834*kappaB*wall_delta**2 + 8.3775804095728*lambdaPu*tilt_polar*tilt_shear + 7.5398223686155*muR*tilt_shear**2`.

Direct differentiation computes both canonical momenta and `Q_p`. Reconstruction independently iterates over the parsed action expressions, substitutes the rigid field embedding, performs the moment reductions, and compares that sum with the claimed coefficient form; the residual is zero for E1–E5. `K_pp` now includes the action-derived EOS Hessian, wall-well Hessian, and P–u cross term. E4 multiplier reaction and E5 Rayleigh loss remain outside `L_eff` in their typed channels.

## Gates and able-to-fail evidence

| gate | production result |
| --- | --- |
| G1 | PASS |
| G10 | NOT_RUN(phase_C) |
| G11 | NOT_RUN(phase_C) |
| G2 | PASS |
| G3 | CLASSIFIED_BY_AXIS2 |
| G4 | NOT_RUN(phase_B) |
| G5 | PASS |
| G6 | PASS_ENDPOINT_MAP_REPORTED |
| G7 | NOT_RUN(phase_B) |
| G8 | NOT_RUN(phase_C) |
| G9 | NOT_RUN(phase_B;tolerance_predeclared) |

G1 passes because every diagonal translated tail is normalizable, while G2 passes because the displayed norms and reduced coefficients are finite in the predeclared ambient-subtracted exterior scheme. G3 is `CLASSIFIED_BY_AXIS2`: the computed P–u determinant supplies the negative internal-mode witness. G5 passes with body-only, symmetric-postulate, and one-sided-asymmetry tags kept distinct.

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
| fat_tail | U1_BASE_ILL_POSED(NO_NORMALIZABLE_ZERO_MODE:bound_phase,polar_tangent,h) |
| gapped_and_localized | U1_BASE_OK |
| open_h_sign | U1_BASE_UNRESOLVED(Mh,cE) |
| unstable_eos | U1_BASE_UNSTABLE(density_EOS,Pu_coupled_helical) |

These controls use the same radial solver/classifier. The fat-tail fixture changes the spatial field domain, the unstable fixture changes the EOS action coefficient, and the unresolved fixture removes the positivity stratum of the OPEN `h` coefficient; none supplies an outcome token.

## Dual-engine agreement

- `axis verdicts derived identically`
- `stage-note action manifests and assembled coverage`
- `parsed action expressions and three remove/re-derive source probes`
- `20 inline dimensional constructions`
- `continuity and momentum co-moving reductions`
- `8 solved radial channels, residuals, and translated norms`
- `action-derived operator Hessians and coupled indicial/degree verdict`
- `force-balance operator, translated source block, Gram/projector algebra`
- `five solved endpoint responses and all evaluated reduced moments`
- `canonical L_eff coefficient monomials and expression-derived dependency sets`
- `substitution reconstruction and no-double-count channel partition`
- `86 per-structure ancestry/native-padding ablations`
- `three-way parity tags and four data-driven outcome classes`

Agreement is on source-derived action coverage, solved ODE residuals, norm values, endpoint BVP coefficients, evaluated moments, canonical additive monomials, and dependency sets extracted from those monomials—not on copied verdict strings.

## Provenance, ancestry, and parity

```text
ACTION + primitive field traces -> balanced exterior family -> D_E -> translated zero mode / quotient
ACTION + primitive field traces + endpoint field BVP -> evaluated moments -> L_eff -> momenta / Q_var
BALANCE + CONSTRAINT(E4) + RAYLEIGH(E5) + RETURN -> total generalized force (not L_eff)
```

The machine graph contains `34` typed nodes and `35` edges, all produced by traversal of the parsed quadratic and rigid-embedding expressions. Forbidden import injection is detected by node-type traversal. `86` nonzero additive structures are actually re-derived with their ancestor removed: certified native structures change from `NONZERO_NATIVE_STRUCTURE` to `ABSENT`, while OPEN `h` remainders stay labeled and receive no native certification.

Parity tags: embedding `BODY_CONJUGATION_ONLY`, symmetric background `BODY_PLUS_AMBIENT_POSTULATE`, one-sided control `ONE_SIDED_ASYMMETRY_MAP`.

## Verdict and Phase-A halt

- Axis 1: **`COMPUTATION_VALID`**.
- Axis 2: **`U1_BASE_UNSTABLE(Pu_coupled_helical)`** in all E1–E5 Phase-A cells.
- The result is conditional on the declared positive coefficient stratum and the field-level core/surface family; the controls show the classifier does not force OK.
- Phase B: `NOT_RUN(upstream)`.
- Phase C: `NOT_RUN(upstream)`.

## Proposed parameter-register rows/edges (not applied)

Proposed OPEN rows: `Mh`, `cE`, `mdot`, `gammaSigma`, `tangentDtN`, `sleeve_core_trace`, `geon_core_bundle`, `throat_surface_functional`, `outer_surface_functional`, `E4_shear_lock`, `E5_rayleigh`, and `return_closure`, with the domains/dimensions/arguments/symmetries in the ledger above.

Proposed edges: `(hbar,m,rho_inf,K_EOS)->density D_E`; `(aB,kappaB)->chiB D_E`; `(m,rho_inf,K_EOS,a)->P D_E`; `(rhoBr,muR,ellg)->u D_E`; `(Mh,cE)->h D_E`; `core traces+surface functional->balanced exterior family`; `endpoint trace systems->N_**,U_**`; `ACTION moments->L_eff`; `E4_shear_lock->F_constraint`; `E5_rayleigh->F_Rayleigh`; `return_closure->F_flux`.

**HALT: Phase A is complete; no Phase-B computation is included.**

