# Path-A 21 emergent G and mass bridge summary

## Verdicts

- P1 force: `G_FREE_PROFILE_FUNCTIONAL_DERIVED_CONDITIONAL_REDUCED_3D`. Form: `F_12 = -[C_F,12/r^2] rhat in the compact reduced-3D drain lane`; coefficient `C_F,12=(m_GNLS*N_infty,3*Q1*Q2/(4*pi))*I_F,12 = m_GNLS*Theta_Q1*Theta_Q2*J1*J2*I_F,12/(4*pi*N_infty,3)`.
- P2 bridge: `MASS_BRIDGE_FORM_NOT_DERIVED`. EP: `EP_NOT_DERIVED`.
- P3 M-collapse: `RETAIN_L_T_M` because `HBAR_FREE_SUBSTRATE_RELATION_MISSING and MASS_BRIDGE_FORM_NOT_DERIVED`.
- P4 G: `NEWTON_G_FORM_NOT_DERIVED` because P2 does not independently derive masses and P1 inverse-square universality remains profile-conditional.

Dual-engine scope: the Python and Mathematica scripts check dimensions and algebra only. They do not prove the non-algebraic P1/P2/P4 source relations.

## P1 force coefficient

Source-equation chain:

1. Parent continuity gives the drain rate: `partial_t rho + partial_i j^i = 0`, `j^i=rho v^i` (`pde.tex:396-406`; `part01_parent_geometry.tex:213-219`).
2. The stationary pressure response comes from the Euler-like identity with `V_conf`, `h(rho)`, and `Q(rho)` retained (`pde.tex:440-451`; `part01_parent_geometry.tex:275-286`).
3. Projected brane continuity is open unless `S_leak` and the topology BC close (`pde.tex:511-539`), and reduction is a separate profile assumption (`pde.tex:541-565`).
4. In the compact reduced-3D far-field lane, the solved profile returns `Q_i=Theta_Qi*J_i/N_infty,3`, `v_i=-Q_i rhat/(4*pi*r^2)`. The pressure/momentum cross term on a control surface around throat 2 gives `F_12=-C_F,12 rhat/r^2`.

Attractiveness: Attractive for positive compressibility, positive N_infty,3, and inward-positive drains Q1,Q2>0 because the pressure/momentum cross term gives F_12 parallel to the external inflow toward the other drain.
Power law: r^-2 only after compact reduced-3D no-leakage behavior; unreduced bulk lane is r^-3.
Profile dependencies: Q(rho), V_conf, R0 geometry, S_leak/topology BC, W_eff, Theta_Q, I_F,12.
Non-closing conditions: compact reduced-3D source, far-field no-leakage, positive compressibility, solved stationary branch.

## P2 mass bridge and EP

Verdict: `MASS_BRIDGE_FORM_NOT_DERIVED`.

- Angular-rate candidate: `m_defect ?= alpha_J*hbar*J_omega/c_gamma^2`.
- Cycle-rate candidate: `m_defect ?= alpha_J*h*J_nu/c_gamma^2 = 2*pi*alpha_J*hbar*J_nu/c_gamma^2`.
- `2*pi` status: H_2PI_RATE_CLASSIFICATION_UNDETERMINED; alpha_J does not absorb 2*pi.
- `alpha_J`: No independent alpha_J profile functional is derived; alpha_H,omega=H_throat/(hbar*J_omega) is specified as a needed future relation, not accepted as mass bridge.
- EP verdict: `EP_NOT_DERIVED` because the accelerated-throat inertial functional and the far-field source mass functional are not separately available with the same normalization.

No row defines `alpha_J := m_defect*c_gamma^2/(hbar*J)`, and no mass formula is accepted by restatement.

## P3 M-collapse

`RETAIN_L_T_M`. INFLOW_MASS_SOURCE_MISSING is sharpened, not closed. The retained base is `L, T, M`.

## P4 G and m-to-G

Verdict: `NEWTON_G_FORM_NOT_DERIVED`.
Conditional algebraic hand-off only: `If future solves prove universal Theta_Q/alpha_J and I_F,12=1, then G_cond=c_gamma^4*m_GNLS*Theta_Q1*Theta_Q2*I_F,12/(4*pi*N_infty,3*alpha_J1*alpha_J2*hbar^2), with N_infty,3=W_eff*rho_infty,4.`.
Width discipline: W_eff named reduction width; not set to a or xi_h/sqrt(2).

This is not an extracted Newton constant because the independent mass bridge is missing and the inverse-square/factorized/universal force conditions are not all closed.

## P5 profile-solve specification

Rows: 35 total; 16 profile-solve, 14 pathA_22, 5 new-physics, 0 known.

| symbol | definition | dimension | frame | source anchor | closes which output | status | residual if absent | downstream consumer |
|---|---|---|---|---|---|---|---|---|
| R0(w) | Stationary throat surface Sigma0(X)=r-R0(w); domain 0<=w<=L0 with mouth R0(0) and bottom/topology BC. | L | 4D-bulk / reduced throat | pde.tex:2515-2518; part01_parent_geometry.tex:447-461 | C_F, W_eff, branch geometry | profile-solve | STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA | pathA_21 P1/P4; pathA_22 scale map |
| psi0(X) | Background matter field on the stationary branch; rho0(X)=abs(psi0)^2 enters current, Q, h(rho), and drain flux. | L^-2 | 4D-bulk | pde.tex:2519-2522; pde.tex:326-406 | C_F, J-value, pressure response | profile-solve | STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA | pathA_21 P1/P5 |
| rho0(X) | rho0=abs(psi0)^2 on the branch; measure d^4X in bulk and reduced d^3x after W_eff integration. | L^-4 bulk; L^-3 reduced-3D after W_eff | 4D-bulk / reduced-3D | pde.tex:431-443; part01_parent_geometry.tex:266-278 | C_F, c_s, W_eff*rho0 | profile-solve | STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA | pathA_21 P1/P4 |
| A_M0(x,w) | Background localized gauge field entering D_t, D_i, v_i, Maxwell/mixed branch data, and c_gamma reduction. | A0: M L^2 T^-2; Ai: M L T^-1 | 4D-bulk / brane | pde.tex:2523-2526; pde.tex:355-416 | brane c_gamma, lambda_gamma, mixed profile data | profile-solve | BRANE_ZERO_MODE_REDUCTION_UNDERIVED | pathA_21 P4; pathA_22 |
| V_conf(X;Sigma0) | Confinement potential on the stationary surface; domain bulk neighborhood of throat; integrand V_conf*rho in L_psi. | M L^2 T^-2 | 4D-bulk | pde.tex:2527-2530; part01_parent_geometry.tex:466-497 | C_F pressure/Bernoulli profile | profile-solve | PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED | pathA_21 P1 |
| Q0(rho0) | Quantum potential -hbar^2/(2m_GNLS) nabla_4^2 sqrt(rho0)/sqrt(rho0), evaluated on the solved branch. | M L^2 T^-2 | 4D-bulk | pde.tex:440-443; part01_parent_geometry.tex:275-286 | C_F pressure/Bernoulli profile | profile-solve | PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED | pathA_21 P1 |
| S_leak | Projected continuity leakage -[W j^w] + int W'(w) j^w dw; domain transverse w boundary plus support of W'. | L^-4 T^-1 projected; L^-3 T^-1 reduced | brane projection / reduced-3D | pde.tex:511-539; part01_parent_geometry.tex:321-330 | r-power, J conservation, C_F | profile-solve | FORCE_POWER_PROFILE_UNDERDETERMINED | pathA_21 P1/P4 |
| W_eff | Named reduction width N_infty,3/rho_infty,4 from the solved brane localization/reduction kernel; measure dw with profile weight, not a or xi_h/sqrt(2). | L | 4D-bulk to reduced-3D | pde.tex:541-565; part01_parent_geometry.tex:298-389 | G, reduced C_F | profile-solve | W_EFF_REDUCTION_UNDERIVED | pathA_21 P4; pathA_22 scale map |
| N_infty,3 | Asymptotic reduced number density int rho0(x,w) chi_N(w) dw = W_eff*rho_infty,4 in the far field. | L^-3 | reduced-3D | pde.tex:496-509; software/stage1_solver/reports/pathA_19_dimensional_foundation.md:20-28 | C_F, conditional G | profile-solve | W_EFF_REDUCTION_UNDERIVED | pathA_21 P1/P4 |
| J | Number-rate flux lim_{S_R} int n_3 v_brane.n dS in a no-leakage steady region, with throat-source sign convention specified. | T^-1 | reduced-3D / brane | pde.tex:396-406; pde.tex:511-539 | C_F, alpha_J candidate, J-value | profile-solve | STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA | pathA_21 P1/P2 |
| Theta_Q | Dimensionless far-field drain factor Theta_Q=-(N_infty,3/J) lim_{R->infty} int_{S_R} v_brane.n dS; fields psi0,R0,A_M0 and leakage BC. | 1 | reduced-3D | pde.tex:511-539; pde.tex:2515-2566 | C_F | profile-solve | PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED | pathA_21 P1/P4 |
| I_F,12 | Dimensionless control-surface pressure/momentum cross integral around throat 2: normalize int_{partial U2} Pi_cross[v1,v2,rho0].n dS by m_GNLS N_infty,3 Q1 Q2/(4*pi r^2). | 1 | reduced-3D | pde.tex:445-451; em_fields.tex:118-124 for legacy Euler pressure form | C_F, attractiveness | profile-solve | PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED | pathA_21 P1 |
| C_F,12 | G-free force coefficient C_F,12=(m_GNLS N_infty,3 Q1 Q2/(4*pi))*I_F,12 with Qi=Theta_Qi*Ji/N_infty,3. | M L^3 T^-2 | reduced-3D | pde.tex:396-451 plus profile rows J,Theta_Q,I_F,12 | C_F | profile-solve | PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED | pathA_21 P1/P4 |
| alpha_H,omega | Would-be dimensionless profile energy ratio H_throat[psi0,A_M0,R0,wall]/(hbar J_omega), with H_throat from the canonical Hamiltonian over the throat domain and bottom BC. | 1 | 4D-bulk / throat | pde.tex:318-391; brane_bulk_ontology.tex:1267-1297 | alpha_J, m_defect bridge | new-physics | BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING | pathA_21 P2; future profile solve |
| J_omega | Angular-rate version of the throat invariant used in alpha_J*hbar*J_omega/c_gamma^2; classification requires source relation to phase/angular rate. | T^-1 | throat / brane | pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60 | 2pi placement, alpha_J candidate | new-physics | H_2PI_RATE_CLASSIFICATION_UNDETERMINED | pathA_21 P2 |
| J_nu | Cycle-rate version of the throat invariant, with h*J_nu=2*pi*hbar*J_nu and alpha_J not absorbing 2*pi. | T^-1 | throat / brane | pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60 | 2pi placement, alpha_J candidate | new-physics | H_2PI_RATE_CLASSIFICATION_UNDETERMINED | pathA_21 P2 |
| M_inertial | Second derivative of the effective moving-throat action with respect to a slow center velocity after integrating fields over the solved throat/support domain. | M | brane / reduced throat | pde.tex:2806-2879 scope statement; no explicit accelerated-throat source equation in cited parents | EP inertial side | new-physics | INERTIAL_PROFILE_RESPONSE_MISSING | pathA_21 P2 EP |
| M_source | Mass normalization inferred from the far-field drain source after C_F factorization, separately from M_inertial and without using Newton G. | M | reduced-3D | pde.tex:396-451; brane_bulk_ontology.tex:1267-1297 | EP source side | new-physics | SOURCE_MASS_PROFILE_NORMALIZATION_MISSING | pathA_21 P2/P4 |
| C_B/C_E | Bulk transverse Maxwell principal coefficient ratio from the localized Maxwell kinetic operator. | L^2 T^-2 | 4D-bulk | pde.tex:355-416; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:39-50 | bulk c_gamma | profile-solve | BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED | pathA_21 P4 |
| lambda_gamma | Observed brane photon/sound ratio c_gamma/c_s from the zero-mode/profile reduction. | 1 | brane | pde.tex:541-565; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-51 | G conditional form, mass bridge c_gamma | profile-solve | BRANE_ZERO_MODE_REDUCTION_UNDERIVED | pathA_21 P2/P4; pathA_22 |
| mu_eta(w) | Wall inertia density in the reduced wall action, integrated over the finite throat axial coordinate. | M L^-1 | reduced throat | pde.tex:2531-2535 | pathA_22 support/scale map | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| T_w(w) | Axial wall tension function in the reduced wall operator. | M L T^-2 | reduced throat | pde.tex:2531-2535 | pathA_22 support/scale map | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| T_Omega(w) | Angular wall stiffness/tension density entering the grouped l=2 wall operator. | M L^-1 T^-2 | reduced throat | pde.tex:2531-2535 | pathA_22 support/scale map | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| K_eta(w) | Wall restoring stiffness density in the reduced wall operator. | M L^-1 T^-2 | reduced throat | pde.tex:2531-2535 | pathA_22 support/scale map | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| varpi_Aalpha | Stable BdG/support frequency for grouped real lane A and mode alpha. | T^-1 | reduced throat | pde.tex:2537-2544; pde.tex:2602-2609 | pathA_22 B_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| c_Aalpha | Wall/support coupling for lane A and mode alpha; enters B_n=sum c_Aalpha^2/varpi_Aalpha^(2+2n). | M^1/2 L^-1/2 T^-2 | reduced throat | pde.tex:2537-2544; pde.tex:2602-2609 | pathA_22 B_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| Omega_U,Ar | Conservative mixed/gauge frequency for U-family mode r in grouped lane A. | T^-1 | reduced throat / mixed gauge | pde.tex:2545-2549; pde.tex:2611-2624 | pathA_22 Z_n/N_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| Omega_W,Ar | Conservative mixed/gauge frequency for W-family mode r in grouped lane A. | T^-1 | reduced throat / mixed gauge | pde.tex:2545-2549; pde.tex:2611-2624 | pathA_22 Z_n/N_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| R_Ar | Mixed U-W coupling satisfying Delta_r=Omega_U^2 Omega_W^2 - R_r^2. | T^-2 | reduced throat / mixed gauge | pde.tex:2545-2549; pde.tex:2611-2624 | pathA_22 Z_n/N_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| g_U,Ar | Wall-to-U mixed/gauge coupling entering Q_r,H_r,P_r. | M^1/2 L^-1/2 T^-2 | reduced throat / mixed gauge | pde.tex:2545-2549; pde.tex:2619-2624 | pathA_22 Z_n/N_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| g_W,Ar | Wall-to-W mixed/gauge coupling entering Q_r,H_r,P_r. | M^1/2 L^-1/2 T^-2 | reduced throat / mixed gauge | pde.tex:2545-2549; pde.tex:2619-2624 | pathA_22 Z_n/N_n moments | pathA_22 | PATHA_22_BRANCH_PACKET_INCOMPLETE | pathA_22 |
| a | Mouth-radius collective moment a(t)=(1/4*pi) int_{S^2} R(Omega,0,t) dOmega; not an invariant reduction width. | L | brane / reduced throat | part01_parent_geometry.tex:503-510; pde.tex:2551-2563 | pathA_22 scale map | pathA_22 | A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT | pathA_22 |
| c_s(branch) | Branch sound speed c_s^2=5K rho0^4/m_GNLS evaluated on the asymptotic/profile state. | L T^-1 | 4D-bulk / brane | pde.tex:342-352; software/stage1_solver/reports/pathA_20_velocity_constants.md:9-17 | lambda_gamma, pathA_22 target | profile-solve | STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA | pathA_21 P4; pathA_22 |
| mhat | Source-map factor in the isotropic normalization law; in the 3D target lane its squared dimension must convert P0 to G*c_s^5/(a^5*c^5). | L^-1 T^-1 M^-1/2 in the 3D target normalization | PN-facing reduced observable | pde.tex:2077-2092; pde.tex:2551-2563 | pathA_22 scale map | pathA_22 | SCALE_MAP_SOURCE_FACTOR_UNDERIVED | pathA_22 |
| chi_Q | Outgoing-normalization scalar retained when the passive/outgoing branch is not fixed to canonical compact DtN. | 1 | PN-facing reduced observable | pde.tex:1980-1996; pde.tex:2053-2082; pde.tex:2551-2552 | pathA_22 outgoing normalization | pathA_22 | OUTGOING_DTN_BRANCH_UNDERIVED | pathA_22 |

## Residuals carried

- `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`: CARRIED_FROM_pathA_20. The drain strength Theta_Q, W_eff, leakage/topology BC, and profile force integral remain symbolic. Source: pde.tex:2515-2566 and 2847-2879 require the realized branch data set before flux and overlap values are known.
- `PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED`: P1_PROFILE_RESIDUAL. P1 provides a G-free profile-functional C_F; the compact reduced-3D r^-2 lane is conditional on the source profile. Source: pde.tex:396-451 gives current, continuity, quantum potential, and Euler-like force balance but no solved finite-throat branch.
- `FORCE_POWER_PROFILE_UNDERDETERMINED`: P1_SCOPE_RESIDUAL. A universal Newtonian r^-2 force is not promoted until the profile solve proves compact reduced-3D no-leakage behavior. Source: pde.tex:511-539 gives projected open-system brane continuity; pde.tex:541-565 distinguishes reduction from projection.
- `BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING`: BLOCKS_ALPHA_J. alpha_J is not an independently derived profile functional; the mass bridge is rejected. Source: part01_parent_geometry.tex:174-219 and pde.tex:326-406 define the action and current; brane_bulk_ontology.tex:1267-1297 gives drainage scaling only.
- `MASS_BRIDGE_FORM_NOT_DERIVED`: P2_VERDICT. m_defect=alpha_J*hbar*J/c_gamma^2 remains a candidate dimensional form only. Source: No action-level, boundary-source, Noether, or Hamiltonian equation was found that maps inflow J to m_defect without restating the target.
- `H_2PI_RATE_CLASSIFICATION_UNDETERMINED`: CARRIED_FROM_pathA_20. Keep J_omega and J_nu separate; h*J_nu=2*pi*hbar*J_nu and alpha_J cannot absorb 2*pi. Source: pde.tex:429-469 gives angular phase variables; circulation/winding sources do not classify the throat inflow J value.
- `INERTIAL_PROFILE_RESPONSE_MISSING`: BLOCKS_EP_INERTIAL_SIDE. m_inertial cannot be matched to the source mass normalization in this step. Source: The cited parent sources do not provide an accelerated-throat kinetic response functional for m_inertial.
- `SOURCE_MASS_PROFILE_NORMALIZATION_MISSING`: BLOCKS_EP_SOURCE_SIDE. m_source is not separately reduced to the same integral and normalization as m_inertial. Source: P1 supplies C_F as a drain-force profile functional, not an independent mass functional.
- `EP_NOT_DERIVED`: P2_VERDICT. Equivalence of m_inertial and m_source is not claimed. Source: The inertial and source masses are not both available as separately sourced profile integrals.
- `HBAR_FREE_SUBSTRATE_RELATION_MISSING`: P3_BLOCKER. The base remains {L,T,M}; INFLOW_MASS_SOURCE_MISSING is sharpened by MASS_BRIDGE_FORM_NOT_DERIVED. Source: pathA_20b retained HBAR_PROVENANCE_UNDETERMINED; no new hbar-free substrate relation appears in pathA_21.
- `NEWTON_G_FORM_NOT_DERIVED`: P4_VERDICT. The m<->G algebraic form is recorded only as a conditional hand-off to the profile solve/pathA_22. Source: P4 requires a G-free universal inverse-square force plus independently derived P2 masses; P2 failed and P1 remains profile-conditional.
- `W_EFF_REDUCTION_UNDERIVED`: P4_PROFILE_RESIDUAL. Use named W_eff; do not set it to a or xi_h/sqrt(2). Source: part01_parent_geometry.tex:298-389 and pde.tex:496-565 define projection/reduction but do not fix an invariant width.
- `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`: CARRIED_FROM_pathA_20b. lambda_gamma remains symbolic in all pathA_21 forms. Source: pde.tex:541-565 and software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-51 keep the observed brane c_gamma/c_s as a profile/reduction residual.

## Harness summary

- Dimensional checks: 16 consistent, 0 inconsistent, 16 total.
- Algebraic checks: 4 consistent, 0 inconsistent, 4 total.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

