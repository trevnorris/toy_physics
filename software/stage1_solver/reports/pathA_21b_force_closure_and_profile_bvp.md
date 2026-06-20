# Path-A 21b force closure and stationary-throat BVP

## Verdicts

- P1b force: `G_FREE_DRAIN_FORCE_FORM_DERIVED_WITH_ATTRACTIVE_SIGN_PROFILE_RESIDUAL`; supersedes pathA_21 `G_FREE_PROFILE_FUNCTIONAL_DERIVED_CONDITIONAL_REDUCED_3D`.
- P2 bridge: `MASS_BRIDGE_FORM_NOT_DERIVED`. EP: `EP_NOT_DERIVED`. Carried verbatim.
- P3 M-collapse: `RETAIN_L_T_M` because `HBAR_FREE_SUBSTRATE_RELATION_MISSING and MASS_BRIDGE_FORM_NOT_DERIVED`. Carried verbatim.
- P4 G: `NEWTON_G_FORM_NOT_DERIVED`. Carried verbatim.

Dual-engine scope: Python and Mathematica check dimensions and algebra only. The derivation is the source-equation chain below; exit 0 is not treated as proof.

## P1b drain field

Continuity source chain: parent continuity gives `partial_t rho + partial_i j^i=0`, `j^i=rho v^i` (`pde.tex:396-406`; `part01_parent_geometry.tex:213-219`). In a stationary no-leakage region outside the localized throat source, Gauss gives the far-field radial solution rather than inserting the power.

- Reduced-3D lane: `int_{S_r} N_infty,3 v.n dS=-Theta_Q J => v_r(r)=-Theta_Q J/(4*pi*N_infty,3*r^2), or -Theta_Q J/(4*pi*r^2*N0(r)) if the asymptotic density is not constant.`.
- Bulk-4D lane: `int_{S_R} rho_infty,4 v.n dS=-Theta_Q4 J => v_R(R)=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3), or -Theta_Q4 J/(2*pi^2*R^3*rho0(R)) if the asymptotic density is not constant.`.
- Inverse-square status: P1_INVERSE_SQUARE_FIELD_ASSUMED_NOT_SOLVED resolved for the compact reduced-3D lane by continuity/Gauss; unreduced bulk lane gives R^-3.

The reduced-3D `r^-2` power is therefore the area power of the enclosing two-sphere. The unreduced four-spatial bulk lane is `R^-3` from the three-sphere area `Omega3 R^3`.

## P1b force

Control surface: `F_12=-int_{partial U2} Pi_cross[v1,v2,rho0,A0,Sigma0].n_2 dS; partial U2 is a small closed reduced-3D surface enclosing throat 2 but excluding throat 1.`

The cross stress used by the option-C solve is

```text
Pi_cross,ij=m_GNLS*rho_asym(v1_i v2_j+v2_i v1_j)+delta_ij P_cross+Pi_Q,cross,ij+Pi_V,cross,ij+Pi_EM,cross,ij. Here P_cross=K[(rho0+drho1+drho2)^5-(rho0+drho1)^5-(rho0+drho2)^5+rho0^5], Pi_Q,cross=cross[(hbar^2/(4m_GNLS))*((partial_i rho partial_j rho)/rho-partial_i partial_j rho)], and partial_j Pi_V,cross,ij=[rho partial_i V_conf]_cross for the selected V_conf branch.
```

The pressure cross term is the explicit EOS cross-difference. The quantum term uses the displayed representative of the Madelung quantum stress, with the cross operation meaning `F[rho0+drho1+drho2]-F[rho0+drho1]-F[rho0+drho2]+F[rho0]`. The confinement term is the branch-selected stress-divergence representative whose divergence equals the cross body force from `rho partial_i V_conf`; absent that representative, the row would fall back to `PI_CROSS_STRESS_TENSOR_UNDERIVED`. `Pi_EM,cross` is retained only when gauge/mixed fields are active. The anchors are the Euler identity and EOS (`pde.tex:342-352`, `pde.tex:440-451`; `part01_parent_geometry.tex:275-286`).

Force form: `F_12=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*I_F,12*rhat_12 in the compact reduced-3D lane, with Qi=Theta_Qi*Ji/N_infty,3.`

Sign verdict: `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`. Bernoulli gives delta h=-(m_GNLS/2)delta(v^2) where V_conf,Q,qA0 are asymptotically fixed. From P=K*rho^5 and h=(5K/4)rho^4, dh/drho=5K*rho^3>0 for a stable K>0,rho>0 branch, so higher entrained speed lowers rho and pressure. The remaining finite-throat traction orientation is the profile sign residual.

## Stationary BVP

### G1 Closed Core

- GNLS: `BVP-G1.GNLS: [-hbar^2/(2m_GNLS)D_i0D_i0+V_conf(X;Sigma0)+h(abs(psi0)^2)+q_*A_00-mu]psi0=0 on the 4D bulk throat exterior/interior domain, with D_i0=partial_i-i*q_*A_i0/hbar.` Anchor: `pde.tex:382-391`; action anchor `pde.tex:318-391`.
- Maxwell: `BVP-G1.Maxwell: partial_M(Z(w)F_0^{MN})+xi^-1 partial^N(partial.A_0)=mu0[J_psi^N(psi0,A0)+J_ext,0^N], plus Bianchi identities and a gauge condition such as partial.A_0=0.` Anchor: `pde.tex:410-426`; source bookkeeping `pde.tex:370-374`.

G1 domain and BCs: the fields live on the 4D bulk throat branch domain with measure `d^4X` and stationary time factored by `psi=e^{-i mu t/hbar}psi0`. Far-field BCs fix the asymptotic density, chemical potential/gauge reference, finite energy, and the Gauss flux lane. Throat BCs require regular finite fields and branch-declared mouth/bottom flux data; the value selector for those branch data is not part of G1.

### Boundary Conditions

| condition | fields/domain | frame and measure | source anchor | status |
|---|---|---|---|---|
| Stationary matter asymptotic | `psi0 -> sqrt(rho_infty,4) exp(i theta_infty)`, `rho0 -> rho_infty,4`; `V_conf`, `Q`, and gauge reference approach branch constants | 4D bulk, `d^4X`; `psi0:L^-2`, `rho0:L^-4` | `pde.tex:382-406`; `pde.tex:2519-2530` | closed equation, branch value supplied by option C |
| Reduced-3D Gauss flux | `int_{S_r} N_infty,3 v.n dS=-Theta_Q J`, so `v_r=-Theta_Q J/(4*pi*N_infty,3*r^2)` | reduced-3D, `dS=r^2 dOmega_2`; `J:T^-1`, `N_infty,3:L^-3` | `pde.tex:396-406`; `pde.tex:511-539` | closed conservation law; `J_VALUE_BRANCH_PARAMETER` remains |
| Bulk-4D Gauss flux | `int_{S_R} rho_infty,4 v.n dS=-Theta_Q4 J`, so `v_R=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3)` | 4D bulk, `dS=R^3 dOmega_3`; `rho_infty,4:L^-4` | `pde.tex:396-406` | closed conservation law; compact reduced lane still branch-selected |
| Throat regularity and mouth flux | finite `psi0`, `rho0`, and `A_M0`; mouth inflow oriented into the throat where the branch declares a drain | 4D bulk/brane mouth, mouth `d^2S`; volumetric flux from brane mouth data | `brane_bulk_ontology.tex:1267-1289`; `part01_parent_geometry.tex:447-461` | branch BC; no `J` value selector |
| Bottom/topology | candidate `R0(L0)=0` or equivalent regular bottom condition | reduced throat, `0<=w<=L0` | `part01_parent_geometry.tex:447-461` | branch assumption, not closure |
| Stationary Maxwell far field and gauge | finite-energy `A_M0`; gauge condition from the gauge-fixed Maxwell equation, e.g. `partial.A_0=0`; optional zero-mode BCs only in the reduced brane lane | 4D bulk, `d^4X`; `A0:ML^2T^-2`, `Ai:MLT^-1` | `pde.tex:355-426`; `pde.tex:541-565` | G1 closed for bulk field; G6 brane cone residual |
| Projection kernel normalization | `int W(w)dw=1`; `rho_brane=int W rho dw`, `j_brane=int W j_xyz dw` | brane projection, `dw`; `W:L^-1` | `pde.tex:496-539`; `part01_parent_geometry.tex:298-330` | formula closed; kernel shape is `W_KERNEL_UNDERDECLARED` |

### Gap Table

| gap | status | equation or residual | note |
|---|---|---|---|
| G1 | CLOSED | `BVP-G1.GNLS plus BVP-G1.Maxwell` | Stationary gauged GNLS and stationary localized Maxwell are transcribed from the parent action for fixed branch geometry and source data. |
| G2 | NAMED RESIDUAL | `R0_FREE_BOUNDARY_CONDITION_UNDERIVED` | The parent supplies Sigma0=r-R0(w), a0=R0(0), and candidate bottom/regularity data, but not an Euler-Lagrange selector for R0(w). |
| G3 | CLOSED AS FUNCTIONAL | `BVP-G3.Pi_cross stress integral` | Pi_cross is written from the Euler stress balance/action terms; its numerical value and attractive orientation are profile integrals. |
| G4 | NAMED RESIDUAL | `J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED` | Continuity and no-leakage conserve flux in a chosen lane but do not select the value of J. |
| G5 | NAMED RESIDUAL | `W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED` | Projection formulas are source-anchored; the kernel shape is not selected by the parent. |
| G6 | NAMED RESIDUAL | `BRANE_ZERO_MODE_REDUCTION_UNDERIVED` | The brane photon cone requires a solved zero-mode/profile reduction; lambda_gamma remains symbolic. |

### G2 Branch Choices

`R0_FREE_BOUNDARY_CONDITION_UNDERIVED`: candidate option-C assumptions are a prescribed analytic `R0(w)`, the parent bottom cap `R0(L0)=0`, an equivalent regular bottom condition, or a future free-boundary stationarity equation. These are branches, not closure.

### G4 Branch Choices

`J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED`: no-leakage conserves the flux in a chosen region, but it does not fix the value. Candidate selectors are sonic choking, regularity at the throat bottom, global topology/outflow balance, or an energy extremum, all explicitly downstream assumptions unless derived.

### G5 Projection Formulas

`rho_brane=int W(w)rho dw, j_brane=int W(w)j_xyz dw, N_infty,3=int chi_N(w)rho_infty,4(w)dw, W_eff=N_infty,3/rho_infty,4 for constant far-field rho.`. The formulas are source-anchored (`pde.tex:496-565`; `part01_parent_geometry.tex:298-390`), but the shape of `W(w)` or `chi_N(w)` is `W_KERNEL_UNDERDECLARED`.

### G6 Brane Photon Residual

`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`: option C must solve the localized Maxwell zero-mode/profile reduction and compute the observed brane photon cone before `lambda_gamma=c_gamma/c_s` can be numerical.

### alpha_H,omega Anchor Fix

pde.tex:318-391 is action-level, not a canonical Hamiltonian; canonical Hamiltonian must be constructed before alpha_H,omega can be used.

## Machine Table

Rows: 35 total; 9 closed profile-solve, 6 branch-residual, 14 pathA_22, 6 new-physics, 0 known.

| symbol | definition | dimension | frame | source anchor | closes which output | status | residual if absent | downstream consumer |
|---|---|---|---|---|---|---|---|---|
| R0(w) | Stationary throat surface Sigma0(X)=r-R0(w). Parent gives a0=R0(0) and candidate bottom/regularity data, but no free-boundary equation selecting R0(w). | L | 4D-bulk / reduced throat | pde.tex:2515-2518; part01_parent_geometry.tex:447-461 | branch geometry for BVP-G1; C_F and W_eff only after a branch selection | branch-residual | R0_FREE_BOUNDARY_CONDITION_UNDERIVED | option C branch realization; pathA_22 scale map |
| psi0(X) | Background matter field solved by BVP-G1.GNLS: [-hbar^2/(2m_GNLS)D_i0D_i0+V_conf(X;Sigma0)+h(abs(psi0)^2)+q_*A_00-mu]psi0=0. | L^-2 | 4D-bulk | pde.tex:382-391; pde.tex:396-406; pde.tex:2519-2522 | stationary density/current profile entering C_F, c_s(branch), and Pi_cross | profile-solve | STATIONARY_GNLS_BVP_NOT_SOLVED | option C; pathA_21b P1b/P5b |
| rho0(X) | rho0=abs(psi0)^2 from BVP-G1.GNLS; it enters h(rho0), P(rho0), Q(rho0), current, and the branch sound speed. | L^-4 bulk; L^-3 reduced-3D after W_eff | 4D-bulk / reduced-3D | pde.tex:431-443; pde.tex:342-352; part01_parent_geometry.tex:266-278 | pressure response, Q0, c_s(branch), reduced density after G5 branch data | profile-solve | STATIONARY_GNLS_BVP_NOT_SOLVED | option C; pathA_21b P1b/P5b |
| A_M0(x,w) | Background localized gauge field solved by BVP-G1.Maxwell: partial_M(Z F_0^{MN})+xi^-1 partial^N(partial.A_0)=mu0 J_tot,0^N with the gauge condition and finite-energy BCs. | A0: M L^2 T^-2; Ai: M L T^-1 | 4D-bulk / brane | pde.tex:355-416; pde.tex:2523-2526; part01_parent_geometry.tex:333-390 | stationary gauge background and mixed-sector branch data; brane lambda_gamma remains G6 residual | profile-solve | STATIONARY_MAXWELL_BVP_NOT_SOLVED | option C; pathA_21b P5b; pathA_22 |
| V_conf(X;Sigma0) | Confinement profile entering BVP-G1.GNLS for a selected Sigma0; parent promotes V_conf(X;a,L) to V_conf(X;Sigma) and gives the smooth-wall variation. | M L^2 T^-2 | 4D-bulk | pde.tex:318-334; pde.tex:2527-2530; part01_parent_geometry.tex:466-497 | stationary GNLS potential and pressure/Bernoulli profile for a chosen branch | profile-solve | V_CONF_BRANCH_PROFILE_UNDERDECLARED | option C; pathA_21b P1b |
| Q0(rho0) | Quantum potential Q0=-(hbar^2/(2m_GNLS))*nabla_4^2(sqrt(rho0))/sqrt(rho0), evaluated from the BVP-G1 density. | M L^2 T^-2 | 4D-bulk | pde.tex:440-443; part01_parent_geometry.tex:275-286 | quantum-stress contribution to Pi_cross and Bernoulli balance | profile-solve | STATIONARY_GNLS_BVP_NOT_SOLVED | option C; pathA_21b P1b |
| S_leak | Exact projected leakage S_leak=-[W j^w]+int W'(w)j^w dw. A compact reduced-3D force lane additionally assumes the far-field no-leakage branch S_leak=0. | L^-4 T^-1 projected; L^-3 T^-1 reduced | brane projection / reduced-3D | pde.tex:511-539; part01_parent_geometry.tex:321-330 | far-field continuity equation, drain r-power in the selected compact lane, and flux conservation | profile-solve | NO_LEAKAGE_BRANCH_BC_UNDERIVED | option C; pathA_21b P1b/P5b |
| W_eff | Reduction width W_eff=N_infty,3/rho_infty,4 after a projection/reduction kernel is selected; formula is declared but the kernel shape is not parent-selected. | L | 4D-bulk to reduced-3D | pde.tex:541-565; part01_parent_geometry.tex:298-389 | G reduction and reduced C_F only after G5 branch realization | branch-residual | W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED | option C branch realization; pathA_22 scale map |
| N_infty,3 | Asymptotic reduced density N_infty,3=int chi_N(w) rho_infty,4(w) dw = W_eff*rho_infty,4 only after the G5 kernel branch is selected. | L^-3 | reduced-3D | pde.tex:496-509; software/stage1_solver/reports/pathA_19_dimensional_foundation.md:20-28 | C_F normalization and conditional G denominator | branch-residual | W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED | option C branch realization; pathA_21b P1b/P4 |
| J | Conserved number-rate flux J in a stationary no-leakage region: int_{S_R} n v.n dS=-J. The parent does not select its value. | T^-1 | reduced-3D / brane | pde.tex:396-406; pde.tex:511-539 | P1b force coefficient after branch-selected value; alpha_J candidate remains rejected | branch-residual | J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED | option C branch realization; pathA_21b P1b; pathA_22 |
| Theta_Q | Dimensionless branch factor relating the mouth/source flux to the far-field Gauss flux; computable only after R0, leakage, and kernel branch choices. | 1 | reduced-3D | pde.tex:511-539; pde.tex:2515-2566 | C_F factorization and far-field source normalization | branch-residual | THETA_Q_BRANCH_REALIZATION_UNDERIVED | option C branch realization; pathA_21b P1b/P4 |
| I_F,12 | Dimensionless Pi_cross control-surface integral I_F,12 defined by BVP-G3: normalize -int_{partial U2} Pi_cross.n dS by m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2). | 1 | reduced-3D | pde.tex:445-451; pde.tex:342-352; part01_parent_geometry.tex:275-286; pathA_21b BVP-G3 | C_F magnitude; attractive orientation only after profile sign integral | profile-solve | PI_CROSS_STRESS_TENSOR_UNDERIVED | option C; pathA_21b P1b |
| C_F,12 | G-free force coefficient C_F,12=(m_GNLS*N_infty,3*Q1*Q2/(4*pi))*I_F,12, with Qi=Theta_Qi*Ji/N_infty,3 from the Gauss-solved lane. | M L^3 T^-2 | reduced-3D | pathA_21b P1b Gauss solve plus BVP-G3 Pi_cross | P1b force form | profile-solve | ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL | pathA_21b P1b; pathA_22 only after mass bridge remains separately derived |
| alpha_H,omega | Would-be dimensionless profile energy ratio H_throat/(hbar*J_omega). pde.tex:318-391 is the action-level source, not a canonical Hamiltonian; the Hamiltonian and boundary mass relation must be constructed separately. | 1 | 4D-bulk / throat | pde.tex:318-391 action-level only; brane_bulk_ontology.tex:1267-1297 scaling only | alpha_J and m_defect bridge remain unclosed | new-physics | CANONICAL_THROAT_HAMILTONIAN_UNCONSTRUCTED / BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING | pathA_22 or later new-physics bridge |
| J_omega | Angular-rate version of the throat invariant used in alpha_J*hbar*J_omega/c_gamma^2; classification requires source relation to phase/angular rate. | T^-1 | throat / brane | pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60 | 2pi placement, alpha_J candidate | new-physics | H_2PI_RATE_CLASSIFICATION_UNDETERMINED | pathA_21 P2 |
| J_nu | Cycle-rate version of the throat invariant, with h*J_nu=2*pi*hbar*J_nu and alpha_J not absorbing 2*pi. | T^-1 | throat / brane | pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60 | 2pi placement, alpha_J candidate | new-physics | H_2PI_RATE_CLASSIFICATION_UNDETERMINED | pathA_21 P2 |
| M_inertial | Second derivative of the effective moving-throat action with respect to a slow center velocity after integrating fields over the solved throat/support domain. | M | brane / reduced throat | pde.tex:2806-2879 scope statement; no explicit accelerated-throat source equation in cited parents | EP inertial side | new-physics | INERTIAL_PROFILE_RESPONSE_MISSING | pathA_21 P2 EP |
| M_source | Mass normalization inferred from the far-field drain source after C_F factorization, separately from M_inertial and without using Newton G. | M | reduced-3D | pde.tex:396-451; brane_bulk_ontology.tex:1267-1297 | EP source side | new-physics | SOURCE_MASS_PROFILE_NORMALIZATION_MISSING | pathA_21 P2/P4 |
| C_B/C_E | Bulk transverse Maxwell principal coefficient ratio. The principal-symbol form is known, but the physical normalization is a calibration/branch datum rather than a throat BVP closure. | L^2 T^-2 | 4D-bulk | pde.tex:355-416; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:39-50 | bulk c_gamma normalization only after calibration/branch choice | branch-residual | BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED | pathA_22 calibration packet |
| lambda_gamma | Observed brane photon/sound ratio c_gamma/c_s. It requires the brane zero-mode/profile reduction and is not closed by P5b. | 1 | brane | pde.tex:541-565; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-51 | G conditional form and mass bridge c_gamma remain symbolic | new-physics | BRANE_ZERO_MODE_REDUCTION_UNDERIVED | option C zero-mode reduction; pathA_22 |
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
| c_s(branch) | Branch sound speed c_s^2=5K*rho0^4/m_GNLS evaluated on the BVP-G1 asymptotic/profile state. | L T^-1 | 4D-bulk / brane | pde.tex:342-352; software/stage1_solver/reports/pathA_20_velocity_constants.md:9-17 | profile sound speed for lambda_gamma and pathA_22 target | profile-solve | STATIONARY_GNLS_BVP_NOT_SOLVED | option C; pathA_22 |
| mhat | Source-map factor in the isotropic normalization law; in the 3D target lane its squared dimension must convert P0 to G*c_s^5/(a^5*c^5). | L^-1 T^-1 M^-1/2 in the 3D target normalization | PN-facing reduced observable | pde.tex:2077-2092; pde.tex:2551-2563 | pathA_22 scale map | pathA_22 | SCALE_MAP_SOURCE_FACTOR_UNDERIVED | pathA_22 |
| chi_Q | Outgoing-normalization scalar retained when the passive/outgoing branch is not fixed to canonical compact DtN. | 1 | PN-facing reduced observable | pde.tex:1980-1996; pde.tex:2053-2082; pde.tex:2551-2552 | pathA_22 outgoing normalization | pathA_22 | OUTGOING_DTN_BRANCH_UNDERIVED | pathA_22 |

## New Residuals

- `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`: P1b_SIGN_RESIDUAL. EOS/compressibility gives the pressure-drop sign, but the full attractive force sign is left to the solved profile integral. Source: pde.tex:445-451 plus the BVP-G3 control-surface integral; the parent does not select the finite-throat normal-orientation sign of I_F,12.
- `R0_FREE_BOUNDARY_CONDITION_UNDERIVED`: G2_BRANCH_REALIZATION_RESIDUAL. Option C must choose or derive the R0 branch selector. Source: part01_parent_geometry.tex:447-461 defines Sigma0 and candidate bottom data; no action variation with respect to R0(w) is supplied.
- `J_VALUE_BRANCH_PARAMETER`: G4_BRANCH_REALIZATION_RESIDUAL. J remains a branch parameter unless a choking, regularity, topology, or energy selector is added downstream. Source: pde.tex:396-406 and pde.tex:511-539 conserve/project flux but do not fix the value.
- `J_SELECTOR_UNDERIVED`: G4_BRANCH_REALIZATION_RESIDUAL. No silent free parameter; option C must enumerate or solve the selector. Source: No separate source-anchored regularity/choking/energy condition selecting J was found in the cited parents.
- `W_KERNEL_UNDERDECLARED`: G5_BRANCH_REALIZATION_RESIDUAL. W_eff and N_infty,3 remain branch-realization data. Source: pde.tex:496-565 and part01_parent_geometry.tex:298-390 define projection/reduction formulas but not a kernel-shape selector.
- `CANONICAL_THROAT_HAMILTONIAN_UNCONSTRUCTED`: ALPHA_H_ANCHOR_FIX. alpha_H,omega cannot be used as a mass bridge restatement. Source: pde.tex:318-391 is action-level; it is not a canonical Hamiltonian construction.

## Carried Negatives

Carried verbatim: `MASS_BRIDGE_FORM_NOT_DERIVED`, `EP_NOT_DERIVED`, `RETAIN_L_T_M`, `NEWTON_G_FORM_NOT_DERIVED`, `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`, `H_2PI_RATE_CLASSIFICATION_UNDETERMINED`, `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`, `HBAR_PROVENANCE_UNDETERMINED`, `HBAR_FREE_SUBSTRATE_RELATION_MISSING`, `W_EFF_REDUCTION_UNDERIVED`.

## Carried Residual Ledger

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

## Harness Summary

- Dimensional checks: 10 consistent, 0 inconsistent, 10 total.
- Algebraic checks: 5 consistent, 0 inconsistent, 5 total.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

