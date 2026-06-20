# Path-A 21c force from Noether stress tensor

## Verdicts

- Corrected P1b -> P1c label: `P1c: G_FREE_NOETHER_STRESS_STRUCTURE_POWER_DERIVED_WITH_FAR_FIELD_MATTER_ATTRACTIVE_SIGN_AND_SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`.
- Leading matter-stress sign: `FORCE_ATTRACTIVE_DERIVED`.
- Full sign verdict: `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`.
- Normalization status: CALIBRATION KNOB, not derived.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

Dual-engine scope: Python and Mathematica check dimensions and algebra only. The derivation is the Noether balance construction plus the explicit traction integrals below; exit 0 is necessary, not sufficient.

## Noether Balance

Noether trace:
- Start from L_psi=(i*hbar/2)(psi^*D_t psi-psi D_t psi^*)-(hbar^2/(2*m_GNLS))(D_i psi)^*(D_i psi)-V_conf*rho-U(rho).
- For an active spatial translation delta psi=-epsilon_i partial_i psi and delta A_M=-epsilon_i partial_i A_M, the canonical identity is partial_t T^0_i+partial_j T^j_i=partial L/partial x_i for explicit backgrounds.
- The matter canonical flux from the spatial-gradient term is reduced on shell with psi=sqrt(rho) exp(i theta), j_i=rho v_i, the phase equation, and P=rho*h-U to Pi_ij^matter=m_GNLS*rho*v_i*v_j+delta_ij*P+sigma_Q,ij.
- The explicit matter background term -V_conf(X;Sigma0)*rho contributes partial L/partial x_i=-rho*partial_i V_conf, so it is f_i^body, not part of Pi_ij.
- The Maxwell action contributes the standard field stress only in lanes where Z(w), gauge fixing, and J_ext backgrounds are proven to vanish/cancel; otherwise their explicit partial L/partial x_i terms are residualized.

Momentum density: `g_i=m_GNLS*rho*v_i`.

Matter stress representative: `Pi_ij^matter=m_GNLS*rho*v_i*v_j+delta_ij*P(rho)+sigma_Q,ij`.

Quantum stress: `sigma_Q,ij=(hbar^2/(4*m_GNLS))*[(partial_i rho partial_j rho)/rho-partial_i partial_j rho]`.

EOS pressure identity: `P=K*rho^5, h=(5K/4)*rho^4, dP=rho*dh`.

Matter balance law: `partial_t g_i+partial_j Pi_ij^matter=q_star*rho*(E_i+v_j B_ij)-rho*partial_i V_conf`.

Explicit-background body forces:
- `f_i^Vconf=-rho*partial_i V_conf`.
- `f_i^Z=-(partial_i Z) F_MN F^MN/(4*mu0) in the Maxwell sector`.
- `f_i^Jext=-A_M*partial_i J_ext^M for explicit external-source gradients`.
- `gauge-fixing/background terms are included only after a selected branch proves cancellation; otherwise residualized`.

Balance-law-vs-Euler check: Using continuity, dP=rho*dh, and partial_j sigma_Q,ij=rho*partial_i Q, the balance law divided by rho reproduces m_GNLS*(partial_t+v.grad)v_i=q_star*(E_i+v_jB_ij)-partial_i(V_conf+h+Q).

Stress representative: Canonical Noether stress is reduced to the displayed Madelung hydrodynamic representative. No Belinfante improvement is used in the accepted matter-stress integral; smooth closed-surface improvements integrate to zero by Gauss/antisymmetry, while singular core/profile pieces are not used as derived normalization.

## Force Integral

Convention: n_hat outward from U2; F_12 is force on defect 2 by defect 1; stationary F_12=-int_boundary Pi_ij n_j dS+int_U2 f_i^body dV.

Control surface: reduced-3D sphere around defect 2 with core scale << a << r12; v2=-Q2*n_hat/(4*pi*a^2), v1 is constant over the sphere to leading order.

Reduced-3D compact lane integral results:

- Convective cross flux: `int Pi_conv,cross.n dS=-(4/3)*m_GNLS*N_infty,3*Q2*v1`.
- Bernoulli/EOS pressure cross flux: `delta P_cross=-m_GNLS*N_infty,3*(v1.v2), so int delta_ij P_cross n_j dS=+(1/3)*m_GNLS*N_infty,3*Q2*v1`.
- Matter cross flux: `int Pi_matter,cross.n dS=-m_GNLS*N_infty,3*Q2*v1`.
- Force result: `F_12^matter=m_GNLS*N_infty,3*Q2*v1=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*rhat_12`.

Bulk lane kept separate:

- `F_12^matter,4D=-(m_GNLS*rho_infty,4*Q1*Q2/(2*pi^2*R12^3))*Rhat_12`.

Structure result: Bilinear Q1*Q2 structure is an integral result from the v1*v2 cross stress plus Bernoulli pressure cross term.

Power result: Reduced compact lane gives r12^-2 because v1 from the carried 4*pi Gauss solve is r12^-2; unreduced bulk gives R12^-3 with Omega_3=2*pi^2.

Normalization: Overall I_F,12^full / Theta_Q / branch-profile normalization remains a CALIBRATION KNOB, not derived.

## Sign

Far-field matter-stress verdict: `FORCE_ATTRACTIVE_DERIVED`.

Full sign verdict: `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`.

Orientation: With rhat_12 from defect 1 to defect 2, the matter-stress force is proportional to -Q1*Q2*rhat_12. Like drains/sources attract; opposite signs repel.

Residual reason: The matter-stress sign is derived target-blind, but pathA_21c does not evaluate the quantum, V_conf body-force, and Maxwell profile pieces that could enter the full normalized control-volume force.

## Calibrate-Predict Ledger

Target-blind predictions:
- force structure: bilinear Q1*Q2 from the stress integral.
- lane power: reduced-3D r^-2 and separate bulk R^-3 from Gauss plus surface measure.
- leading matter-stress sign: attractive for like drains, repulsive for opposite signs.

Calibration knobs:
- overall dimensionless normalization class: I_F,12^full / Theta_Q / branch-profile data.

Counts: predictions=3; knobs=1.

Guardrail: No observable is both calibrated-to and predicted; knobs < independent target-blind predictions. The full sign remains residual rather than being hidden in the normalization knob.

## Residuals

- `VCONF_BODY_FORCE_RESIDUAL`: BODY_FORCE_PROFILE_RESIDUAL. The exterior hydrodynamic matter-stress integral is derived, but the full core/body-force normalization is not called derived. Source: The action contains -V_conf(X;Sigma0)*rho, so f_i^Vconf=-rho*partial_i V_conf. The selected finite throat profile is not solved in pathA_21c.
- `QUANTUM_STRESS_FAR_FIELD_RESIDUAL`: PROFILE_DERIVATIVE_RESIDUAL. Quantum stress is not used to tune or flip the derived matter-stress sign. Source: sigma_Q,ij is derived and its divergence check is machine-verified, but its cross surface integral needs the density derivative profile near the throat branch.
- `MAXWELL_Z_GAUGE_JEXT_CANCELLATION_RESIDUAL`: MAXWELL_BODY_FORCE_RESIDUAL. Maxwell stress is not included in the force coefficient until the profile/zero-mode branch proves the cancellation. Source: The localized Maxwell action has explicit Z(w), gauge fixing, and -A_M J_ext^M terms. The compact matter-drain lane does not prove their cross terms vanish or cancel.
- `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`: FULL_SIGN_RESIDUAL. The full sign is not claimed as derived; the leading compact reduced-3D matter-stress sign is a target-blind far-field result. Source: The convective plus Bernoulli-pressure matter stress gives an attractive like-drain far-field sign, but quantum, V_conf body-force, and Maxwell profile pieces are not all evaluated.

## Carried Items

- pathA_21b G1 stationary BVP and BC table.
- reduced-3D Gauss solve v_r=-Theta_Q J/(4*pi*N_infty,3*r^2).
- bulk-4D Gauss solve v_R=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3).
- pathA_21b I_F,12^full definition carried as calibration/profile knob.
- MASS_BRIDGE_FORM_NOT_DERIVED.
- EP_NOT_DERIVED.
- RETAIN_L_T_M.
- NEWTON_G_FORM_NOT_DERIVED.
- pathA_20/20b residuals carried unchanged.

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

- Dimensional checks: 9 consistent, 0 inconsistent, 9 total.
- Algebraic checks: 9 consistent, 0 inconsistent, 9 total.

