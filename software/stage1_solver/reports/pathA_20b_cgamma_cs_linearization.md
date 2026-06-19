# Path-A 20b c_gamma vs c_s coupled-linearization summary

## Verdicts

- `bulk_verdict`: `C_GAMMA_BULK_UNDERDETERMINED` with `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`.
- `brane_verdict`: `C_GAMMA_RATIO_STILL_UNDERDETERMINED` with `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`; sub-residual `BRANE_PHOTON_CONE_REQUIRES_PROFILE`.
- `lambda_gamma`: No pure number derived. If c_bulk is independent of rho0, lambda_gamma is proportional to rho0^-2; a pure number or equality would require the missing source equation.
- `pathA_21`: brane_verdict only; lambda_gamma remains symbolic.

No `EQUALS` verdict is claimed. The required source equation `C_B/C_E=5*K*rho0^4/m_GNLS` was not found in the cited parent-action sources.

## L1 background validity

- Status: `LEGAL_WITH_EXPLICIT_NEUTRALIZING_EXTERNAL_SOURCE`.
- `psi0`: psi0=sqrt(rho0)*exp(-i*mu*t/hbar) with uniform rho0 and v_b0=0.
- `A_M0`: A_M0=0, so F_MN0=0.
- `J_psi0`: J_psi0^0=q_star*rho0 and J_psi0^i=0.
- `J_ext0`: J_ext0^0=-q_star*rho0 and J_ext0^i=0.
- Neutrality: J_tot0^M=J_psi0^M+J_ext0^M=0, making the Maxwell background equation 0=0.
- Caveat: Without this explicit neutralizer the correct stop residual would be HOMOGENEOUS_CHARGE_NEUTRALITY_UNSPECIFIED.

Source anchor: pde.tex:370-374 permits explicit external/background sourcing; pde.tex:912-925 gives delta J^0=q_star*delta rho and the London current term.

## L1b coupled principal symbol

- Variables: `(delta rho, delta theta, delta A_M)`.
- Phonon block: `det P_ph=hbar*(omega^2 - (5*K*rho0^4/m_GNLS)*k^2)`.
- Gauge transverse block: `P_T=C_E*omega^2-C_B*k^2=C_E*(omega^2-c_bulk^2*k^2), c_bulk^2=C_B/C_E`.
- Coupled determinant: `det P = P_ph * P_T^2 for the physical transverse gauge polarizations after lower-derivative off-diagonal terms are dropped from the principal symbol`.
- Off-diagonal principal status: `VANISHES_ON_HOMOGENEOUS_NEUTRALIZED_BACKGROUND`.
- Physical photon cone: The cone is read from transverse field-strength modes. Gauge-fixing consistency is not used as speed evidence.

Lower-order/gapped terms, separated from the cone:
- delta J^0=q_star*delta rho
- delta J^i contains rho0*(hbar/m_GNLS)*grad(delta theta)
- delta J^i contains -(q_star/m_GNLS)*rho0*delta A^i, a London/plasma term
- background-gradient current terms vanish on the homogeneous background or are lower order

## L2 dispersions

- Phonon: omega^2=c_s0^2*k^2, c_s0^2=5*K*rho0^4/m_GNLS; quantum pressure gives the usual k^4 Bogoliubov correction but not c_gamma.
- Gauge: omega^2=c_bulk^2*k^2 plus possible lower-order London/plasma shifts; c_bulk^2=C_B/C_E.
- Massless/transverse branch status: BULK_PRINCIPAL_TRANSVERSE_BRANCH_ESTABLISHED; lower-order gapped/longitudinal branches are labeled separately and do not set the cone.

Machine checks confirm `[c_s]=[c_gamma]=L T^-1` and `C_B/C_E` has speed-squared dimension. These checks are non-evidentiary for equality.

## L3 two-layer verdict

- Bulk: `C_GAMMA_BULK_UNDERDETERMINED` because `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`. Conditional formula: `c_bulk/c_s0=sqrt((C_B/C_E)*m_GNLS/(5*K*rho0^4))`.
- Brane: `C_GAMMA_RATIO_STILL_UNDERDETERMINED` because `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`; profile sub-residual `BRANE_PHOTON_CONE_REQUIRES_PROFILE`.
- Rho dependence: No pure number derived. If c_bulk is independent of rho0, lambda_gamma is proportional to rho0^-2; a pure number or equality would require the missing source equation.

## L4 implications

- Standing-wave ceiling: The standing-wave ceiling remains c=c_gamma, not c_s unless the missing brane verdict later proves equality.
- Tail factor: R_tail=(c/c_s)^3=lambda_gamma^3; conditionally (c_bulk/c_s0)^3, not set to 1.
- Brane handoff: The brane localization/profile question is a blocking sub-residual for pathA_21, not an afterthought.

## Negative control

- Independent symbols: `c_bulk, c_s0`.
- Forbidden forced equality: `C_B/C_E=5*K*rho0^4/m_GNLS`.
- Result: `FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION`.

## Residuals carried

- `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`: BULK_VERDICT_RESIDUAL. bulk_verdict is C_GAMMA_BULK_UNDERDETERMINED. Carry c_bulk^2=C_B/C_E and c_bulk/c_s0=sqrt((C_B/C_E)*m_GNLS/(5*K*rho0^4)); do not set lambda_gamma=1. Source: part01_parent_geometry.tex:225-247 and pde.tex:357-416 give F_MN F^MN with no source-derived equation C_B/C_E=5*K*rho0^4/m_GNLS; pathA_20 showed Z(w) is an overall principal prefactor.
- `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING`: BLOCKS_BULK_EQUALS_C_S. An EQUALS verdict is forbidden unless a later source equation identifies the gauge kinetic metric with the acoustic metric. Source: part01_parent_geometry.tex:191-203 derives the acoustic speed from the EOS, while part01_parent_geometry.tex:225-247 gives the Maxwell parent metric separately; em_fields.tex:160-178,482-504,692-705 is legacy acoustic reuse, not parent-action evidence.
- `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`: BRANE_VERDICT_RESIDUAL. brane_verdict is C_GAMMA_RATIO_STILL_UNDERDETERMINED. pathA_21 consumes this brane verdict and keeps lambda_gamma symbolic. Source: part01_parent_geometry.tex:333-389 and pde.tex:541-565 give projection and zero-mode reduction as controlled assumptions; pde.tex:749-763 and 903-931 keep A_w, J^w, and F_muw alive in the microscopic linearized problem.
- `BRANE_PHOTON_CONE_REQUIRES_PROFILE`: BRANE_SUB_RESIDUAL. If the observed photon is a localized mixed-sector mode rather than a strict far-field zero mode, its cone must be computed from the solved profile and reduction kernel. Source: part01_parent_geometry.tex:944-946 and 1502-1511 state that the matched Maxwell/mixed reduction and actual nonlinear branch realization remain branch/profile tasks.

## Harness summary

- Dimensional checks: 11 consistent, 0 inconsistent, 11 total.
- Algebraic checks: 7 consistent, 0 inconsistent, 7 total.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

