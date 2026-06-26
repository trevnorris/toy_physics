# PathA-34 cross-ell unification result

Computed verdict: `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.

The dual-engine headline derives the outgoing DtN fingerprints for ell=0, ell=1, and ell=2, assembles the scalar/dipole residuals, and reruns the Gate-4 quadrupole non-regression check. The baseline does not pass: the native fixed-P0 nullspace moves `T0(0)`/`T1(0)` and no derived Gate-5 selector removes that freedom.

## Sector provenance

- ell=0 map: `handoff section 8.2 / Gate-2 collective (delta_a,delta_L) reduction`.
- ell=1 map: `handoff section 9.4 harmonic reduction; angular stiffness K_eta+2*T_Omega`.
- ell=2 map: `handoff sections 9.4 and 10.2-10.3 grouped-P2 port kernel`.
- Native nullspace dimension: `8`; return-moving nullity: `2`; moves return: `True`; selector present: `False`.
- Named constraints collected: `ell0_collective_gate2_stiffness, ell1_section9_4_harmonic_stiffness, ell2_section10_port_kernel`.
- Return dofs untouched by named constraints: `['Z0_ret', 'Z1_ret']`; pathA_29 premise: `Z_is_premise=True, boundary_dof=none`.
- Selector equation control: `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` -> `CROSS_L_RESIDUAL_PREDICTION` using `['Z0_ret = K0c', 'Z1_ret = K_eta + 2*T_Omega']`.

## Residual class

- `A0`: `I*M0*Z0_ret*a*omega/(c_s*(K0c + Z0_ret))`.
- `A1`: `I*D1*Z1_ret*a**3*omega**3/(2*c_s**3*(K_eta + 2*T_Omega + Z1_ret))`.
- `epsilon0_eff`: `Z0_ret/K0c`; `epsilon1_eff`: `Z1_ret/(K_eta + 2*T_Omega)`.

## Dimensional check

- Dimensional verdict: `NO_FAIL`; `dimensional_ok=True`.
- Sourced corruption verdict: `FAIL_DIMENSIONAL`.
- Unconstrained carrier corruption verdict: `NO_FAIL`.

## Dual-engine agreement

- Status: `pass`; max symbolic delta `0.0`; max numeric delta `5.551115123125783e-17`.

## Probe verdicts

- R1 port-kernel dependency: P0 changes `True`; rank row changes `True`.
- 3a decouple knobs: `FAIL_DECOUPLED` / `CROSS_L_RESIDUAL_PREDICTION`.
- 3b injected null detector: `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` / `CROSS_L_RESIDUAL_PREDICTION`; baseline native: `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.
- 3c wrong sign: `FAIL_EPSILON_MISMATCH` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.
- 3d perfect return: `FAIL_OVERCANCEL` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.
- 3e break Gate 4: `FAIL_QUAD_REGRESSION` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.
- 3f corrupt dimension: `FAIL_DIMENSIONAL` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`; free carrier `NO_FAIL`.
- 3g assert not derive: `FAIL_TAUTOLOGICAL` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.
- 3h no consistent return: `FAIL_NO_CONSISTENT_RETURN` / `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`.

## Earned vs deferred

Earned: the `-(ell+1)/Lambda_l^out` fingerprints, raw radiative orders, residual form/sign/order conditional on a positive bounded branch, and the Gate-4 non-regression. Deferred: the scalar/dipole return magnitude and nonzero prediction, because the native nullspace leaves `epsilon_eff` free at the linear Gate-5 level.
