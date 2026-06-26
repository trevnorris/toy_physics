BREATHING_CALIBRATED

## Operator, BCs, Inner Product

Modal equation: `mu_eta*q_tt - d_w(T_w*d_w q) + K_eta*q = S_0^(psi)+f_0^ext`.
`L0=mu_eta^(-1)*(-d_w(T_w*d_w .) + K_eta .)`, `K_eta=Tw*beta**2`.
Collective BCs: `{'alpha_a': 'alpha_a(0)=1, alpha_a(L0)=0', 'alpha_L': 'alpha_L(0)=0, alpha_L(L0)=rAL'}`.
`g`-lane BCs: `g(0)=0, T_w*g'(L0)=0`.
Inner product: `<f,g>_mu=4*pi*int_0^L0 mu_eta f g dw`.

## Profiles And Projection

`alpha_a=sinh(L0*beta - beta*w)/sinh(L0*beta)`.
`alpha_L=rAL*sinh(beta*w)/sinh(L0*beta)`.
Profile provenance: `derived`.
Residuals and BCs: `{'L0_alpha_a': '0', 'L0_alpha_L': '0', 'alpha_a_mouth': '1', 'alpha_a_cap': '0', 'alpha_L_mouth': '0', 'alpha_L_cap': 'rAL'}`.

M integrands:
- `aa`: `muEta*sinh(L0*beta - beta*w)**2/sinh(L0*beta)**2`
- `aL`: `muEta*rAL*sinh(beta*w)*sinh(L0*beta - beta*w)/sinh(L0*beta)**2`
- `LL`: `muEta*rAL**2*sinh(beta*w)**2/sinh(L0*beta)**2`

K integrands:
- `aa`: `Tw*beta**2*cosh(2*L0*beta - 2*beta*w)/sinh(L0*beta)**2`
- `aL`: `-Tw*beta**2*rAL*cosh(L0*beta - 2*beta*w)/sinh(L0*beta)**2`
- `LL`: `Tw*beta**2*rAL**2*cosh(2*beta*w)/sinh(L0*beta)**2`

`M_AB={'aa': '-2*pi*muEta*(L0*beta*tanh(L0*beta) - sinh(L0*beta)**2)/(beta*sinh(L0*beta)**2*tanh(L0*beta))', 'aL': '2*pi*muEta*rAL*(L0*beta - tanh(L0*beta))/(beta*sinh(L0*beta)*tanh(L0*beta))', 'LL': '-2*pi*muEta*rAL**2*(L0*beta*tanh(L0*beta) - sinh(L0*beta)**2)/(beta*sinh(L0*beta)**2*tanh(L0*beta))'}`, `det(M)=-4*pi**2*muEta**2*rAL**2*(L0*beta - sinh(L0*beta))*(L0*beta + sinh(L0*beta))/(beta**2*sinh(L0*beta)**2)`.
`K_AB={'aa': '4*pi*Tw*beta/tanh(L0*beta)', 'aL': '-4*pi*Tw*beta*rAL/sinh(L0*beta)', 'LL': '4*pi*Tw*beta*rAL**2/tanh(L0*beta)'}`, `det(K)=16*pi**2*Tw**2*beta**2*rAL**2`.

## Dynamical EOM

`M_AB*Qddot^B + K_AB*Q^B = F_A^(HF)` with `Q=(delta_a, delta_L)`.
Expanded rows: `['M_aa*delta_a_ddot + M_aL*delta_L_ddot + K_aa*delta_a + K_aL*delta_L = F_a', 'M_aL*delta_a_ddot + M_LL*delta_L_ddot + K_aL*delta_a + K_LL*delta_L = F_L']`.
`K_AB` provenance: `operator_projection_not_static_Hessian`.

## Truncation Consistency

Generalized Galerkin basis: `B={alpha_a, alpha_L, g_1..g_N}; generalized problem K v = omega^2 M v`.
At `beta_L0=1.85`, `N=16`: `o_1=0.993109102589`, `o_2=0.98776369936`, `min(omega_1^2,omega_2^2)=3.42251944599`, `gap=2.22787035351`.
`beta_from_R0={'status': 'derived_from_frozen_constant_coeff_packet', 'R0_reference': 'Gate-1 frozen finite throat: L0=37/20, Tw=1, K_eta=1', 'beta': 1.0, 'L0': 1.85, 'beta_L0': 1.85, 'provenance': 'calibrated wall constitutive packet; R0 geometry alone does not derive Tw or K_eta'}`.
`beta_window={'beta_L0_min_in_sweep': 0.1, 'beta_L0_max_in_sweep': 3.0, 'predicate': 'o_1,o_2 >= 0.9 and min(omega_1^2,omega_2^2)>0'}`.

Full beta sweep:
- beta_L0=0.1, o_1=0.999999930575, o_2=0.992719635542, min_omega12_sq=0.0100000001662, pass=True
- beta_L0=0.2, o_1=0.999998890781, o_2=0.992680172759, min_omega12_sq=0.0400000026607, pass=True
- beta_L0=0.5, o_1=0.999957103167, o_2=0.992401710672, min_omega12_sq=0.250000103913, pass=True
- beta_L0=1, o_1=0.999337943286, o_2=0.99137629996, min_omega12_sq=1.00000166202, pass=True
- beta_L0=1.85, o_1=0.993109102589, o_2=0.98776369936, min_omega12_sq=3.42251944599, pass=True
- beta_L0=2, o_1=0.990850612592, o_2=0.986841633425, min_omega12_sq=4.00002655491, pass=True
- beta_L0=3, o_1=0.963700535077, o_2=0.978084610117, min_omega12_sq=9.00013411608, pass=True
- beta_L0=5, o_1=0.859847180331, o_2=0.945331526733, min_omega12_sq=25.001026871, pass=False
- beta_L0=10, o_1=0.631305307872, o_2=0.813703842156, min_omega12_sq=100.015797523, pass=False
- beta_L0=18, o_1=0.46589994191, o_2=0.641585934163, min_omega12_sq=324.14391923, pass=False
- beta_L0=30, o_1=0.347789348428, o_2=0.483518176405, min_omega12_sq=900.70444825, pass=False
- beta_L0=50, o_1=0.266889377513, o_2=0.318493859849, min_omega12_sq=2501.66782695, pass=False

Convergence in N: `[{'beta_L0': 1.85, 'N': 4, 'profile': 'alpha_a=sinh(beta*(L0-w))/sinh(beta*L0)', 'o_1': 0.9930092905690113, 'o_2': 0.9871004922686136, 'omega1_squared': 3.4236392243245994, 'omega2_squared': 13.327931790649671, 'omega3_squared': 43.329698015773054, 'min_omega12_squared': 3.4236392243245994, 'gap': 2.2510444003150876, 'pass': True, 'rank_deficient_basis': False, 'mass_condition': 3143.999988790218}, {'beta_L0': 1.85, 'N': 8, 'profile': 'alpha_a=sinh(beta*(L0-w))/sinh(beta*L0)', 'o_1': 0.993097805037466, 'o_2': 0.9876694141071264, 'omega1_squared': 3.422652656604015, 'omega2_squared': 13.296755443232552, 'omega3_squared': 42.950858398158616, 'min_omega12_squared': 3.422652656604015, 'gap': 2.2301758561723917, 'pass': True, 'rank_deficient_basis': False, 'mass_condition': 23898.670290241742}, {'beta_L0': 1.85, 'N': 12, 'profile': 'alpha_a=sinh(beta*(L0-w))/sinh(beta*L0)', 'o_1': 0.9931068863770757, 'o_2': 0.9877444260283081, 'omega1_squared': 3.4225458658290386, 'omega2_squared': 13.29349403082082, 'omega3_squared': 42.91558377486189, 'min_omega12_squared': 3.4225458658290386, 'gap': 2.228314818915371, 'pass': True, 'rank_deficient_basis': False, 'mass_condition': 79940.76565409845}, {'beta_L0': 1.85, 'N': 16, 'profile': 'alpha_a=sinh(beta*(L0-w))/sinh(beta*L0)', 'o_1': 0.9931091025894195, 'o_2': 0.9877636993599208, 'omega1_squared': 3.4225194459861923, 'omega2_squared': 13.292692442596563, 'omega3_squared': 42.90708785376419, 'min_omega12_squared': 3.4225194459861923, 'gap': 2.2278703535085196, 'pass': True, 'rank_deficient_basis': False, 'mass_condition': 188973.08392677506}]`.

## Structure

Legacy `E_geom`: `1/2*kappa*(delta_L-chi*delta_a)^2 + 1/2*sigmaA*delta_a^2 + 1/2*sigmaL*delta_L^2`.
Full Hessian: `{'aa': 'chi**2*kappa + sigmaA', 'aL': '-chi*kappa', 'LL': 'kappa + sigmaL'}`.
Pattern check: `symmetric full-rank ratio-plus-support Hessian with negative off-diagonal for chi>0`, M_posdef=`True`, K_structure_ok=`True`, K_offdiag_negative=`True`.
Computed structure provenance: `structure_from_computed_matrices=True`.
Able-to-fail probe: `{'non_posdef_M_probe': {'mutation': 'M_aa -> -M_aa', 'M_posdef': False, 'leading_minor_positive': False, 'det_positive': False}, 'sign_flipped_K_probe': {'mutation': 'K_aL -> -K_aL', 'K_offdiag_negative': False, 'K_structure_ok': False}}`.

## Hellmann-Feynman Force

Measure: `action measure 4*pi*int dw; not mu_eta-weighted`.
Distributed unsimplified: `{'F_a': '4*pi*Integral(Vp0*rhoStar*sinh(L0*beta - beta*w)/(ellC*sinh(L0*beta)), (w, 0, L0))', 'F_L': '4*pi*Integral(Vp0*rAL*rhoStar*sinh(beta*w)/(ellC*sinh(L0*beta)), (w, 0, L0))'}`.
Legacy unsimplified: `{'F_a': '-4*pi*Integral(-Vp0*rhoStar*sinh(L0*beta - beta*w)/(ellC*sinh(L0*beta)), (w, 0, L0))', 'F_L': '-4*pi*Integral(-Vp0*rAL*rhoStar*sinh(beta*w)/(ellC*sinh(L0*beta)), (w, 0, L0))'}`.
Simplified differences: `{'F_a': '0', 'F_L': '0'}`.

## Static-Dynamic Limit

`Qddot -> 0 in the same M_AB Qddot + K_AB Q = F_A system gives K_AB Q = F_A; after the legacy pattern comparison this is the static dE_geom/dQ=0 limit.`.

## Counterfactual Guard

`{'calibrations_frozen': True, 'alpha_L_frozen': True, 'refit_forbidden': True, 'degenerate': {'profile': 'alpha_a=0', 'M_det': '0', 'M_posdef': False, 'wrong_o_k': {'o_1': 0.9690040196623307, 'o_2': 0.2226896627819251, 'rank_deficient_basis': True}, 'F_a_dist': '0', 'F_a_legacy_frozen': '4*pi*Vp0*rhoStar*(cosh(L0*beta) - 1)/(beta*ellC*sinh(L0*beta))', 'hf_mismatch': True, 'fails': True}, 'nontrivial': {'profile': 'alpha_a=1', 'M_det': '-8*pi**2*muEta**2*rAL**2*(L0**2*beta**2*tanh(L0*beta) - L0*beta*sinh(L0*beta)**2 + 2*sinh(L0*beta)**2*tanh(L0*beta) - 4*cosh(L0*beta)*tanh(L0*beta) + 4*tanh(L0*beta))/(beta**2*sinh(L0*beta)**2*tanh(L0*beta))', 'M_posdef': True, 'wrong_o_k': {'o_1': 0.9999999999999994, 'o_2': 0.9738471876732767}, 'F_a_dist': '4*pi*L0*Vp0*rhoStar/ellC', 'F_a_legacy_frozen': '4*pi*Vp0*rhoStar*(cosh(L0*beta) - 1)/(beta*ellC*sinh(L0*beta))', 'hf_mismatch': True, 'overlap_passes': True, 'fails': True}}`.

## Engine Agreement

`{'status': 'pass', 'engine_agreement': True, 'digest_matches': True, 'mathematica_yaml': 'software/stage1_solver/_scratch/pathA_31_mathematica_results.yaml', 'mathematica_route': 'native_DSolveValue_BVP_Integrate_plus_NIntegrate_generalized_Eigensystem', 'checks': {'alpha_a': True, 'alpha_L': True, 'M_aa': True, 'M_aL': True, 'M_LL': True, 'K_aa': True, 'K_aL': True, 'K_LL': True, 'M_det': True, 'K_det': True, 'force_dist_a': True, 'force_dist_L': True, 'force_legacy_a': True, 'force_legacy_L': True, 'wrong_zero_M_det': True, 'wrong_const_M_det': True, 'wrong_zero_F_a': True, 'wrong_const_F_a': True, 'legacy_H_aa': True, 'legacy_H_aL': True, 'legacy_H_LL': True, 'structure_M_posdef': True, 'structure_K_symmetric': True, 'structure_K_offdiag_negative': True, 'structure_K_structure_ok': True, 'structure_K_rank': True, 'structure_probe_M_posdef': True, 'structure_probe_K_structure_ok': True}, 'numeric_checks': {'beta_sweep': True, 'selected_overlap': True, 'wrong_zero_overlap': True, 'wrong_const_overlap': True}, 'max_numeric_abs_delta': 5.684341886080801e-13}`.

## Dimensional Check

`{'status': 'pass', 'dimension_order': '{M,L,T}', 'dimensions': {'Q=(delta_a,delta_L)': '(0, 1, 0)', 'muEta': '(1, -1, 0)', 'Tw': '(1, 1, -2)', 'K_eta': '(1, -1, -2)', 'beta': '(0, -1, 0)', 'source_density': '(1, 0, -2)'}, 'checks': {'M_AB_dimension_mass': True, 'kinetic_energy_dimension': True, 'K_AB_from_Tw_dimension': True, 'K_AB_from_Keta_dimension': True, 'potential_energy_dimension': True, 'HF_force_dimension': True, 'perturbed_Keta_dimension_fails': True, 'perturbed_alpha_L_dimension_fails': True}}`.

## Reduction Certificate

`{'ell0_restriction': {'Y00': '1/(2*sqrt(pi)); int_S2 Y00^2 dOmega = 1', 'eta': 'eta(w,t)=eta_00(w,t)*Y00', 'T_Omega_term': 'drops because l(l+1)=0'}, 'background': {'gate1_reference': 'straight constant-coefficient finite throat, rho0=rho_star', 'R0_endpoint_data': 'R0(0)=a0, R0(L0)=0; frozen Gate-1 packet has L0=37/20 and Tw=K_eta=1', 'domain': '0 <= w <= L0', 'rho0': 'rhoStar'}, 'source_linearization': {'Sigma': 'r-R0(w)-eta', 'delta_V_conf': '-(V_wall_prime_0/ell_c)*eta', 'action_force_density': 'rhoStar*V_wall_prime_0/ellC, conjugate to eta', 'action_measure': 'int f_Sigma alpha_A sqrt(gamma0) dw dOmega; straight-reference normalization gives 4*pi int dw, not the mu_eta inner product'}, 'input_provenance': {'alpha_a': 'derived: harmonic lifting L0 alpha_a=0 with alpha_a(0)=1, alpha_a(L0)=0 from Gate-1 endpoint data', 'alpha_L': 'derived: harmonic lifting L0 alpha_L=0 with alpha_L(0)=0, alpha_L(L0)=rAL from straight-throat length endpoint data', 'muEta': 'calibration', 'Tw': 'calibration', 'K_eta': 'calibration_tied_to_beta_squared_Tw', 'beta_from_R0': 'beta=sqrt(K_eta/Tw)=1 from frozen constant-coefficient packet; geometry alone does not derive it'}, 'phonon_limit_caveat': 'No BdG k^4 matter-sector term enters this gate; Gate-1 caveat remains deferred under k*xi << 1 if matter dynamics are restored.', 'forbidden_fit_flags': {'K_from_static_hessian': False, 'M_or_K_typed_from_legacy_values': False, 'full_matrix_fit': False, 'chain_rule_HF_identity_used': False, 'counterfactual_flags_hardcoded': False}}`.

## Files

- `sympy_engine`: `software/stage1_solver/tools/pathA_31_scalar_breathing_sympy.py`
- `mathematica_engine`: `software/stage1_solver/tools/pathA_31_scalar_breathing.wl`
- `report`: `software/stage1_solver/reports/pathA_31_scalar_breathing.md`
- `results_yaml`: `software/stage1_solver/reports/pathA_31_results.yaml`
- `feed_note`: `research/pde_ledger/notes/stages/moving_throat_pde_pathA_31_scalar_breathing_result.md`
- `sympy_scratch_yaml`: `software/stage1_solver/_scratch/pathA_31_sympy_results.yaml`
- `mathematica_scratch_yaml`: `software/stage1_solver/_scratch/pathA_31_mathematica_results.yaml`
- `sympy_expr_export`: `software/stage1_solver/_scratch/pathA_31_sympy_exprs.wl`
