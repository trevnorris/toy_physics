# Adversarial Audit Status

Generated: 2026-06-12T07:04:29Z
Project: pde_ledger

Authoritative consult record: `BATCHING_DECISIONS.md`.

| Scope | Counts |
|---|---|
| all candidates | audited=148 provenance_built=766 scanned=7 verdict_logged=1 |
| dry-run candidates | none |
| binding verdict fields populated | 0 |
| dedup canonicals / aliases | canonical=915 aliases=7 alias_groups=6 |

## Phase B Dedup And Families

- Dedup state: 915 canonical entries, 7 aliased entries.
- Family map: `redteam_adversarial/provenance/_family_map.yaml`
- Families: 563; singletons: 404; unmapped canonicals: 0.

Non-singleton family groupings:

| Family | Members | Parameters | Values | Stages |
|---|---:|---|---|---|
| fam_0005_a_f1 | 4 | A_F1 | - | 079, 080, 083, 089 |
| fam_0007_a_k | 2 | A_K | - | 054, 232 |
| fam_0008_a_t | 8 | A_T | -4.27263956256927, Sqrt[(9/20)(p/(1-sFormula/4))] | 147, 148, 151, 152 |
| fam_0012_a_w | 2 | A_w | - | 026, 111 |
| fam_0014_b_0 | 3 | B_0 | - | 011, 020, 023 |
| fam_0020_b_star | 2 | B_star | - | 184, 191 |
| fam_0028_c_mix | 6 | C_mix | - | 051, 081, 083, 084, 240, 242 |
| fam_0030_c_res | 3 | C_res | 0.994418836451529, 0.994418836451529348 | 067, 068, 069 |
| fam_0039_d2_over_d0 | 2 | D2_over_D0 | - | 170, 173 |
| fam_0040_d_0 | 2 | D_0 | 0 | 022, 025 |
| fam_0041_d_01 | 2 | D_01 | 0 | 173, 226 |
| fam_0051_delta | 6 | Delta | 84071/400 | 012, 178, 201, 204, 222, 247 |
| fam_0053_deltak_ax | 2 | DeltaK_ax | 0 | 029, 037 |
| fam_0055_delta_0 | 9 | Delta_0 | 0 | 058, 061, 066, 068, 069, 075, 083, 251 |
| fam_0058_delta_q | 10 | Delta_Q | 3 | 101, 106, 158, 192, 193, 194, 195, 196, 197 |
| fam_0061_delta_branch | 2 | Delta_branch | - | 191 |
| fam_0069_e_edge | 2 | E_edge | - | 250 |
| fam_0074_e_sub | 2 | E_sub | - | 248 |
| fam_0083_f_tr | 2 | F_tr | - | 045, 048 |
| fam_0087_g_w | 2 | G_W | - | 176, 247 |
| fam_0089_g_fail | 3 | G_fail | - | 061, 063 |
| fam_0093_gamma_5 | 10 | Gamma_5 | - | 021, 022, 023, 097, 099, 102, 104, 106, 251 |
| fam_0098_h_even | 2 | H_even | - | 226 |
| fam_0105_i_3 | 2 | I_3 | - | 192, 201 |
| fam_0109_i_f | 4 | I_f | 1/3 | 070, 071, 077 |
| fam_0112_j_1 | 2 | J_1 | - | 070 |
| fam_0115_k | 2 | K | - | 222 |
| fam_0123_k_sigma | 2 | K_Sigma | - | 019 |
| fam_0126_k_compat | 2 | K_compat | 24.473754879291 | 226, 228 |
| fam_0129_k_eta | 2 | K_eta | - | 001, 236 |
| fam_0132_k_geom | 3 | K_geom | 3, 3/4 | 091, 092 |
| fam_0135_k_m | 3 | K_m | - | 071, 073 |
| fam_0136_k_norm | 2 | K_norm | - | 012, 223 |
| fam_0137_k_q | 4 | K_q | - | 115, 118, 166 |
| fam_0138_k_req | 3 | K_req | 8 | 026 |
| fam_0139_k_s | 4 | K_s | - | 114, 115, 137, 138 |
| fam_0140_k_turn | 4 | K_turn | 2.73855812 | 243, 253 |
| fam_0141_kbar_0 | 3 | Kbar_0 | - | 001, 097 |
| fam_0145_l_w | 10 | L_W | - | 116, 117, 119, 121, 164, 165 |
| fam_0146_l_over_a | 13 | L_over_a | 20, 37/20, sp.sqrt(12***2/sp.pi**2-1) | 072, 073, 075, 078, 121, 131, 139, 146, 232 |
| fam_0147_l_r | 2 | L_r | - | 211, 214 |
| fam_0149_lambda_0 | 3 | Lambda_0 | 27/20 | 184, 189, 242 |
| fam_0156_m | 3 | M | - | 222, 228 |
| fam_0160_m_mix | 4 | M_mix | - | 036, 037, 043, 190 |
| fam_0161_m_q | 9 | M_q | - | 134, 135, 137, 139, 141, 142, 143, 145 |
| fam_0174_n_q | 7 | N_Q | 1, 54/5 | 023, 025, 030, 100, 195, 232 |
| fam_0183_omega_q | 4 | Omega_Q | - | 088, 105, 193, 195 |
| fam_0184_omega_u | 2 | Omega_U | - | 222 |
| fam_0185_omega_w | 2 | Omega_W | - | 222 |
| fam_0187_p0_target | 20 | P0_target | 38, 54, 54/5 | 019, 022, 025, 038, 177, 189, 193, 195, 197, 222, 223, 225 |
| fam_0188_p_2 | 7 | P_2 | 2 | 110, 192, 193, 194, 195, 196, 197 |
| fam_0194_p_res | 4 | P_res | 1.005612 | 068, 069 |
| fam_0198_pe_fail_chi | 2 | Pe_fail_chi | - | 089 |
| fam_0200_pe_req | 9 | Pe_req | 0 | 059, 061, 075, 078, 079, 082, 089, 090 |
| fam_0204_pi_can | 3 | Pi_can | 4.651033550168867 | 155, 156, 157 |
| fam_0209_pi_star | 14 | Pi_star | 1.50882951349316 | 034, 130, 131, 134, 135, 144, 150, 152, 155, 156, 180, 191, 198, 238 |
| fam_0210_pi_tr | 2 | Pi_tr | - | 240, 244 |
| fam_0214_r_u | 2 | R_U | - | 044, 045 |
| fam_0220_r_q | 5 | R_q | 1/4 | 138, 142, 144, 145 |
| fam_0222_r_target | 3 | R_target | - | 035, 047, 181 |
| fam_0223_r_tr | 2 | R_tr | 94 | 045 |
| fam_0228_s_can | 6 | S_can | 4.651033550168867, 4.651033550168876 | 155, 158, 163, 168 |
| fam_0229_s_q | 2 | S_q | - | 134, 150 |
| fam_0233_s_star | 2 | S_star | - | 152, 155 |
| fam_0237_sigma0_can | 6 | Sigma0_can | 4.651033550168867, 4.651033550168876 | 155, 156, 158, 163, 168 |
| fam_0238_sigma0_star | 2 | Sigma0_star | - | 155, 156 |
| fam_0239_sigma_0 | 3 | Sigma_0 | - | 135, 140, 154 |
| fam_0251_t_x | 2 | T_X | - | 070, 073 |
| fam_0252_t_can | 2 | T_can | 4.651033550168867, 4.651033550168876 | 155, 158 |
| fam_0263_theta_sigma | 2 | Theta_sigma | - | 140 |
| fam_0264_theta_w | 10 | Theta_w | - | 076, 077, 084, 166, 232 |
| fam_0266_u_obs | 2 | U_obs | 0.14313458 | 245 |
| fam_0268_upsilon_lat | 2 | Upsilon_lat | - | 253 |
| fam_0269_upsilon_lat_sess | 2 | Upsilon_lat_sess | - | 253 |
| fam_0274_v_known | 2 | V_known | 1.181909222592 | 222, 223 |
| fam_0280_w_wall | 2 | W_wall | - | 066, 070 |
| fam_0283_xi_1 | 9 | Xi_1 | - | 177, 225, 226, 227, 234, 243, 245 |
| fam_0285_xi_t | 2 | Xi_T | - | 099, 123 |
| fam_0287_xi_micro | 2 | Xi_micro | - | 060, 062 |
| fam_0291_y_q_cons | 4 | Y_Q_cons | - | 091, 092, 093, 095 |
| fam_0298_yhat_q_cons | 3 | Yhat_Q_cons | - | 097, 098, 099 |
| fam_0301_a | 2 | a | - | 021, 222 |
| fam_0302_a0 | 2 | a0 | - | 246, 247 |
| fam_0303_a_0 | 4 | a_0 | - | 109, 110, 111, 246 |
| fam_0309_alpha | 2 | alpha | - | 028, 029 |
| fam_0312_alpha_crit | 2 | alpha_crit | 0 | 028, 031 |
| fam_0317_assert_nonzero | 2 | assert_nonzero | - | 005, 012 |
| fam_0320_b_over_a | 2 | b_over_a | - | 024 |
| fam_0332_branch_sign | 2 | branch_sign | - | 125, 128 |
| fam_0333_c0 | 3 | c0 | - | 088, 096 |
| fam_0336_c_contact | 2 | c_contact | - | 088, 090 |
| fam_0337_c_geom | 3 | c_geom | 73 | 094, 096 |
| fam_0339_c_pole | 3 | c_pole | 1/4 | 093, 095, 098 |
| fam_0344_chi_0_star | 2 | chi_0_star | - | 186, 191 |
| fam_0345_chi_q | 30 | chi_Q | 1, 9 | 100, 101, 102, 103, 104, 105, 106, 107, 108, 110, 111, 112, 113, 192, 194, 195, 196, 197 |
| fam_0347_chi_peak | 2 | chi_peak | - | 250 |
| fam_0348_chi_raw | 2 | chi_raw | - | 250 |
| fam_0363_delta | 2 | delta | - | 038, 229 |
| fam_0380_ell | 2 | ell | - | 074 |
| fam_0381_engine_disagreement | 2 | engine_disagreement | - | 169, 252 |
| fam_0382_eps_2 | 3 | eps_2 | 0 | 093, 096 |
| fam_0385_eps_blk | 4 | eps_blk | 0, 1 | 081, 085, 086, 087 |
| fam_0389_eps_split_coefficient | 5 | eps_split_coefficient | - | 181, 182, 183, 184, 189 |
| fam_0390_epsilon | 3 | epsilon | 1 | 047, 241, 242 |
| fam_0391_epsilon_2 | 9 | epsilon_2 | - | 091, 092, 093, 094, 095, 096, 097, 098, 099 |
| fam_0393_epsilon_eta | 3 | epsilon_eta | - | 236, 237, 245 |
| fam_0399_f_lat | 4 | f_lat | - | 252, 253 |
| fam_0406_g_uv | 3 | g_UV | 0.14313458 | 245, 250 |
| fam_0407_g_frak | 2 | g_frak | - | 165, 167 |
| fam_0409_g_nat | 3 | g_nat | - | 122, 124 |
| fam_0412_g_star | 14 | g_star | 0.758035078944663, 4.651033550168876 | 122, 126, 128, 130, 152, 154, 155, 156, 163, 164, 165, 168, 169 |
| fam_0415_gamma_0 | 7 | gamma_0 | - | 114, 115, 116, 160, 161, 162 |
| fam_0417_gamma_w | 3 | gamma_W | 1/9 | 112, 159 |
| fam_0418_gamma_crit | 2 | gamma_crit | - | 251 |
| fam_0420_gamma_lattice_red | 2 | gamma_lattice_red | - | 253 |
| fam_0431_kappa | 7 | kappa | sp.sqrt(2)/((n+sp.Rational(1,2) | 024, 026, 070, 074, 222, 227 |
| fam_0433_kappa_0 | 3 | kappa_0 | - | 026, 032, 229 |
| fam_0434_kappa_0_sq | 2 | kappa_0_sq | - | 033 |
| fam_0444_lambda_20 | 2 | lambda_20 | 1 | 173, 224 |
| fam_0445_lambda_b | 2 | lambda_B | - | 222 |
| fam_0446_lambda_l | 5 | lambda_L | - | 247 |
| fam_0449_lambda_w | 2 | lambda_W | - | 222, 247 |
| fam_0452_lambda_mu | 3 | lambda_mu | 1 | 076, 086, 087 |
| fam_0459_matched_fingerprint_value | 28 | matched_fingerprint_value | - | 009, 021, 022, 024, 100, 101, 102, 103, 104, 105, 106, 108, 110, 111, 112, 115, 117, 158, 159, 160, 168, 172, 189, 194, 195, 196, 197, 220 |
| fam_0461_mathcal | 2 | mathcal | - | 186, 253 |
| fam_0462_mathfrak | 4 | mathfrak | - | 119, 123, 126, 144 |
| fam_0463_mathfrak_g | 2 | mathfrak_g | - | 119, 120 |
| fam_0469_mhat_0 | 5 | mhat_0 | 1 | 022, 032, 189, 195 |
| fam_0470_mhat_ang | 2 | mhat_ang | - | 024, 025 |
| fam_0474_moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit | 2 | moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit | - | 252, 253 |
| fam_0475_moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit | 2 | moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit | - | 252, 253 |
| fam_0476_mu | 2 | mu | - | 236, 239 |
| fam_0478_mu_eta | 6 | mu_eta | - | 001, 250, 252, 253 |
| fam_0488_paragraph | 2 | paragraph | - | 017, 021 |
| fam_0491_q_eta | 4 | q_eta | - | 192, 235, 237, 238 |
| fam_0493_r_f1 | 18 | r_F1 | 1.77799353547498, 1/4, 4.651033550168876, Sqrt[4107-100*Pi^2]/(10*Pi), sp.sqrt(4107-100*sp.pi**2)/(10*sp.pi), sp.sqrt(sp.Rational(12,1)/sp.pi**2 | 122, 123, 124, 128, 148, 154, 155, 156, 158, 161, 162, 163, 165, 168, 230, 246 |
| fam_0500_rho_alpha | 4 | rho_alpha | 4/3 | 085, 090, 098, 240 |
| fam_0502_s_0 | 2 | s_0 | - | 251, 253 |
| fam_0503_s_c | 3 | s_c | - | 251, 253 |
| fam_0506_sigma_pi | 2 | sigma_Pi | - | 129, 132 |
| fam_0517_stale_output | 5 | stale_output | - | 031, 055, 195, 200, 207 |
| fam_0522_t_cross | 2 | t_cross | - | 250, 251 |
| fam_0526_tautological_check | 2 | tautological_check | - | 062, 169 |
| fam_0532_u_2 | 2 | u_2 | 1/9 | 170, 172 |
| fam_0534_varpi | 2 | varpi | - | 222 |
| fam_0540_x | 2 | x | - | 055, 057 |
| fam_0543_x_star_exp | 2 | x_star_exp | - | 127 |
| fam_0544_x_star_slab | 2 | x_star_slab | - | 127 |
| fam_0549_xi_star | 2 | xi_star | - | 126 |
| fam_0552_zeta_max | 3 | zeta_max | 2.46752922945601, 71 | 082, 084, 098 |
| fam_0553_zeta_max_f1 | 2 | zeta_max_F1 | - | 098 |
| fam_0555_zeta_req | 4 | zeta_req | - | 048, 085, 240 |
| chain_aspect_37_20 | 44 | L_over_a, epsilon_r, eta, g_star, r_F1 | 0.758035078944663, 1.77799353547498, 1/20, 1/4, 20, 37/20, 4.651033550168876, Sqrt[4107-100*Pi^2]/(10*Pi), sp.sqrt(12***2/sp.pi**2-1), sp.sqrt(4107-100*sp.pi**2)/(10*sp.pi), sp.sqrt(sp.Rational(12,1)/sp.pi**2 | 072, 073, 075, 078, 121, 122, 123, 124, 126, 128, 130, 131, 139, 146, 148, 152, 154, 155, 156, 158, 161, 162, 163, 164, 165, 168, 169, 230, 232, 246 |
| chain_barrier_222_224 | 7 | DeltaV_req, T_quad, V_known, barP0_compat, epsilon_barrier, lambda_20, lambda_21, lambda_22 | 1, 1.181909222592 | 173, 222, 223, 224 |
| chain_calibration_245_253 | 22 | K_turn, Upsilon_lat, epsilon_eta, f_lat, gamma_lattice_legacy, gamma_lattice_red, lambda_L, m_s, mu_eta | 2.73855812 | 001, 236, 237, 243, 245, 247, 250, 252, 253 |
| chain_chi_Q_norm | 43 | N_Q, chi_Q, mhat_0, sigma_Q_can, xi_Q | 1, 54/5, 9 | 022, 023, 025, 030, 032, 100, 101, 102, 103, 104, 105, 106, 107, 108, 110, 111, 112, 113, 189, 192, 194, 195, 196, 197, 232 |
| chain_quad_54_5 | 42 | Gamma_5, N_Q, P0_target, mhat_0 | 1, 38, 54, 54/5 | 019, 021, 022, 023, 025, 030, 032, 038, 097, 099, 100, 102, 104, 106, 177, 189, 193, 195, 197, 222, 223, 225, 232, 251 |
| chain_sigma0_transport | 10 | Pi_can, S_can, Sigma0_can, T_can | 4.651033550168867, 4.651033550168876 | 155, 156, 157, 158, 163, 168 |
| chain_wall_action | 34 | K_Sigma, K_eta, P0_target, T_Omega, T_w, mhat_0, mu_eta | 1, 38, 54, 54/5 | 001, 019, 022, 025, 032, 038, 177, 189, 193, 195, 197, 222, 223, 225, 236, 250, 252, 253 |

## Phase C Batch Assignment

batch assignment pending (Step 5)
