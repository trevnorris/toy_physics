# Numbering Reconciliation Scan -- Moving-Throat PDE (stages 001-253)

**Date:** <FILL: scan date>  
**Status:** READ-ONLY discovery enumeration. No files edited. Canonical numbering = scripts/ filenames + paper/stages/*.tex + redteam/MANIFEST.yaml.

## Method

The derivation has 253 canonical stages. During an EM-extension realignment, content blocks were renumbered and a global stage-number bump missed stragglers, leaving STALE stage-number labels in notes, .wl/.py banner/Print strings, and (notably) some published paper-card titles.

**Canonical numbering is GROUND TRUTH**, defined by (priority order): scripts/moving_throat_pde_stageNNN_<stem>_sympy_audit.py -> mathematica/...stageNNN_<stem>_mathematica_audit.wl -> paper/stages/stage_NNN.tex (label `\label{stage:NNN}`) -> redteam/MANIFEST.yaml.

Script **stems are invariant** across the realignment, so reconciliation is **STEM-KEYED**. The canonical stem->NNN map was built by scanning every script/.wl filename. **No stem maps to more than one canonical number** (verified -- zero ambiguous stems), which makes stem-keying safe. **NEVER apply a uniform per-file or per-batch offset** -- observed offsets are inconsistent (+17 dominant but also +1,+2,...,+22,+40,+51,+68,+215,-31,-67), so a blind sweep would mint confidently-wrong references.

### Coverage / map provenance
- 242 of 253 stages have a script-derived stem. The remaining **11 are tex-only** (status/checkpoint stages, no script): 103, 113, 120, 124, 128, 132, 136, 141, 145, 149, 153 -- for these the filename number IS canonical.
- Map is regenerable from filenames; a representative slice is below.

### Compact stem->canonical map (representative slice; full map = the filename set)

| canonical NNN | stem |
| --- | --- |
| 001 | `geometry_lift` |
| 002 | `breathing_reduction` |
| 003 | `bdg` |
| 004 | `projected_maxwell_bundle_index` |
| 005 | `projected_maxwell_covariant` |
| 006 | `projected_maxwell_vector` |
| 007 | `projection_reduction_comparison` |
| 008 | `projected_maxwell_extension` |
| 009 | `projected_maxwell_near_throat` |
| 010 | `projected_maxwell_push_bundle_master` |
| 089 | `family1_minimal_isotropic_verdict` |
| 090 | `updated_reduced_status` |
| 091 | `grouped_p2_static_geometry_derivation` |
| 092 | `dynamic_geometry_obstruction` |
| 093 | `grouped_p2_status_update` |
| 094 | `isotropic_geometry_decoupling` |
| 095 | `second_order_geometry_contamination` |
| 096 | `geometry_lane_check_verdict` |
| 097 | `single_normalization_defect` |
| 098 | `family1_support_is_automatic` |
| 099 | `reduced_finish_line` |
| 100 | `outgoing_normalization_factorization` |
| 101 | `natural_source_map_reduction` |
| 102 | `higher_odd_irrelevance` |
| 104 | `outgoing_dtn_fingerprint` |
| 105 | `chiQ_fix_from_outgoing_dtn` |
| 106 | `canonical_outgoing_reduced_closure` |
| 107 | `general_dtn_deformation` |
| 108 | `robustness_classes` |
| 109 | `linearized_branch_selection` |
| 110 | `robin_outlet_model` |
| 111 | `mixed_sidechannel_pole` |
| 112 | `hybrid_robin_mixed_compensation` |
| 114 | `concrete_core_schur` |
| 115 | `core_balance_compensation` |
| 116 | `dn_mixed_tube_realization` |
| 117 | `outlet_core_status` |
| 118 | `parent_core` |
| 119 | `parent_balance` |
| 121 | `geometric_r_selection` |
| 122 | `mouth_source_compensation_test` |
| 123 | `parent_normalized_branch_values` |
| 125 | `positive_source_theorem` |
| 168 | `off_bundle_slippage` |
| 169 | `no_linear_p2_scalar_slippage` |
| 170 | `linear_grouped_outlet_map` |
| 171 | `microscopic_grouped_obstructions` |
| 172 | `physical_slope_collapse` |
| 173 | `axisymmetric_loading_mismatch` |
| 174 | `static_self_similarity` |
| 175 | `wall_normalized_load_shape` |
| 176 | `outgoing_load_factorization` |
| 177 | `weak_axisymmetric_outgoing_slippage` |
| 178 | `outgoing_port_coloading` |
| 179 | `transfer_shape` |
| 180 | `effective_transfer_shape_collapse` |
| 181 | `coherent_tracking_defect` |
| 182 | `microscopic_coherent_slippage` |
| 183 | `triangular_normal_form` |
| 184 | `branch_invariant_coordinates` |
| 185 | `microscopic_monomials` |
| 186 | `similarity_orbit_closure` |
| 187 | `orbit_quotient_closure` |
| 188 | `branch_observables_completion` |
| 189 | `transfer_shape_prefactor_compiler` |
| 190 | `direct_defect_vs_dressing_split` |
| 191 | `minimal_pde_data_packet` |
| 192 | `orbit_quotient_projectors` |
| 193 | `isotropic_grouped_p2_target_surface` |
| 194 | `outgoing_l2_fingerprint_and_deformation_algebra` |
| 195 | `source_map_reduction_of_canonical_outgoing_branch` |
| 196 | `higher_odd_irrelevance_theorem` |
| 197 | `conditional_packetA_closure_theorem` |
| 198 | `exact_finite_orbit_law` |
| 199 | `pairwise_orbit_transport_law` |
| 200 | `reference_free_home_stretch_theorem` |
| 201 | `explicit_realization_compiler_and_canonical_orbit_projection` |
| 202 | `free_quintuple_target_graph` |
| 203 | `free_quintuple_scalar_closure_slice_and_crossing_theorem` |
| 204 | `explicit_log_ray_sweep_and_scalar_root_predictor` |
| 218 | `full_support_cardinality_5_completion_and_local_mixed_ray_se` |
| 219 | `one_port_mixed_bundle_static_kernel_and_square_law_suppressi` |
| 220 | `dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survi` |
| 221 | `resonance_linewidth_tradeoff_dispersive_no_free_lunch_theore` |
| 222 | `concrete_finite_throat_primitive_branch_pole_census_and_resi` |
| 223 | `5pn_isotropic_target_surface_primitive_branch_compatibility_` |
| 224 | `pde_branch_packet_compiler_weak_axisymmetric_ceiling_transpo` |
| 225 | `microscopic_xi1_compiler_first_order_conservative_compensati` |
| 226 | `strict_5pn_even_gate_package_surviving_mixed_corridor_and_pu` |
| 227 | `pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_` |
| 228 | `numerator_denominator_split_of_the_pure_transfer_corridor_an` |
| 229 | `selected_branch_numerator_denominator_signature_and_softenin` |
| 230 | `selected_branch_classifier_to_dynamic_window_compiler_and_st` |
| 231 | `continuum_placement_pullback_of_the_selected_branch_dynamic_` |
| 232 | `known_5pn_data_injection_and_current_branch_verdict` |
| 233 | `rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_ga` |
| 234 | `direct_branch_observable_static_gate_and_the_two_observable_` |
| 235 | `rigid_mouth_packet_projectors_static_blind_dressing_line_and` |
| 236 | `rigid_mouth_microscopic_dependent_plane_projectors_equal_dri` |
| 237 | `actual_branch_dressing_compiler_finite_static_blind_curve_an` |
| 238 | `physical_branch_transfer_shape_compiler_packet_factorization` |
| 239 | `rigid_mouth_physical_normal_form_exact_physical_to_microscop` |
| 240 | `selected_branch_loading_ratio_from_the_minimal_isotropic_qua` |
| 241 | `exact_primitive_ranking_on_the_selected_twin_support_branch` |
| 242 | `actual_twin_support_placement_and_coherent_orbit_lock_compil` |

(Stages outside these slices follow the same one-to-one filename rule.)

## CLASS A -- stale `stageNNN_<stem>` citations (deterministic, stem-keyed)

**25 findings. ALL HIGH confidence** (stem maps one-to-one to canonical; cited number != canonical and no real file exists for the cited token).

| file | line | stale citation | -> canonical | cited | conf |
| --- | --- | --- | --- | --- | --- |
| `notes/stages/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md` | 402 | `stage244_isotropic_grouped_p2_target_surface_sympy_audit` | **193** | 244 | high |
| `notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md` | 497 | `stage244_pure_transfer_load_factor_outgoing_rigidity_sieve_a` | **227** | 244 | high |
| `notes/stages/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md` | 486 | `stage253_one_port_mixed_bundle_static_kernel_and_square_law_` | **219** | 253 | high |
| `notes/stages/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md` | 408 | `stage252_rigid_mouth_packet_projectors_static_blind_dressing` | **235** | 252 | high |
| `notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md` | 421 | `stage242_minimal_pde_data_packet_sympy_audit` | **191** | 242 | high |
| `notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md` | 422 | `stage242_minimal_pde_data_packet_sympy_audit_output` | **191** | 242 | high |
| `notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` | 36 | `stage252_full_support_cardinality_5_completion_and_local_mix` | **218** | 252 | high |
| `notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` | 39 | `stage252_full_support_cardinality_5_completion_and_local_mix` | **218** | 252 | high |
| `notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md` | 478 | `stage240_5pn_isotropic_target_surface_primitive_branch_compa` | **223** | 240 | high |
| `notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md` | 484 | `stage245_numerator_denominator_split_of_the_pure_transfer_co` | **228** | 245 | high |
| `notes/stages/moving_throat_pde_stage096_geometry_lane_check_verdict.md` | 12 | `stage147_geometry_lane_check_verdict_sympy_audit` | **096** | 147 | high |
| `notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md` | 410 | `stage243_orbit_quotient_projectors_sympy_audit` | **192** | 243 | high |
| `notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md` | 411 | `stage243_orbit_quotient_projectors_sympy_audit_output` | **192** | 243 | high |
| `notes/stages/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md` | 496 | `stage253_rigid_mouth_microscopic_dependent_plane_projectors_` | **236** | 253 | high |
| `notes/stages/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md` | 437 | `stage251_direct_branch_observable_static_gate_and_the_two_ob` | **234** | 251 | high |
| `notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md` | 59 | `stage251_reference_free_home_stretch_theorem_sympy_audit` | **200** | 251 | high |
| `notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md` | 62 | `stage251_reference_free_home_stretch_theorem_mathematica_aud` | **200** | 251 | high |
| `notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md` | 551 | `stage253_free_quintuple_target_graph_sympy_audit` | **202** | 253 | high |
| `notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md` | 552 | `stage253_free_quintuple_target_graph_sympy_audit_output` | **202** | 253 | high |
| `notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md` | 458 | `stage239_concrete_finite_throat_primitive_branch_pole_census` | **222** | 239 | high |
| `notes/stages/moving_throat_pde_stage090_updated_reduced_status.md` | 10 | `stage141_updated_reduced_status_sympy_audit` | **090** | 141 | high |
| `notes/stages/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md` | 45 | `stage237_free_quintuple_scalar_closure_slice_and_crossing_th` | **203** | 237 | high |
| `notes/stages/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md` | 48 | `stage237_free_quintuple_scalar_closure_slice_and_crossing_th` | **203** | 237 | high |
| `notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md` | 480 | `stage246_selected_branch_numerator_denominator_signature_and` | **229** | 246 | high |
| `notes/stages/review/PHYSICS_COEVOLUTION_AUDIT.md` | 248 | `stage170_outgoing_load_factorization` | **176** | 170 | high |

## CLASS B -- notes self-title / H1 drift

**19 findings. ALL HIGH confidence** (target = the note file's own canonical stage from its filename stem).

| file | line | title says | -> canonical (own) | conf | title text |
| --- | --- | --- | --- | --- | --- |
| `notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md` | 2 | Stage 099 | **116** | high | # Moving-Throat PDE — Stage 99: Finite D/N Mixed-Tube Realization |
| `notes/stages/moving_throat_pde_stage114_concrete_core_schur.md` | 2 | Stage 097 | **114** | high | # Moving-Throat PDE — Stage 97: Concrete Two-Channel Core Outlet Model |
| `notes/stages/moving_throat_pde_stage119_parent_balance_family.md` | 2 | Stage 221 | **119** | high | # Moving-Throat PDE — Stage 221: One-Parameter Parent Compensation Fam |
| `notes/stages/moving_throat_pde_stage069_final_reduced_verdict.md` | 2 | Stage 052 | **069** | high | # Moving-Throat PDE — Stage 52: Final Reduced Verdict for the Support/ |
| `notes/stages/moving_throat_pde_stage068_resonance_thresholds.md` | 2 | Stage 051 | **068** | high | # Moving-Throat PDE — Stage 51: Resonance-Corrected Thresholds for the |
| `notes/stages/moving_throat_pde_stage183_triangular_normal_form.md` | 2 | Stage 251 | **183** | high | # Moving-Throat PDE — Stage 251: Exact Triangular Normal Form of the C |
| `notes/stages/moving_throat_pde_stage172_physical_slope_collapse.md` | 3 | Stage 240 | **172** | high | # Moving-Throat PDE — Stage 240: Collapse of the Linear Grouped Outlet |
| `notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md` | 2 | Stage 240 | **223** | high | # Moving-Throat PDE — Stage 240: `5`PN Isotropic Target Surface, Primi |
| `notes/stages/moving_throat_pde_stage125_positive_source_theorem.md` | 2 | Stage 227 | **125** | high | # Moving-Throat PDE — Stage 227: Positive Local Mouth-Source Theorem |
| `notes/stages/moving_throat_pde_stage126_positive_source_families.md` | 2 | Stage 228 | **126** | high | # Moving-Throat PDE — Stage 228: Explicit Positive Source Families and |
| `notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md` | 2 | Stage 237 | **186** | high | # Moving-Throat PDE — Stage 237: Exact Microscopic Similarity Orbit an |
| `notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md` | 2 | Stage 249 | **232** | high | # Moving-Throat PDE — Stage 249: Known 5PN Data Injection and Current  |
| `notes/stages/moving_throat_pde_stage181_coherent_tracking_defect.md` | 2 | Stage 249 | **181** | high | # Moving-Throat PDE — Stage 249: Coherent Tracking-Branch Weak-Axisymm |
| `notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md` | 2 | Stage 050 | **067** | high | # Moving-Throat PDE — Stage 50: Exact Sech–Gaussian Coherence Resonanc |
| `notes/stages/moving_throat_pde_stage185_microscopic_monomials.md` | 2 | Stage 253 | **185** | high | # Moving-Throat PDE — Stage 253: Direct Microscopic Monomial Coordinat |
| `notes/stages/moving_throat_pde_stage118_parent_core_extraction.md` | 2 | Stage 220 | **118** | high | # Moving-Throat PDE — Stage 220: Parent-Action Extraction of Core Para |
| `notes/stages/moving_throat_pde_stage120_core_parameter_status.md` | 2 | Stage 222 | **120** | high | # Moving-Throat PDE — Stage 222: Core-Parameter Extraction Status |
| `notes/stages/moving_throat_pde_stage115_core_balance_compensation.md` | 2 | Stage 098 | **115** | high | # Moving-Throat PDE — Stage 98: Exact Core-Balance Compensation Theore |
| `notes/stages/moving_throat_pde_stage117_outlet_core_status.md` | 2 | Stage 219 | **117** | high | # Moving-Throat PDE — Stage 219: Concrete Outlet-Core Status |

## CLASS C -- script / paper-card self-banner drift

**155 script-banner findings (122 distinct scripts) + 34 published paper-card self-titles = 189 total. ALL HIGH confidence** (banner/title identifies the OWN stage; target = filename/label canonical). Script-banner subset is uniformly +17; paper-card titles uniformly -17 (band 091-124). Offset uniformity is descriptive only -- fixes MUST still be file-keyed, NOT swept.

### C.1 Script self-banners (.py / .wl)

| file | line | banner says | -> canonical (own) | conf | banner text |
| --- | --- | --- | --- | --- | --- |
| `scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py` | 35 | Stage 177 | **194** | high | banner("STAGE 177 — EXACT OUTGOING l=2 DTN FINGERPRINT, FIXI |
| `scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py` | 208 | Stage 177 | **194** | high | banner("STAGE 177 LEDGER") |
| `scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py` | 28 | Stage 189 | **206** | high | banner("STAGE 189 — CERTIFIED RAY RANKING AND LOCAL BRACKETI |
| `scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py` | 228 | Stage 189 | **206** | high | banner("STAGE 189 SYMPY AUDIT PASSED") |
| `scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py` | 35 | Stage 171 | **188** | high | banner("STAGE 171 — BRANCH-OBSERVABLE COMPLETION AND THE EXA |
| `scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py` | 198 | Stage 171 | **188** | high | banner("STAGE 171 LEDGER") |
| `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py` | 5 | Stage 008 | **025** | high | SymPy audit for Stage 8 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py` | 31 | Stage 031 | **048** | high | banner("STAGE 31 — SUPPORT COMPENSATION THEOREM AUDIT") |
| `scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py` | 25 | Stage 037 | **054** | high | banner("STAGE 37 — ROBIN-COMPLIANCE SOFTENING") |
| `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py` | 182 | Stage 013 | **030** | high | banner("STAGE 13 AUDIT COMPLETE") |
| `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py` | 39 | Stage 029 | **046** | high | banner("STAGE 29 — TRACKING-BRANCH BOUNDS AUDIT") |
| `scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` | 41 | Stage 197 | **214** | high | banner("STAGE 197 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OP |
| `scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` | 356 | Stage 197 | **214** | high | banner("STAGE 197 SYMPY AUDIT COMPLETED SUCCESSFULLY") |
| `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py` | 51 | Stage 024 | **041** | high | banner("STAGE 24 — RANK-2 SUPPORT COMPLETION") |
| `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py` | 55 | Stage 026 | **043** | high | banner("STAGE 26 — CONTINUUM EXTRACTION OF THE ACTUAL SUPPOR |
| `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py` | 32 | Stage 030 | **047** | high | banner("STAGE 30 — COHERENT KERNEL MAP AUDIT") |
| `scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py` | 51 | Stage 025 | **042** | high | banner("STAGE 25 — SELECTED-MODE NORMALIZATION UNDER RANK-2  |
| `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py` | 62 | Stage 022 | **039** | high | banner("STAGE 22 — SPLIT-U SECTOR: EXACT CONTINUUM REDUCTION |
| `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py` | 51 | Stage 027 | **044** | high | banner("STAGE 27 — CONTINUUM-SELECTED RANK-2 CLOSURE") |
| `scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py` | 5 | Stage 010 | **027** | high | SymPy audit for Stage 10 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py` | 35 | Stage 167 | **184** | high | banner("STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES") |
| `scripts/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.py` | 35 | Stage 191 | **208** | high | banner("STAGE 191 — PAIRWISE MIXED RAYS AND OFF-DIAGONAL HES |
| `scripts/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.py` | 219 | Stage 191 | **208** | high | banner("STAGE 191 SYMPY AUDIT COMPLETED SUCCESSFULLY") |
| `scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py` | 40 | Stage 182 | **199** | high | banner("STAGE 182 — EXACT PAIRWISE ORBIT-TRANSPORT LAW AND T |
| `scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py` | 403 | Stage 182 | **199** | high | banner("STAGE 182 LEDGER") |
| `scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py` | 35 | Stage 187 | **204** | high | banner("STAGE 187 — EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT P |
| `scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py` | 225 | Stage 187 | **204** | high | banner("STAGE 187 SYMPY AUDIT PASSED") |
| `scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py` | 36 | Stage 021 | **038** | high | banner("STAGE 21 — DIMENSIONLESS CONTINUUM PLACEMENT AUDIT") |
| `scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py` | 224 | Stage 021 | **038** | high | banner("STAGE 21 AUDIT COMPLETE") |
| `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py` | 45 | Stage 042 | **059** | high | banner("STAGE 42 — OPERATOR-SELECTED RESIDUAL BOUNDS") |
| `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py` | 136 | Stage 201 | **218** | high | banner("STAGE 201 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MI |
| `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py` | 49 | Stage 047 | **064** | high | banner("STAGE 47 — PARENT EQUILIBRIUM SOURCE/SUPPORT ALIGNME |
| `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py` | 176 | Stage 047 | **064** | high | banner("STAGE 47 AUDIT PASSED") |
| `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` | 51 | Stage 020 | **037** | high | banner("STAGE 20 — CONTINUUM-KERNEL EXTRACTION AUDIT") |
| `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` | 222 | Stage 020 | **037** | high | banner("STAGE 20 AUDIT COMPLETE") |
| `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py` | 5 | Stage 006 | **023** | high | SymPy audit for Stage 6 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py` | 35 | Stage 193 | **210** | high | banner("STAGE 193 — THREE-COORDINATE MIXED-SIMPLEX AND THE C |
| `scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py` | 256 | Stage 193 | **210** | high | banner("STAGE 193 SYMPY AUDIT COMPLETED SUCCESSFULLY") |
| `scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py` | 42 | Stage 166 | **183** | high | banner("STAGE 166 — TRIANGULAR NORMAL FORM OF THE COHERENT D |
| `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py` | 5 | Stage 007 | **024** | high | SymPy audit for Stage 7 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py` | 53 | Stage 023 | **040** | high | banner("STAGE 23 — GENERALIZED SELECTED-BRANCH NORMALIZATION |
| `scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py` | 25 | Stage 036 | **053** | high | banner("STAGE 36 — EXACT OVERLAP-BOOST WINDOW") |
| `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py` | 45 | Stage 049 | **066** | high | banner("STAGE 49 — DIMENSIONLESS WALL FIGURE OF MERIT") |
| `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py` | 101 | Stage 049 | **066** | high | banner("STAGE 49 AUDIT PASSED") |
| `scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py` | 51 | Stage 032 | **049** | high | banner("STAGE 32 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PH |
| `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py` | 34 | Stage 033 | **050** | high | banner("STAGE 33 — PHYSICAL ZETA VS ZETA_REQ") |
| `scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py` | 5 | Stage 034 | **051** | high | SymPy audit for Stage 34 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py` | 63 | Stage 034 | **051** | high | banner("STAGE 34 — EXACT TRACKING-BRANCH PRODUCT") |
| `scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py` | 35 | Stage 179 | **196** | high | banner("STAGE 179 — EXACT HIGHER-ODD IRRELEVANCE THEOREM") |
| `scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py` | 189 | Stage 179 | **196** | high | banner("STAGE 179 LEDGER") |
| `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py` | 51 | Stage 048 | **065** | high | banner("STAGE 48 — THIN-WALL CONFINEMENT BRANCH") |
| `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py` | 185 | Stage 048 | **065** | high | banner("STAGE 48 AUDIT PASSED") |
| `scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py` | 49 | Stage 176 | **193** | high | banner("STAGE 176 — EXACT ISOTROPIC GROUPED-P2 TARGET SURFAC |
| `scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py` | 164 | Stage 176 | **193** | high | banner("STAGE 176 LEDGER") |
| `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py` | 5 | Stage 053 | **070** | high | SymPy audit for Stage 53: |
| `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py` | 24 | Stage 053 | **070** | high | banner("STAGE 53 — EXPLICIT GNLS WALL-SHELL REDUCTION") |
| `scripts/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.py` | 35 | Stage 181 | **198** | high | banner("STAGE 181 — EXACT FINITE ORBIT LAW FOR THE DEPENDENT |
| `scripts/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.py` | 240 | Stage 181 | **198** | high | banner("STAGE 181 LEDGER") |
| `scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` | 35 | Stage 178 | **195** | high | banner("STAGE 178 — EXACT SOURCE-MAP REDUCTION OF THE CANONI |
| `scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` | 189 | Stage 178 | **195** | high | banner("STAGE 178 LEDGER") |
| `scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py` | 5 | Stage 071 | **088** | high | SymPy audit for Stage 71. |
| `scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py` | 32 | Stage 162 | **179** | high | banner("STAGE 162 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM") |
| `scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py` | 22 | Stage 044 | **061** | high | banner("STAGE 44 — MICROSCOPIC GAIN THRESHOLDS") |
| `scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py` | 31 | Stage 039 | **056** | high | banner("STAGE 39 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYM |
| `scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py` | 65 | Stage 052 | **069** | high | banner("STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT") |
| `scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` | 35 | Stage 200 | **217** | high | banner("STAGE 200 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR  |
| `scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py` | 76 | Stage 183 | **200** | high | banner("STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEORE |
| `scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py` | 439 | Stage 183 | **200** | high | banner("STAGE 183 LEDGER") |
| `scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py` | 42 | Stage 196 | **213** | high | banner("STAGE 196 — FOUR-COORDINATE MIXED SIMPLEX AND THE SU |
| `scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py` | 402 | Stage 196 | **213** | high | banner("STAGE 196 SYMPY AUDIT COMPLETED SUCCESSFULLY") |
| `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py` | 37 | Stage 051 | **068** | high | banner("STAGE 51 — RESONANCE-CORRECTED THRESHOLDS") |
| `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py` | 108 | Stage 014 | **031** | high | banner("STAGE 14 AUDIT COMPLETE") |
| `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py` | 5 | Stage 009 | **026** | high | SymPy audit for Stage 9 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py` | 22 | Stage 043 | **060** | high | banner("STAGE 43 — ENTROPIC SOURCE MICROCLOSURE") |
| `scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py` | 6 | Stage 035 | **052** | high | SymPy audit for Stage 35 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py` | 42 | Stage 035 | **052** | high | banner("STAGE 35 — EXACT BRANCH-PRODUCT REGIME CLASSIFICATIO |
| `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py` | 5 | Stage 054 | **071** | high | SymPy audit for Stage 54: |
| `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py` | 24 | Stage 054 | **071** | high | banner("STAGE 54 — CANONICAL TANH-WALL BRANCH") |
| `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py` | 31 | Stage 046 | **063** | high | banner("STAGE 46 — PARENT-OVERLAP THRESHOLD THEOREM") |
| `scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py` | 36 | Stage 170 | **187** | high | banner("STAGE 170 — EXACT ORBIT–QUOTIENT CLOSURE") |
| `scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py` | 35 | Stage 192 | **209** | high | banner("STAGE 192 — EXACT PAIRWISE RATIO OPTIMIZER AND MIXED |
| `scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py` | 148 | Stage 192 | **209** | high | banner("STAGE 192 SYMPY AUDIT COMPLETED SUCCESSFULLY") |
| `scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py` | 5 | Stage 011 | **028** | high | SymPy audit for Stage 11 of the moving-throat PDE program. |
| `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py` | 31 | Stage 045 | **062** | high | banner("STAGE 45 — PARENT-ACTION PROJECTION OF THE MICROSCOP |
| `scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py` | 49 | Stage 185 | **202** | high | banner("STAGE 185 — EXACT FREE-QUINTUPLE TARGET GRAPH AND TH |
| `scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py` | 255 | Stage 185 | **202** | high | banner("STAGE 185 AUDIT COMPLETE") |
| `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py` | 33 | Stage 041 | **058** | high | banner("STAGE 41 — COUPLED SUPPORT/SOURCE OPERATOR") |
| `scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py` | 65 | Stage 186 | **203** | high | banner("STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND  |
| `scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py` | 364 | Stage 186 | **203** | high | banner("STAGE 186 SYMPY AUDIT PASSED") |
| `scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py` | 25 | Stage 038 | **055** | high | banner("STAGE 38 — EXPLICIT LOWEST-LANE REACHABILITY") |
| `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py` | 33 | Stage 040 | **057** | high | banner("STAGE 40 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP") |
| `scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py` | 35 | Stage 173 | **190** | high | banner("STAGE 173 — DIRECT DEFECT VS DRESSING SPLIT, SUPPORT |
| `scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py` | 269 | Stage 173 | **190** | high | banner("STAGE 173 LEDGER") |
| `scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py` | 35 | Stage 175 | **192** | high | banner("STAGE 175 — EXACT ORBIT/QUOTIENT PROJECTORS AND THE  |
| `scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py` | 196 | Stage 175 | **192** | high | banner("STAGE 175 LEDGER") |
| `scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py` | 41 | Stage 184 | **201** | high | banner("STAGE 184 — EXPLICIT REALIZATION COMPILER AND CANONI |
| `scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py` | 263 | Stage 184 | **201** | high | banner("STAGE 184 LEDGER") |
| `scripts/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py` | 33 | Stage 163 | **180** | high | banner("STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE") |
| `scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py` | 35 | Stage 190 | **207** | high | banner("STAGE 190 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERT |
| `scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py` | 244 | Stage 190 | **207** | high | banner("STAGE 190 SYMPY AUDIT PASSED") |
| `scripts/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.py` | 42 | Stage 180 | **197** | high | banner("STAGE 180 — EXACT CONDITIONAL PACKET-A CLOSURE THEOR |
| `scripts/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.py` | 212 | Stage 180 | **197** | high | banner("STAGE 180 LEDGER") |
| `mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl` | 26 | Stage 044 | **061** | high | banner["STAGE 044 — MICROSCOPIC GAIN THRESHOLDS"]; |
| `mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl` | 26 | Stage 046 | **063** | high | banner["STAGE 046 — PARENT THRESHOLDS"]; |
| `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl` | 31 | Stage 019 | **036** | high | banner["STAGE 019 — SUPPORT-FEASIBILITY FRONTIER"]; |
| `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl` | 68 | Stage 183 | **200** | high | banner["STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEORE |
| `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl` | 304 | Stage 183 | **200** | high | banner["STAGE 183 LEDGER"]; |
| `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` | 61 | Stage 007 | **024** | high | banner["STAGE 007 — OVERLAP ISOTROPY"]; |
| `mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl` | 26 | Stage 039 | **056** | high | banner["STAGE 039 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASY |
| `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl` | 148 | Stage 013 | **030** | high | banner["STAGE 13 AUDIT COMPLETE"]; |
| `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl` | 33 | Stage 023 | **040** | high | banner["STAGE 023 — GENERALIZED SELECTED-BRANCH NORMALIZATIO |
| `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl` | 26 | Stage 051 | **068** | high | banner["STAGE 051 — RESONANCE-CORRECTED THRESHOLDS"]; |
| `mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl` | 32 | Stage 035 | **052** | high | banner["STAGE 035 — NON-TWIN ASYMMETRY THRESHOLD"]; |
| `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl` | 26 | Stage 045 | **062** | high | banner["STAGE 045 — PARENT ACTION GAIN"]; |
| `mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl` | 26 | Stage 031 | **048** | high | banner["STAGE 031 — SUPPORT COMPENSATION THEOREM"]; |
| `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl` | 32 | Stage 043 | **060** | high | banner["STAGE 043 — ENTROPIC MICROCLOSURE"]; |
| `mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl` | 73 | Stage 186 | **203** | high | banner["STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND  |
| `mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl` | 333 | Stage 186 | **203** | high | banner["STAGE 186 MATHEMATICA AUDIT PASSED"]; |
| `mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl` | 26 | Stage 011 | **028** | high | banner["STAGE 011 — LOADED PROFILE SELECTION"]; |
| `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl` | 26 | Stage 018 | **035** | high | banner["STAGE 018 — DIMENSIONLESS NORMALIZATION LOCUS"]; |
| `mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl` | 32 | Stage 049 | **066** | high | banner["STAGE 049 — DIMENSIONLESS WALL FIGURE OF MERIT"]; |
| `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl` | 26 | Stage 167 | **184** | high | banner["STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES"]; |
| `mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl` | 32 | Stage 036 | **053** | high | banner["STAGE 036 — EXACT OVERLAP-BOOST WINDOW"]; |
| `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl` | 37 | Stage 006 | **023** | high | banner["STAGE 006 — FULL GROUPED BUNDLE"]; |
| `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl` | 26 | Stage 053 | **070** | high | banner["STAGE 053 — EXPLICIT GNLS WALL-SHELL REDUCTION"]; |
| `mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl` | 26 | Stage 163 | **180** | high | banner["STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE"]; |
| `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl` | 32 | Stage 041 | **058** | high | banner["STAGE 041 — COUPLED SUPPORT/SOURCE OPERATOR"]; |
| `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl` | 36 | Stage 032 | **049** | high | banner["STAGE 32 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PH |
| `mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl` | 32 | Stage 038 | **055** | high | banner["STAGE 038 — EXPLICIT LOWEST-LANE REACHABILITY"]; |
| `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl` | 26 | Stage 030 | **047** | high | banner["STAGE 030 — COHERENT KERNEL MAP"]; |
| `mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl` | 33 | Stage 021 | **038** | high | banner["STAGE 021 — DIMENSIONLESS CONTINUUM PLACEMENT"]; |
| `mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl` | 33 | Stage 025 | **042** | high | banner["STAGE 025 — RANK-2 SELECTED-MODE NORMALIZATION"]; |
| `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl` | 44 | Stage 042 | **059** | high | banner["STAGE 042 — OPERATOR BRANCH RESIDUAL BOUNDS"]; |
| `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl` | 32 | Stage 033 | **050** | high | banner["STAGE 033 — PHYSICAL ZETA VS ZETA_REQ"]; |
| `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl` | 33 | Stage 026 | **043** | high | banner["STAGE 026 — CONTINUUM EXTRACTION OF THE ACTUAL SUPPO |
| `mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl` | 32 | Stage 037 | **054** | high | banner["STAGE 037 — ROBIN-COMPLIANCE SOFTENING"]; |
| `mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl` | 58 | Stage 048 | **065** | high | banner["STAGE 048 — THIN-WALL CONFINEMENT BRANCH"]; |
| `mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl` | 26 | Stage 054 | **071** | high | banner["STAGE 054 — CANONICAL TANH-WALL BRANCH"]; |
| `mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl` | 32 | Stage 166 | **183** | high | banner["STAGE 166 — TRIANGULAR NORMAL FORM OF THE COHERENT D |
| `mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl` | 32 | Stage 029 | **046** | high | banner["STAGE 029 — TRACKING-BRANCH BOUNDS"]; |
| `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl` | 35 | Stage 004 | **021** | high | banner["STAGE 004 — MAXWELL + MIXED-SECTOR REDUCTION"]; |
| `mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl` | 26 | Stage 017 | **034** | high | banner["STAGE 017 — SOFTENING-DEPTH NORMAL FORM"]; |
| `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl` | 26 | Stage 047 | **064** | high | banner["STAGE 047 — EQUILIBRIUM ALIGNMENT"]; |
| `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl` | 33 | Stage 022 | **039** | high | banner["STAGE 022 — SPLIT-U SECTOR"]; |
| `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl` | 35 | Stage 204 | **221** | high | banner["STAGE 204 — RESONANCE / LINEWIDTH TRADEOFF"]; |
| `mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl` | 40 | Stage 020 | **037** | high | banner["STAGE 020 — CONTINUUM-KERNEL EXTRACTION"]; |
| `mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl` | 33 | Stage 052 | **069** | high | banner["STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT"]; |
| `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl` | 26 | Stage 016 | **033** | high | banner["STAGE 016 — MICROSCOPIC NORMALIZATION EQUATION"]; |
| `mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl` | 26 | Stage 162 | **179** | high | banner["STAGE 162 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM"] |
| `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl` | 32 | Stage 034 | **051** | high | banner["STAGE 034 — LOWEST TWIN CRITERION"]; |
| `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl` | 33 | Stage 027 | **044** | high | banner["STAGE 027 — CONTINUUM-SELECTED RANK-2 CLOSURE"]; |
| `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl` | 26 | Stage 005 | **022** | high | banner["STAGE 005 — GROUPED P2 NORMALIZATION BRIDGE"]; |
| `mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl` | 26 | Stage 170 | **187** | high | banner["STAGE 170 — EXACT ORBIT-QUOTIENT CLOSURE"]; |
| `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl` | 32 | Stage 040 | **057** | high | banner["STAGE 040 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP" |
| `mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl` | 33 | Stage 024 | **041** | high | banner["STAGE 024 — RANK-2 SUPPORT COMPLETION"]; |

### C.2 Published paper-card self-titles (paper/stages/stage_NNN.tex; section Stage MMM != label/filename)

| file | line | title says | -> canonical (label/filename) | conf |
| --- | --- | --- | --- | --- |
| `paper/stages/stage_091.tex` | 1 | Stage 108 | **091** | high |
| `paper/stages/stage_092.tex` | 1 | Stage 109 | **092** | high |
| `paper/stages/stage_093.tex` | 1 | Stage 110 | **093** | high |
| `paper/stages/stage_094.tex` | 1 | Stage 111 | **094** | high |
| `paper/stages/stage_095.tex` | 1 | Stage 112 | **095** | high |
| `paper/stages/stage_096.tex` | 1 | Stage 113 | **096** | high |
| `paper/stages/stage_097.tex` | 1 | Stage 114 | **097** | high |
| `paper/stages/stage_098.tex` | 1 | Stage 115 | **098** | high |
| `paper/stages/stage_099.tex` | 1 | Stage 116 | **099** | high |
| `paper/stages/stage_100.tex` | 1 | Stage 117 | **100** | high |
| `paper/stages/stage_101.tex` | 1 | Stage 118 | **101** | high |
| `paper/stages/stage_102.tex` | 1 | Stage 119 | **102** | high |
| `paper/stages/stage_103.tex` | 1 | Stage 120 | **103** | high |
| `paper/stages/stage_104.tex` | 1 | Stage 121 | **104** | high |
| `paper/stages/stage_105.tex` | 1 | Stage 122 | **105** | high |
| `paper/stages/stage_106.tex` | 1 | Stage 123 | **106** | high |
| `paper/stages/stage_107.tex` | 1 | Stage 124 | **107** | high |
| `paper/stages/stage_108.tex` | 1 | Stage 125 | **108** | high |
| `paper/stages/stage_109.tex` | 1 | Stage 126 | **109** | high |
| `paper/stages/stage_110.tex` | 1 | Stage 127 | **110** | high |
| `paper/stages/stage_111.tex` | 1 | Stage 128 | **111** | high |
| `paper/stages/stage_112.tex` | 1 | Stage 129 | **112** | high |
| `paper/stages/stage_113.tex` | 1 | Stage 130 | **113** | high |
| `paper/stages/stage_114.tex` | 1 | Stage 131 | **114** | high |
| `paper/stages/stage_115.tex` | 1 | Stage 132 | **115** | high |
| `paper/stages/stage_116.tex` | 1 | Stage 133 | **116** | high |
| `paper/stages/stage_117.tex` | 1 | Stage 134 | **117** | high |
| `paper/stages/stage_118.tex` | 1 | Stage 135 | **118** | high |
| `paper/stages/stage_119.tex` | 1 | Stage 136 | **119** | high |
| `paper/stages/stage_120.tex` | 1 | Stage 137 | **120** | high |
| `paper/stages/stage_121.tex` | 1 | Stage 138 | **121** | high |
| `paper/stages/stage_122.tex` | 1 | Stage 139 | **122** | high |
| `paper/stages/stage_123.tex` | 1 | Stage 140 | **123** | high |
| `paper/stages/stage_124.tex` | 1 | Stage 141 | **124** | high |

## CLASS D -- bare-prose body cross-references (CONTENT-dependent, NOT stem-keyed)

These are `Stage NNN` mentions in note/script BODY prose referencing a DIFFERENT stage's deliverable (or this stage's own role using a pre-realignment number). They require CONTENT mapping and are **NOT safe for a deterministic sweep**. The raw pool is large: **121 high-signal self-role body cross-refs** (across 57 notes) of the form `Stage N <keeps|closed|reduced|insertion|...>`, plus ~1,200 looser forward `Stage N` mentions project-wide. Most resolve by recognizing N as a stale OLD number for an UPSTREAM stage (target = that stage's canonical, typically N-17 in the band) -- but MUST be content-confirmed, never swept.

**121 flagged Class-D candidates below** (best-guess target = needs content map):

| file | line | cites Stage | target | conf | context |
| --- | --- | --- | --- | --- | --- |
| `moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md` | 453 | 242 | needs content map | medium | Stage 242 reduced the first honest mixed-ray search to  |
| `moving_throat_pde_stage169_no_linear_p2_scalar_slippage.md` | 5 | 253 | needs content map | medium | Stage 253 reduced the first-order off-bundle defect to  |
| `moving_throat_pde_stage175_wall_normalized_load_shape.md` | 5 | 242 | needs content map | medium | Stage 242 reduced the remaining linear grouped `2.5`PN  |
| `moving_throat_pde_stage175_wall_normalized_load_shape.md` | 382 | 241 | needs content map | medium | Stage 241 already proved that, on the even-preserving b |
| `moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md` | 22 | 243 | needs content map | medium | Stage 243 already reduced the live same-charge corridor |
| `moving_throat_pde_stage171_microscopic_grouped_obstructions.md` | 100 | 238 | needs content map | medium | Stage 238 already gave |
| `moving_throat_pde_stage171_microscopic_grouped_obstructions.md` | 476 | 238 | needs content map | medium | Stage 238 already reduced the linear grouped-anisotropy |
| `moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md` | 149 | 244 | needs content map | medium | Stage 244 already gave the first-order defect field |
| `moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md` | 252 | 241 | needs content map | medium | Stage 241 already showed that on the weak-axisymmetric  |
| `moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md` | 364 | 244 | needs content map | medium | Stage 244 reduced the remaining grouped weak-axisymmetr |
| `moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md` | 56 | 239 | needs content map | medium | Stage 239 already defined the logarithmic scalar |
| `moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md` | 301 | 239 | needs content map | medium | Stage 239 already identified the turning-point criterio |
| `moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md` | 14 | 252 | needs content map | medium | This is the first post-Stage 252 insertion of actual re |
| `moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md` | 22 | 252 | needs content map | medium | Stage 252 closed the local mixed-ray search sieve. |
| `moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md` | 41 | 253 | needs content map | medium | So Stage 253 keeps the same-charge mixed corridor alive |
| `moving_throat_pde_stage173_axisymmetric_loading_mismatch.md` | 5 | 240 | needs content map | medium | Stage 240 reduced the linear grouped outlet problem to  |
| `moving_throat_pde_stage173_axisymmetric_loading_mismatch.md` | 191 | 240 | needs content map | medium | Stage 240 already translated the one-parameter hidden-e |
| `moving_throat_pde_stage183_triangular_normal_form.md` | 6 | 250 | needs content map | medium | Stage 250 reduced the coherent weak-axisymmetric groupe |
| `moving_throat_pde_stage183_triangular_normal_form.md` | 85 | 250 | needs content map | medium | Stage 250 already gave the exact tracking/nontracking s |
| `moving_throat_pde_stage183_triangular_normal_form.md` | 142 | 250 | needs content map | medium | Stage 250 already gave |
| `moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md` | 14 | 243 | needs content map | medium | Stage 243 closed the pairwise problem exactly. Every mo |
| `moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md` | 550 | 243 | needs content map | medium | Stage 243 already proved that the pairwise boundary pro |
| `moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md` | 284 | 246 | needs content map | medium | Stage 246 already reduced the odd observable closure to |
| `moving_throat_pde_stage174_static_self_similarity.md` | 5 | 241 | needs content map | medium | Stage 241 reduced the remaining linear grouped `2.5`PN  |
| `moving_throat_pde_stage174_static_self_similarity.md` | 367 | 241 | needs content map | medium | Stage 241 already proved that on the even-preserving br |
| `moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md` | 40 | 253 | needs content map | medium | and Stage 253 already gave their exact microscopic depe |
| `moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md` | 61 | 253 | needs content map | medium | Stage 253 already proved that on the rigid-mouth slice  |
| `moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md` | 104 | 239 | needs content map | medium | Stage 239 already fixed the corrected nontracking branc |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 14 | 246 | needs content map | medium | Stage 246 closed the full support-`<=3` local search. |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 136 | 246 | needs content map | medium | So `\(F_{ijk}\)` is exactly the Stage 246 closed simple |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 155 | 246 | needs content map | medium | So `\(F_{ijl}\)` is exactly the Stage 246 closed simple |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 174 | 246 | needs content map | medium | So `\(F_{ikl}\)` is exactly the Stage 246 closed simple |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 193 | 246 | needs content map | medium | So `\(F_{jkl}\)` is exactly the Stage 246 closed simple |
| `moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md` | 688 | 246 | needs content map | medium | Stage 246 already proved that the support-`<=3` search  |
| `moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md` | 25 | 251 | needs content map | medium | Stage 251 already proved that on the direct coherent br |
| `moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md` | 356 | 251 | needs content map | medium | Stage 251 already showed that first-order same-charge s |
| `moving_throat_pde_stage172_physical_slope_collapse.md` | 7 | 239 | needs content map | medium | Stage 239 reduced the weak grouped-lane outlet problem  |
| `moving_throat_pde_stage172_physical_slope_collapse.md` | 120 | 239 | needs content map | medium | Stage 239 reduced the direct grouped outlet obstruction |
| `moving_throat_pde_stage172_physical_slope_collapse.md` | 265 | 238 | needs content map | medium | Stage 238 already gave the direct outlet map |
| `moving_throat_pde_stage167_bundle_transport_tangent_compensation.md` | 5 | 251 | needs content map | medium | Stage 251 reduced the unresolved first-order branch tra |
| `moving_throat_pde_stage167_bundle_transport_tangent_compensation.md` | 114 | 250 | needs content map | medium | Stage 250 already fixed the lower-branch transport laws |
| `moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` | 14 | 249 | needs content map | medium | Stage 249 already finished the global support-`<=4` sea |
| `moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` | 287 | 251 | needs content map | medium | But this is no longer a continuum ambiguity. Stage 251  |
| `moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md` | 367 | 251 | needs content map | medium | - The only new support-`5` content was the unique inter |
| `moving_throat_pde_stage187_orbit_quotient_closure.md` | 5 | 237 | needs content map | medium | Stage 237 already showed that the coherent weak-axisymm |
| `moving_throat_pde_stage187_orbit_quotient_closure.md` | 108 | 237 | needs content map | medium | Stage 237 already constructed a five-parameter multipli |
| `moving_throat_pde_stage180_effective_transfer_shape_collapse.md` | 5 | 247 | needs content map | medium | Stage 247 reduced the remaining weak-axisymmetric group |
| `moving_throat_pde_stage180_effective_transfer_shape_collapse.md` | 61 | 247 | needs content map | medium | Stage 247 already gave the exact portwise factorization |
| `moving_throat_pde_stage180_effective_transfer_shape_collapse.md` | 362 | 241 | needs content map | medium | Stage 241 already fixed the grouped lane pattern of the |
| `moving_throat_pde_stage180_effective_transfer_shape_collapse.md` | 388 | 247 | needs content map | medium | Stage 247 already reduced the last grouped defect to th |
| `moving_throat_pde_stage164_microscopic_log_channels.md` | 5 | 248 | needs content map | medium | Stage 248 reduced the first off-family defect of the co |
| `moving_throat_pde_stage197_conditional_packetA_closure_theorem.md` | 30 | 246 | needs content map | medium | Stage 246 reduced the observable normalization finish l |
| `moving_throat_pde_stage197_conditional_packetA_closure_theorem.md` | 227 | 245 | needs content map | medium | Stage 245 already gave the exact isotropic DtN deformat |
| `moving_throat_pde_stage197_conditional_packetA_closure_theorem.md` | 327 | 247 | needs content map | medium | Because Stage 247 already proved that every extra isotr |
| `moving_throat_pde_stage182_microscopic_coherent_slippage.md` | 5 | 249 | needs content map | medium | Stage 249 reduced the weak-axisymmetric grouped defect  |
| `moving_throat_pde_stage182_microscopic_coherent_slippage.md` | 102 | 249 | needs content map | medium | but Stage 249 already showed that \(\zeta\) drops out o |
| `moving_throat_pde_stage182_microscopic_coherent_slippage.md` | 216 | 249 | needs content map | medium | enter the coherent branch only through \(\zeta\), and S |
| `moving_throat_pde_stage182_microscopic_coherent_slippage.md` | 347 | 249 | needs content map | medium | Stage 249 already showed that the tracking factor is |
| `moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md` | 24 | 244 | needs content map | medium | Stage 244 already reduced the live same-charge corridor |
| `moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md` | 100 | 244 | needs content map | medium | Stage 244 already reduced the strict same-charge surviv |
| `moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md` | 14 | 249 | needs content map | medium | Stage 249 closed the full support-`<=4` local search. |
| `moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md` | 91 | 249 | needs content map | medium | Let the imported Stage 249 closed-simplex intervals be |
| `moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md` | 327 | 250 | needs content map | medium | Stage 250 keeps only the first exact gate needed before |
| `moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md` | 494 | 249 | needs content map | medium | Stage 249 closed the support-`<=4` ledger. |
| `moving_throat_pde_stage192_orbit_quotient_projectors.md` | 108 | 242 | needs content map | medium | So the finite Packet B of Stage 242 is not abstract. It |
| `moving_throat_pde_stage134_family1_mouth_fixedpoint.md` | 70 | 233 | needs content map | medium | Stage 233 already fixed the exact source-shape compensa |
| `moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md` | 442 | 252 | needs content map | medium | Stage 252 already said that the first static same-charg |
| `moving_throat_pde_stage168_off_bundle_slippage.md` | 216 | 248 | needs content map | medium | Stage 248 already split the mouth-bias variation as |
| `moving_throat_pde_stage168_off_bundle_slippage.md` | 286 | 248 | needs content map | medium | Stage 248 already gave |
| `moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md` | 416 | 248 | needs content map | medium | Stage 248 already reduced every primitive quadruple int |
| `moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md` | 433 | 246 | needs content map | medium | Stage 246 already closed the support-`<=3` search with  |
| `moving_throat_pde_stage166_bundle_inversion_four_drifts.md` | 5 | 250 | needs content map | medium | Stage 250 reduced the explicit lower compensated Family |
| `moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md` | 95 | 249 | needs content map | medium | Let the imported Stage 249 closed-simplex intervals be |
| `moving_throat_pde_stage122_mouth_source_compensation_test.md` | 34 | 221 | needs content map | medium | Stage 221 already showed that the compensated canonical |
| `moving_throat_pde_stage200_reference_free_home_stretch_theorem.md` | 411 | 242 | needs content map | medium | Stage 242 reduced the endgame to two packets, but not y |
| `moving_throat_pde_stage186_similarity_orbit_closure.md` | 6 | 253 | needs content map | medium | Stage 253 reduced the coherent weak-axisymmetric proble |
| `moving_throat_pde_stage186_similarity_orbit_closure.md` | 60 | 253 | needs content map | medium | Stage 253 already reduced the three branch-adapted coor |
| `moving_throat_pde_stage186_similarity_orbit_closure.md` | 351 | 251 | needs content map | medium | Stage 251 already proved |
| `moving_throat_pde_stage202_free_quintuple_target_graph.md` | 14 | 252 | needs content map | medium | Stage 252 already gave an exact realization compiler: |
| `moving_throat_pde_stage202_free_quintuple_target_graph.md` | 497 | 252 | needs content map | medium | Stage 252 already told us how to project any state onto |
| `moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md` | 19 | 249 | needs content map | medium | Stage 249 already showed that the explicit `5`PN suppor |
| `moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md` | 82 | 239 | needs content map | medium | Stage 239 already gives the first-order observable comp |
| `moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md` | 183 | 241 | needs content map | medium | Stage 241 already transported the primitive-family same |
| `moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md` | 303 | 249 | needs content map | medium | Stage 249 already showed that the explicit `5`PN suppor |
| `moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md` | 55 | 241 | needs content map | medium | At the base point `\(\boldsymbol\ell_\circ\)`, Stage 24 |
| `moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md` | 469 | 241 | needs content map | medium | Stage 241 closed the primitive certified ray table and  |
| `moving_throat_pde_stage184_branch_invariant_coordinates.md` | 188 | 251 | needs content map | medium | Nothing essentially new had to be added. Stage 251 alre |
| `moving_throat_pde_stage184_branch_invariant_coordinates.md` | 334 | 251 | needs content map | medium | Since Stage 251 already proved |
| `moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md` | 15 | 238 | needs content map | medium | Stage 238 reduced the linear dynamic corridor to one ex |
| `moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md` | 36 | 239 | needs content map | medium | So Stage 239 keeps the dynamic same-charge idea alive,  |
| `moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md` | 48 | 238 | needs content map | medium | Stage 238 already showed that the only linear corridor  |
| `moving_throat_pde_stage165_exact_branch_drifts.md` | 5 | 249 | needs content map | medium | Stage 249 reduced the first off-family defect to the ex |
| `moving_throat_pde_stage188_branch_observables_completion.md` | 216 | 238 | needs content map | medium | Stage 238 already fixed the tangent quotient packet as |
| `moving_throat_pde_stage188_branch_observables_completion.md` | 457 | 238 | needs content map | medium | Stage 238 already reduced the coherent weak-axisymmetri |
| `moving_throat_pde_stage181_coherent_tracking_defect.md` | 6 | 248 | needs content map | medium | Stage 248 reduced the remaining grouped weak-axisymmetr |
| `moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md` | 165 | 247 | needs content map | medium | So Stage 247 already gave two exact global statements: |
| `moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md` | 529 | 247 | needs content map | medium | Stage 247 already proved the global inequalities |
| `moving_throat_pde_stage176_outgoing_load_factorization.md` | 5 | 243 | needs content map | medium | Stage 243 reduced the remaining linear grouped weak-axi |
| `moving_throat_pde_stage176_outgoing_load_factorization.md` | 178 | 243 | needs content map | medium | Stage 243 already showed that if the conservative shape |
| `moving_throat_pde_stage176_outgoing_load_factorization.md` | 274 | 243 | needs content map | medium | Stage 243 reduced the remaining linear grouped defect t |
| `moving_throat_pde_stage185_microscopic_monomials.md` | 6 | 252 | needs content map | medium | Stage 252 reduced the coherent weak-axisymmetric contin |
| `moving_throat_pde_stage185_microscopic_monomials.md` | 92 | 250 | needs content map | medium | Stage 250 already gave the microscopic drift variables |
| `moving_throat_pde_stage185_microscopic_monomials.md` | 322 | 251 | needs content map | medium | Stage 251 already proved the exact triangular observabl |
| `moving_throat_pde_stage185_microscopic_monomials.md` | 462 | 252 | needs content map | medium | Stage 252 reduced the continuation point to the drifts  |
| `moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md` | 14 | 238 | needs content map | medium | Stage 238 already gave the exact primitive free-directi |
| `moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md` | 206 | 238 | needs content map | medium | Stage 238 already fixed the dependent graph exponents i |
| `moving_throat_pde_stage179_transfer_shape_theorem.md` | 5 | 246 | needs content map | medium | Stage 246 reduced the remaining weak-axisymmetric group |
| `moving_throat_pde_stage179_transfer_shape_theorem.md` | 238 | 246 | needs content map | medium | Stage 246 already gave |
| `moving_throat_pde_stage179_transfer_shape_theorem.md` | 425 | 246 | needs content map | medium | Stage 246 already reduced the remaining grouped weak-ax |
| `moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md` | 14 | 238 | needs content map | medium | Stage 238 reduced the graph-lifted home-stretch problem |
| `moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md` | 70 | 238 | needs content map | medium | and Stage 238 already proved that this graph lift lies  |
| `moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md` | 395 | 238 | needs content map | medium | Stage 238 reduced the reduced closure search to a scala |
| `moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md` | 533 | 253 | needs content map | medium | Stage 253 reduced the target orbit from an abstract sim |
| `moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md` | 21 | 245 | needs content map | medium | Stage 245 already reduced the surviving same-charge pur |
| `moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md` | 561 | 247 | needs content map | medium | Stage 247 reduced the first support-cardinality-`4` aud |
| `moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md` | 584 | 247 | needs content map | medium | - Stage 247 closed the boundary and canonical-screen au |
| `moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md` | 518 | 252 | needs content map | medium | - Stage 252 closed the local mixed-ray search sieve. |
| `moving_throat_pde_stage178_outgoing_port_coloading_theorem.md` | 5 | 245 | needs content map | medium | Stage 245 reduced the remaining weak-axisymmetric group |
| `moving_throat_pde_stage178_outgoing_port_coloading_theorem.md` | 103 | 243 | needs content map | medium | Stage 243 already showed that on conservative-shape-pre |
| `moving_throat_pde_stage178_outgoing_port_coloading_theorem.md` | 410 | 245 | needs content map | medium | Stage 245 already showed that the remaining weak-axisym |
| `moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md` | 248 | 245 | needs content map | medium | Stage 245 already showed that the most general first is |

**Special Class-D items from the residual cross-check:** `mathematica/...stage227_..._mathematica_audit.wl:154,243,244` say `Stage-242 K compatibility` / `M4 Stage-243 transfer basis` -- stems for those cross-refs map to canonical 225/226 per the P4-53 residual note; needs Codex content confirmation (the .wl banner itself is already canonical STAGE 227).

## Logged-residual cross-check (P4-53 / P4-54)

### P4-53 (batch VII.1, stages 219-230)
- **(a) VII.1 notes multi-epoch self-number drift** (log: 219 self-numbers 253; 221/222/225 self-number 238/239/242; 221's notes call today's 219 'Stage 253'): **CONFIRMED PRESENT, RECLASSIFIED as Class D (body prose), not title.** The H1 TITLES of 219/221/222/225 are already CANONICAL; stale self-numbers live in BODY prose (stage219 lines 14/22 'Stage 252', line 41 'Stage 253'; stage222 lines 15/36 'Stage 238/239'). Captured in Class D. Each note's appendix FOOTER also carries a stale `stageMMM_<stem>` script citation -> captured in **Class A** (219->cites 253, 222->239, etc.).
- **(b) cosmetic .wl labels:** stage221 .wl:35 banner `STAGE 204` (canon 221) -- **CONFIRMED PRESENT** (Class C). stage227 .wl:154/243/244 `Stage-242/Stage-243` -- **CONFIRMED PRESENT** (Class D cross-ref).
- **RENUMBER findings the tracker says were APPLIED** (221: 238->221/237->220; 225: 240/241/242->223/224/225; informational 220/224/226/230): banners for 220/224/226/230/233/239/242 **CONFIRMED already canonical** (fixed). Residual = notes-body + the two .wl cosmetics above.

### P4-54 (batch VII.2, stages 231-242)
- **Notes-title drift 232 'Stage 249', 234 'Stage 251', 235 'Stage 251/252', 236:** stage232 .md H1 = 'Stage 249' -- **CONFIRMED PRESENT** (Class B). stage234/235/236 .md H1 lines are blank-padded; their stale self-numbers surface in the appendix-footer `stageMMM_<stem>` citations (234->251, 235->252, 236->253) -- **CONFIRMED PRESENT** (Class A). H1-line title drift for 234/235/236 was NOT found on the title line -- residual is real but sits in Class A footers + Class D body; **slightly EXPANDED** vs the 'notes-title' wording.
- **In-file .wl cosmetics the tracker says were APPLIED** (233: 188/223/224->239/240/241; 239: STAGE 222->239 + Stage221->238; 242: STAGE 225->242): **CONFIRMED already canonical/fixed** (233/239/242 banners read 233/239/242).

**Net:** every logged residual is CONFIRMED present (or confirmed-fixed where the log said APPLIED). The scan EXPANDS them -- the VII.1/VII.2 drift is not only titles/banners but also Class-A footer citations and Class-D body self-references, and the same pattern recurs project-wide.

## SUMMARY

| class | description | count | confidence |
| --- | --- | --- | --- |
| A | stale `stageNNN_<stem>` citations | 25 | HIGH (stem-keyed) |
| B | notes self-title / H1 drift | 19 | HIGH (own canonical) |
| C.1 | script self-banner drift (.py/.wl, 122 scripts) | 155 | HIGH (own canonical) |
| C.2 | paper-card self-title drift (091-124) | 34 | HIGH (label/filename) |
| D | body prose cross-refs (high-signal flagged) | 121 | MEDIUM/LOW (content map) |

- **Total tabulated stale-label findings: 354** (plus a looser ~1,200-mention Class-D pool not individually tabulated).
- **HIGH-confidence, deterministic, safe for stem/file-keyed Codex fix: 233** (A+B+C.1+C.2 = 25+19+155+34).
- **Needs content review (Class D): 121 flagged candidates** + the broad ~1,200-mention pool.
- **Stage-range spread: ~stage 021 -> 242 (canonical) -- PROJECT-WIDE**, not confined to VII.1/VII.2.

### Surprises / cautions
1. **Published paper cards are NOT clean for their own titles.** Cards stage_091..stage_124 carry section titles offset +17 from their label/filename (e.g. stage_100.tex titled 'Stage~117'). Labels/filenames ARE canonical; visible titles are stale -- affects the rendered paper, in-scope.
2. **Within the +17 band, fixes are still NOT uniform per file.** In stage_091/stage_100.tex the body prose is already canonical while only the title is stale; in stage_124.tex the body 'Stage~141' is ALSO stale. Per-occurrence, not per-file.
3. **Offsets are genuinely inconsistent** (+17 dominant, also +1,+2,+22,+40,+51,+68,+215,-31,-67), confirming a uniform sweep is unsafe.
4. **False-positive class to avoid:** longer stems extending a shorter canonical stem (e.g. `stage132_mouth_boundary_layer_status` is a REAL tex-only stage; do NOT 'correct' it to 129's `mouth_boundary_layer`). The Class-A scan filtered these by checking the cited file exists.
5. **Class A/B/C are mechanically resolvable today; Class D is the real residual** motivating the user's stem-lookup post-253 cleanup -- needs Codex-applied + Claude-reviewed content checks.
