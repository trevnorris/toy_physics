# Phase A Claim-Label Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `001, 002, 003, 004, 005, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 023, 024, 025, 026, 027, 028, 029, 030, 031, 032, 033, 034, 035, 036, 037, 038, 039, 040, 041, 042, 043, 044, 045, 046, 047, 048, 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 082, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 093, 094, 095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253`
- Source files: `- stage: '001'
  role: paper_stage_tex
  path: paper/stages/stage_001.tex
- stage: '001'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage001_geometry_lift.md
- stage: '001'
  role: sympy_script
  path: scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py
- stage: '001'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl
- stage: '002'
  role: paper_stage_tex
  path: paper/stages/stage_002.tex
- stage: '002'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage002_breathing_reduction.md
- stage: '002'
  role: sympy_script
  path: scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py
- stage: '002'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl
- stage: '003'
  role: paper_stage_tex
  path: paper/stages/stage_003.tex
- stage: '003'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage003_bdg_coupling.md
- stage: '003'
  role: sympy_script
  path: scripts/moving_throat_pde_stage003_bdg_sympy_audit.py
- stage: '003'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl
- stage: '004'
  role: paper_stage_tex
  path: paper/stages/stage_004.tex
- stage: '004'
  role: sympy_script
  path: scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py
- stage: '004'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl
- stage: '005'
  role: paper_stage_tex
  path: paper/stages/stage_005.tex
- stage: '005'
  role: sympy_script
  path: scripts/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py
- stage: '005'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl
- stage: '006'
  role: paper_stage_tex
  path: paper/stages/stage_006.tex
- stage: '006'
  role: sympy_script
  path: scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py
- stage: '006'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl
- stage: '007'
  role: paper_stage_tex
  path: paper/stages/stage_007.tex
- stage: '007'
  role: sympy_script
  path: scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py
- stage: '007'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl
- stage: 008
  role: paper_stage_tex
  path: paper/stages/stage_008.tex
- stage: 008
  role: sympy_script
  path: scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py
- stage: 008
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl
- stage: 009
  role: paper_stage_tex
  path: paper/stages/stage_009.tex
- stage: 009
  role: sympy_script
  path: scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py
- stage: 009
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl
- stage: '010'
  role: paper_stage_tex
  path: paper/stages/stage_010.tex
- stage: '010'
  role: sympy_script
  path: scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py
- stage: '010'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl
- stage: '011'
  role: paper_stage_tex
  path: paper/stages/stage_011.tex
- stage: '011'
  role: sympy_script
  path: scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py
- stage: '011'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl
- stage: '012'
  role: paper_stage_tex
  path: paper/stages/stage_012.tex
- stage: '012'
  role: sympy_script
  path: scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py
- stage: '012'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl
- stage: '013'
  role: paper_stage_tex
  path: paper/stages/stage_013.tex
- stage: '013'
  role: sympy_script
  path: scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py
- stage: '013'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl
- stage: '014'
  role: paper_stage_tex
  path: paper/stages/stage_014.tex
- stage: '014'
  role: sympy_script
  path: scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py
- stage: '014'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl
- stage: '015'
  role: paper_stage_tex
  path: paper/stages/stage_015.tex
- stage: '015'
  role: sympy_script
  path: scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py
- stage: '015'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl
- stage: '016'
  role: paper_stage_tex
  path: paper/stages/stage_016.tex
- stage: '016'
  role: sympy_script
  path: scripts/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.py
- stage: '016'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl
- stage: '017'
  role: paper_stage_tex
  path: paper/stages/stage_017.tex
- stage: '017'
  role: sympy_script
  path: scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py
- stage: '017'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl
- stage: 018
  role: paper_stage_tex
  path: paper/stages/stage_018.tex
- stage: 018
  role: sympy_script
  path: scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py
- stage: 018
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl
- stage: 019
  role: paper_stage_tex
  path: paper/stages/stage_019.tex
- stage: 019
  role: sympy_script
  path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
- stage: 019
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl
- stage: '020'
  role: paper_stage_tex
  path: paper/stages/stage_020.tex
- stage: '020'
  role: sympy_script
  path: scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py
- stage: '020'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl
- stage: '021'
  role: paper_stage_tex
  path: paper/stages/stage_021.tex
- stage: '021'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md
- stage: '021'
  role: sympy_script
  path: scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py
- stage: '021'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl
- stage: '022'
  role: paper_stage_tex
  path: paper/stages/stage_022.tex
- stage: '022'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
- stage: '022'
  role: sympy_script
  path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
- stage: '022'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl
- stage: '023'
  role: paper_stage_tex
  path: paper/stages/stage_023.tex
- stage: '023'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md
- stage: '023'
  role: sympy_script
  path: scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py
- stage: '023'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl
- stage: '024'
  role: paper_stage_tex
  path: paper/stages/stage_024.tex
- stage: '024'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage024_overlap_isotropy.md
- stage: '024'
  role: sympy_script
  path: scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py
- stage: '024'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl
- stage: '025'
  role: paper_stage_tex
  path: paper/stages/stage_025.tex
- stage: '025'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
- stage: '025'
  role: sympy_script
  path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
- stage: '025'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl
- stage: '026'
  role: paper_stage_tex
  path: paper/stages/stage_026.tex
- stage: '026'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
- stage: '026'
  role: sympy_script
  path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
- stage: '026'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl
- stage: '027'
  role: paper_stage_tex
  path: paper/stages/stage_027.tex
- stage: '027'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage027_nonconstant_axial_family.md
- stage: '027'
  role: sympy_script
  path: scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py
- stage: '027'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl
- stage: 028
  role: paper_stage_tex
  path: paper/stages/stage_028.tex
- stage: 028
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage028_loaded_profile_selection.md
- stage: 028
  role: sympy_script
  path: scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py
- stage: 028
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl
- stage: 029
  role: paper_stage_tex
  path: paper/stages/stage_029.tex
- stage: 029
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage029_dynamic_loading.md
- stage: 029
  role: sympy_script
  path: scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py
- stage: 029
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl
- stage: '030'
  role: paper_stage_tex
  path: paper/stages/stage_030.tex
- stage: '030'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
- stage: '030'
  role: sympy_script
  path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
- stage: '030'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl
- stage: '031'
  role: paper_stage_tex
  path: paper/stages/stage_031.tex
- stage: '031'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
- stage: '031'
  role: sympy_script
  path: scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py
- stage: '031'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl
- stage: '032'
  role: paper_stage_tex
  path: paper/stages/stage_032.tex
- stage: '032'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage032_source_map_from_mode_integrals.md
- stage: '032'
  role: sympy_script
  path: scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py
- stage: '032'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl
- stage: '033'
  role: paper_stage_tex
  path: paper/stages/stage_033.tex
- stage: '033'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage033_microscopic_normalization_equation.md
- stage: '033'
  role: sympy_script
  path: scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py
- stage: '033'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl
- stage: '034'
  role: paper_stage_tex
  path: paper/stages/stage_034.tex
- stage: '034'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage034_softening_depth_normal_form.md
- stage: '034'
  role: sympy_script
  path: scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py
- stage: '034'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl
- stage: '035'
  role: paper_stage_tex
  path: paper/stages/stage_035.tex
- stage: '035'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
- stage: '035'
  role: sympy_script
  path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
- stage: '035'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl
- stage: '036'
  role: paper_stage_tex
  path: paper/stages/stage_036.tex
- stage: '036'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage036_support_feasibility_frontier.md
- stage: '036'
  role: sympy_script
  path: scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py
- stage: '036'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl
- stage: '037'
  role: paper_stage_tex
  path: paper/stages/stage_037.tex
- stage: '037'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md
- stage: '037'
  role: sympy_script
  path: scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py
- stage: '037'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl
- stage: 038
  role: paper_stage_tex
  path: paper/stages/stage_038.tex
- stage: 038
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
- stage: 038
  role: sympy_script
  path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
- stage: 038
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl
- stage: 039
  role: paper_stage_tex
  path: paper/stages/stage_039.tex
- stage: 039
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage039_split_u_sector.md
- stage: 039
  role: sympy_script
  path: scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py
- stage: 039
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl
- stage: '040'
  role: paper_stage_tex
  path: paper/stages/stage_040.tex
- stage: '040'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage040_generalized_selected_branch.md
- stage: '040'
  role: sympy_script
  path: scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py
- stage: '040'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl
- stage: '041'
  role: paper_stage_tex
  path: paper/stages/stage_041.tex
- stage: '041'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage041_rank2_support_completion.md
- stage: '041'
  role: sympy_script
  path: scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py
- stage: '041'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl
- stage: '042'
  role: paper_stage_tex
  path: paper/stages/stage_042.tex
- stage: '042'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage042_rank2_selected_mode_normalization.md
- stage: '042'
  role: sympy_script
  path: scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py
- stage: '042'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl
- stage: '043'
  role: paper_stage_tex
  path: paper/stages/stage_043.tex
- stage: '043'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage043_support_direction_extraction.md
- stage: '043'
  role: sympy_script
  path: scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py
- stage: '043'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl
- stage: '044'
  role: paper_stage_tex
  path: paper/stages/stage_044.tex
- stage: '044'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md
- stage: '044'
  role: sympy_script
  path: scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py
- stage: '044'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl
- stage: '045'
  role: paper_stage_tex
  path: paper/stages/stage_045.tex
- stage: '045'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md
- stage: '045'
  role: sympy_script
  path: scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py
- stage: '045'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl
- stage: '046'
  role: paper_stage_tex
  path: paper/stages/stage_046.tex
- stage: '046'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md
- stage: '046'
  role: sympy_script
  path: scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py
- stage: '046'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl
- stage: '047'
  role: paper_stage_tex
  path: paper/stages/stage_047.tex
- stage: '047'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage047_coherent_kernel_map.md
- stage: '047'
  role: sympy_script
  path: scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py
- stage: '047'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl
- stage: 048
  role: paper_stage_tex
  path: paper/stages/stage_048.tex
- stage: 048
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage048_support_compensation_theorem.md
- stage: 048
  role: sympy_script
  path: scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py
- stage: 048
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl
- stage: 049
  role: paper_stage_tex
  path: paper/stages/stage_049.tex
- stage: 049
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage049_dn_overlap_zeta.md
- stage: 049
  role: sympy_script
  path: scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py
- stage: 049
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl
- stage: '050'
  role: paper_stage_tex
  path: paper/stages/stage_050.tex
- stage: '050'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md
- stage: '050'
  role: sympy_script
  path: scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py
- stage: '050'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl
- stage: '051'
  role: paper_stage_tex
  path: paper/stages/stage_051.tex
- stage: '051'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage051_lowest_twin_criterion.md
- stage: '051'
  role: sympy_script
  path: scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py
- stage: '051'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl
- stage: '052'
  role: paper_stage_tex
  path: paper/stages/stage_052.tex
- stage: '052'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage052_nontwin_asymmetry_threshold.md
- stage: '052'
  role: sympy_script
  path: scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py
- stage: '052'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl
- stage: '053'
  role: paper_stage_tex
  path: paper/stages/stage_053.tex
- stage: '053'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage053_overlap_boost_window.md
- stage: '053'
  role: sympy_script
  path: scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py
- stage: '053'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl
- stage: '054'
  role: paper_stage_tex
  path: paper/stages/stage_054.tex
- stage: '054'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage054_robin_softening_support_lane.md
- stage: '054'
  role: sympy_script
  path: scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py
- stage: '054'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl
- stage: '055'
  role: paper_stage_tex
  path: paper/stages/stage_055.tex
- stage: '055'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage055_explicit_lowest_lane_reachability.md
- stage: '055'
  role: sympy_script
  path: scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py
- stage: '055'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl
- stage: '056'
  role: paper_stage_tex
  path: paper/stages/stage_056.tex
- stage: '056'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md
- stage: '056'
  role: sympy_script
  path: scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py
- stage: '056'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl
- stage: '057'
  role: paper_stage_tex
  path: paper/stages/stage_057.tex
- stage: '057'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage057_physical_parameter_map.md
- stage: '057'
  role: sympy_script
  path: scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py
- stage: '057'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl
- stage: 058
  role: paper_stage_tex
  path: paper/stages/stage_058.tex
- stage: 058
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage058_coupled_support_source_operator.md
- stage: 058
  role: sympy_script
  path: scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py
- stage: 058
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl
- stage: 059
  role: paper_stage_tex
  path: paper/stages/stage_059.tex
- stage: 059
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage059_operator_branch_residual_bounds.md
- stage: 059
  role: sympy_script
  path: scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py
- stage: 059
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl
- stage: '060'
  role: paper_stage_tex
  path: paper/stages/stage_060.tex
- stage: '060'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage060_entropic_microclosure.md
- stage: '060'
  role: sympy_script
  path: scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py
- stage: '060'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl
- stage: '061'
  role: paper_stage_tex
  path: paper/stages/stage_061.tex
- stage: '061'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage061_microscopic_gain_thresholds.md
- stage: '061'
  role: sympy_script
  path: scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py
- stage: '061'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl
- stage: '062'
  role: paper_stage_tex
  path: paper/stages/stage_062.tex
- stage: '062'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage062_parent_action_gain.md
- stage: '062'
  role: sympy_script
  path: scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py
- stage: '062'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl
- stage: '063'
  role: paper_stage_tex
  path: paper/stages/stage_063.tex
- stage: '063'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage063_parent_thresholds.md
- stage: '063'
  role: sympy_script
  path: scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py
- stage: '063'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl
- stage: '064'
  role: paper_stage_tex
  path: paper/stages/stage_064.tex
- stage: '064'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage064_equilibrium_alignment.md
- stage: '064'
  role: sympy_script
  path: scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py
- stage: '064'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl
- stage: '065'
  role: paper_stage_tex
  path: paper/stages/stage_065.tex
- stage: '065'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage065_thin_wall_confinement.md
- stage: '065'
  role: sympy_script
  path: scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py
- stage: '065'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl
- stage: '066'
  role: paper_stage_tex
  path: paper/stages/stage_066.tex
- stage: '066'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage066_wall_figure_of_merit.md
- stage: '066'
  role: sympy_script
  path: scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py
- stage: '066'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl
- stage: '067'
  role: paper_stage_tex
  path: paper/stages/stage_067.tex
- stage: '067'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md
- stage: '067'
  role: sympy_script
  path: scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py
- stage: '067'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl
- stage: 068
  role: paper_stage_tex
  path: paper/stages/stage_068.tex
- stage: 068
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage068_resonance_thresholds.md
- stage: 068
  role: sympy_script
  path: scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py
- stage: 068
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl
- stage: 069
  role: paper_stage_tex
  path: paper/stages/stage_069.tex
- stage: 069
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage069_final_reduced_verdict.md
- stage: 069
  role: sympy_script
  path: scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py
- stage: 069
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl
- stage: '070'
  role: paper_stage_tex
  path: paper/stages/stage_070.tex
- stage: '070'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage070_gnls_wall_shell.md
- stage: '070'
  role: sympy_script
  path: scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py
- stage: '070'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl
- stage: '071'
  role: paper_stage_tex
  path: paper/stages/stage_071.tex
- stage: '071'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage071_tanh_wall_branch.md
- stage: '071'
  role: sympy_script
  path: scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py
- stage: '071'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl
- stage: '072'
  role: paper_stage_tex
  path: paper/stages/stage_072.tex
- stage: '072'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md
- stage: '072'
  role: sympy_script
  path: scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py
- stage: '072'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl
- stage: '073'
  role: paper_stage_tex
  path: paper/stages/stage_073.tex
- stage: '073'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage073_family1_geometry_map.md
- stage: '073'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl
- stage: '074'
  role: paper_stage_tex
  path: paper/stages/stage_074.tex
- stage: '074'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage074_family1_healing_lock.md
- stage: '074'
  role: sympy_script
  path: scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py
- stage: '074'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl
- stage: '075'
  role: paper_stage_tex
  path: paper/stages/stage_075.tex
- stage: '075'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage075_family1_threshold_window.md
- stage: '075'
  role: sympy_script
  path: scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py
- stage: '075'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl
- stage: '076'
  role: paper_stage_tex
  path: paper/stages/stage_076.tex
- stage: '076'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md
- stage: '076'
  role: sympy_script
  path: scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py
- stage: '076'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl
- stage: '077'
  role: paper_stage_tex
  path: paper/stages/stage_077.tex
- stage: '077'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage077_family1_theta_extraction.md
- stage: '077'
  role: sympy_script
  path: scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py
- stage: '077'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl
- stage: 078
  role: paper_stage_tex
  path: paper/stages/stage_078.tex
- stage: 078
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md
- stage: 078
  role: sympy_script
  path: scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py
- stage: 078
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl
- stage: 079
  role: paper_stage_tex
  path: paper/stages/stage_079.tex
- stage: 079
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md
- stage: 079
  role: sympy_script
  path: scripts/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.py
- stage: 079
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl
- stage: 080
  role: paper_stage_tex
  path: paper/stages/stage_080.tex
- stage: 080
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md
- stage: 080
  role: sympy_script
  path: scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py
- stage: 080
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl
- stage: 081
  role: paper_stage_tex
  path: paper/stages/stage_081.tex
- stage: 081
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md
- stage: 081
  role: sympy_script
  path: scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py
- stage: 081
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl
- stage: 082
  role: paper_stage_tex
  path: paper/stages/stage_082.tex
- stage: 082
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md
- stage: 082
  role: sympy_script
  path: scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py
- stage: 082
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl
- stage: 083
  role: paper_stage_tex
  path: paper/stages/stage_083.tex
- stage: 083
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md
- stage: 083
  role: sympy_script
  path: scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py
- stage: 083
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl
- stage: 084
  role: paper_stage_tex
  path: paper/stages/stage_084.tex
- stage: 084
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md
- stage: 084
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl
- stage: 085
  role: paper_stage_tex
  path: paper/stages/stage_085.tex
- stage: 085
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage085_quadrupole_demand_cancellation.md
- stage: 085
  role: sympy_script
  path: scripts/moving_throat_pde_stage085_quadrupole_demand_cancellation_sympy_audit.py
- stage: 085
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage085_quadrupole_demand_cancellation_mathematica_audit.wl
- stage: 086
  role: paper_stage_tex
  path: paper/stages/stage_086.tex
- stage: 086
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage086_family1_loading_ratio_window.md
- stage: 086
  role: sympy_script
  path: scripts/moving_throat_pde_stage086_family1_loading_ratio_window_sympy_audit.py
- stage: 086
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_mathematica_audit.wl
- stage: 087
  role: paper_stage_tex
  path: paper/stages/stage_087.tex
- stage: 087
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md
- stage: 087
  role: sympy_script
  path: scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py
- stage: 087
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl
- stage: 088
  role: paper_stage_tex
  path: paper/stages/stage_088.tex
- stage: 088
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage088_loading_ratio_from_minimal_module.md
- stage: 088
  role: sympy_script
  path: scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py
- stage: 088
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl
- stage: 089
  role: paper_stage_tex
  path: paper/stages/stage_089.tex
- stage: 089
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md
- stage: 089
  role: sympy_script
  path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
- stage: 089
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl
- stage: 090
  role: paper_stage_tex
  path: paper/stages/stage_090.tex
- stage: 090
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage090_updated_reduced_status.md
- stage: 090
  role: sympy_script
  path: scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py
- stage: 090
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl
- stage: 091
  role: paper_stage_tex
  path: paper/stages/stage_091.tex
- stage: 091
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md
- stage: 091
  role: sympy_script
  path: scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py
- stage: 091
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl
- stage: 092
  role: paper_stage_tex
  path: paper/stages/stage_092.tex
- stage: 092
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage092_dynamic_geometry_obstruction.md
- stage: 092
  role: sympy_script
  path: scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py
- stage: 092
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl
- stage: 093
  role: paper_stage_tex
  path: paper/stages/stage_093.tex
- stage: 093
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage093_grouped_p2_status_update.md
- stage: 093
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl
- stage: 094
  role: paper_stage_tex
  path: paper/stages/stage_094.tex
- stage: 094
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage094_isotropic_geometry_decoupling.md
- stage: 094
  role: sympy_script
  path: scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py
- stage: 094
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl
- stage: 095
  role: paper_stage_tex
  path: paper/stages/stage_095.tex
- stage: 095
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md
- stage: 095
  role: sympy_script
  path: scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py
- stage: 095
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl
- stage: 096
  role: paper_stage_tex
  path: paper/stages/stage_096.tex
- stage: 096
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage096_geometry_lane_check_verdict.md
- stage: 096
  role: sympy_script
  path: scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py
- stage: 096
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl
- stage: 097
  role: paper_stage_tex
  path: paper/stages/stage_097.tex
- stage: 097
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage097_single_normalization_defect.md
- stage: 097
  role: sympy_script
  path: scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py
- stage: 097
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl
- stage: 098
  role: paper_stage_tex
  path: paper/stages/stage_098.tex
- stage: 098
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md
- stage: 098
  role: sympy_script
  path: scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py
- stage: 098
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl
- stage: 099
  role: paper_stage_tex
  path: paper/stages/stage_099.tex
- stage: 099
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage099_reduced_finish_line.md
- stage: 099
  role: sympy_script
  path: scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py
- stage: 099
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl
- stage: '100'
  role: paper_stage_tex
  path: paper/stages/stage_100.tex
- stage: '100'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md
- stage: '100'
  role: sympy_script
  path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
- stage: '100'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl
- stage: '101'
  role: paper_stage_tex
  path: paper/stages/stage_101.tex
- stage: '101'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage101_natural_source_map_reduction.md
- stage: '101'
  role: sympy_script
  path: scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py
- stage: '101'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl
- stage: '102'
  role: paper_stage_tex
  path: paper/stages/stage_102.tex
- stage: '102'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage102_higher_odd_irrelevance.md
- stage: '102'
  role: sympy_script
  path: scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py
- stage: '102'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl
- stage: '103'
  role: paper_stage_tex
  path: paper/stages/stage_103.tex
- stage: '103'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage103_reduced_25pn_conditional_closure.md
- stage: '104'
  role: paper_stage_tex
  path: paper/stages/stage_104.tex
- stage: '104'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
- stage: '104'
  role: sympy_script
  path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
- stage: '104'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl
- stage: '105'
  role: paper_stage_tex
  path: paper/stages/stage_105.tex
- stage: '105'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
- stage: '105'
  role: sympy_script
  path: scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py
- stage: '105'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl
- stage: '106'
  role: paper_stage_tex
  path: paper/stages/stage_106.tex
- stage: '106'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
- stage: '106'
  role: sympy_script
  path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
- stage: '106'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl
- stage: '107'
  role: paper_stage_tex
  path: paper/stages/stage_107.tex
- stage: '107'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md
- stage: '107'
  role: sympy_script
  path: scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py
- stage: '107'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl
- stage: '108'
  role: paper_stage_tex
  path: paper/stages/stage_108.tex
- stage: '108'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage108_robustness_classes.md
- stage: '108'
  role: sympy_script
  path: scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py
- stage: '108'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl
- stage: '109'
  role: paper_stage_tex
  path: paper/stages/stage_109.tex
- stage: '109'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage109_linearized_branch_selection.md
- stage: '109'
  role: sympy_script
  path: scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py
- stage: '109'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl
- stage: '110'
  role: paper_stage_tex
  path: paper/stages/stage_110.tex
- stage: '110'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage110_robin_outlet_model.md
- stage: '110'
  role: sympy_script
  path: scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py
- stage: '110'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl
- stage: '111'
  role: paper_stage_tex
  path: paper/stages/stage_111.tex
- stage: '111'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage111_mixed_sidechannel_pole.md
- stage: '111'
  role: sympy_script
  path: scripts/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.py
- stage: '111'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl
- stage: '112'
  role: paper_stage_tex
  path: paper/stages/stage_112.tex
- stage: '112'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md
- stage: '112'
  role: sympy_script
  path: scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py
- stage: '112'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl
- stage: '113'
  role: paper_stage_tex
  path: paper/stages/stage_113.tex
- stage: '113'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage113_outlet_model_status.md
- stage: '114'
  role: paper_stage_tex
  path: paper/stages/stage_114.tex
- stage: '114'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage114_concrete_core_schur.md
- stage: '114'
  role: sympy_script
  path: scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py
- stage: '114'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl
- stage: '115'
  role: paper_stage_tex
  path: paper/stages/stage_115.tex
- stage: '115'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage115_core_balance_compensation.md
- stage: '115'
  role: sympy_script
  path: scripts/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py
- stage: '115'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl
- stage: '116'
  role: paper_stage_tex
  path: paper/stages/stage_116.tex
- stage: '116'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md
- stage: '116'
  role: sympy_script
  path: scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py
- stage: '116'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl
- stage: '117'
  role: paper_stage_tex
  path: paper/stages/stage_117.tex
- stage: '117'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage117_outlet_core_status.md
- stage: '117'
  role: sympy_script
  path: scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py
- stage: '117'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl
- stage: '118'
  role: paper_stage_tex
  path: paper/stages/stage_118.tex
- stage: '118'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage118_parent_core_extraction.md
- stage: '118'
  role: sympy_script
  path: scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py
- stage: '118'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl
- stage: '119'
  role: paper_stage_tex
  path: paper/stages/stage_119.tex
- stage: '119'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage119_parent_balance_family.md
- stage: '119'
  role: sympy_script
  path: scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py
- stage: '119'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl
- stage: '120'
  role: paper_stage_tex
  path: paper/stages/stage_120.tex
- stage: '120'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage120_core_parameter_status.md
- stage: '121'
  role: paper_stage_tex
  path: paper/stages/stage_121.tex
- stage: '121'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage121_geometric_r_selection.md
- stage: '121'
  role: sympy_script
  path: scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py
- stage: '121'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl
- stage: '122'
  role: paper_stage_tex
  path: paper/stages/stage_122.tex
- stage: '122'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md
- stage: '122'
  role: sympy_script
  path: scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py
- stage: '122'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl
- stage: '123'
  role: paper_stage_tex
  path: paper/stages/stage_123.tex
- stage: '123'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md
- stage: '123'
  role: sympy_script
  path: scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py
- stage: '123'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl
- stage: '124'
  role: paper_stage_tex
  path: paper/stages/stage_124.tex
- stage: '124'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage124_core_branch_status.md
- stage: '125'
  role: paper_stage_tex
  path: paper/stages/stage_125.tex
- stage: '125'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage125_positive_source_theorem.md
- stage: '125'
  role: sympy_script
  path: scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py
- stage: '125'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl
- stage: '126'
  role: paper_stage_tex
  path: paper/stages/stage_126.tex
- stage: '126'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage126_positive_source_families.md
- stage: '126'
  role: sympy_script
  path: scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py
- stage: '126'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl
- stage: '127'
  role: paper_stage_tex
  path: paper/stages/stage_127.tex
- stage: '127'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage127_penetration_families.md
- stage: '127'
  role: sympy_script
  path: scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py
- stage: '127'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl
- stage: '128'
  role: paper_stage_tex
  path: paper/stages/stage_128.tex
- stage: '128'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage128_mouth_source_bias_status.md
- stage: '129'
  role: paper_stage_tex
  path: paper/stages/stage_129.tex
- stage: '129'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md
- stage: '129'
  role: sympy_script
  path: scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py
- stage: '129'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl
- stage: '130'
  role: paper_stage_tex
  path: paper/stages/stage_130.tex
- stage: '130'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage130_mouth_bias_map.md
- stage: '130'
  role: sympy_script
  path: scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py
- stage: '130'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl
- stage: '131'
  role: paper_stage_tex
  path: paper/stages/stage_131.tex
- stage: '131'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md
- stage: '131'
  role: sympy_script
  path: scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py
- stage: '131'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl
- stage: '132'
  role: paper_stage_tex
  path: paper/stages/stage_132.tex
- stage: '132'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md
- stage: '133'
  role: paper_stage_tex
  path: paper/stages/stage_133.tex
- stage: '133'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage133_coupled_mouth_fixedpoint.md
- stage: '133'
  role: sympy_script
  path: scripts/moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py
- stage: '133'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl
- stage: '134'
  role: paper_stage_tex
  path: paper/stages/stage_134.tex
- stage: '134'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md
- stage: '134'
  role: sympy_script
  path: scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py
- stage: '134'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl
- stage: '135'
  role: paper_stage_tex
  path: paper/stages/stage_135.tex
- stage: '135'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage135_outlet_consistent_mouth_closure.md
- stage: '135'
  role: sympy_script
  path: scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py
- stage: '135'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl
- stage: '136'
  role: paper_stage_tex
  path: paper/stages/stage_136.tex
- stage: '136'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md
- stage: '137'
  role: paper_stage_tex
  path: paper/stages/stage_137.tex
- stage: '137'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md
- stage: '137'
  role: sympy_script
  path: scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py
- stage: '137'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl
- stage: '138'
  role: paper_stage_tex
  path: paper/stages/stage_138.tex
- stage: '138'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md
- stage: '138'
  role: sympy_script
  path: scripts/moving_throat_pde_stage138_normalized_mouth_gain_family_sympy_audit.py
- stage: '138'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl
- stage: '139'
  role: paper_stage_tex
  path: paper/stages/stage_139.tex
- stage: '139'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md
- stage: '139'
  role: sympy_script
  path: scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py
- stage: '139'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl
- stage: '140'
  role: paper_stage_tex
  path: paper/stages/stage_140.tex
- stage: '140'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md
- stage: '140'
  role: sympy_script
  path: scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py
- stage: '140'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl
- stage: '141'
  role: paper_stage_tex
  path: paper/stages/stage_141.tex
- stage: '141'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage141_mouth_gain_status.md
- stage: '142'
  role: paper_stage_tex
  path: paper/stages/stage_142.tex
- stage: '142'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md
- stage: '142'
  role: sympy_script
  path: scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py
- stage: '142'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl
- stage: '143'
  role: paper_stage_tex
  path: paper/stages/stage_143.tex
- stage: '143'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage143_equal_normalized_singular_limit.md
- stage: '143'
  role: sympy_script
  path: scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py
- stage: '143'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl
- stage: '144'
  role: paper_stage_tex
  path: paper/stages/stage_144.tex
- stage: '144'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md
- stage: '144'
  role: sympy_script
  path: scripts/moving_throat_pde_stage144_unique_regular_canonical_branch_sympy_audit.py
- stage: '144'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl
- stage: '145'
  role: paper_stage_tex
  path: paper/stages/stage_145.tex
- stage: '145'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage145_mouth_branch_selection_status.md
- stage: '146'
  role: paper_stage_tex
  path: paper/stages/stage_146.tex
- stage: '146'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage146_positive_deformation_expansion.md
- stage: '146'
  role: sympy_script
  path: scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py
- stage: '146'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl
- stage: '147'
  role: paper_stage_tex
  path: paper/stages/stage_147.tex
- stage: '147'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage147_first_order_rigidity_kernel.md
- stage: '147'
  role: sympy_script
  path: scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py
- stage: '147'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl
- stage: '148'
  role: paper_stage_tex
  path: paper/stages/stage_148.tex
- stage: '148'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage148_representative_positive_families.md
- stage: '148'
  role: sympy_script
  path: scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py
- stage: '148'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl
- stage: '149'
  role: paper_stage_tex
  path: paper/stages/stage_149.tex
- stage: '149'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage149_mouth_rigidity_status.md
- stage: '150'
  role: paper_stage_tex
  path: paper/stages/stage_150.tex
- stage: '150'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage150_full_profile_residual.md
- stage: '150'
  role: sympy_script
  path: scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py
- stage: '150'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl
- stage: '151'
  role: paper_stage_tex
  path: paper/stages/stage_151.tex
- stage: '151'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md
- stage: '151'
  role: sympy_script
  path: scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py
- stage: '151'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl
- stage: '152'
  role: paper_stage_tex
  path: paper/stages/stage_152.tex
- stage: '152'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage152_family1_actual_correction.md
- stage: '152'
  role: sympy_script
  path: scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py
- stage: '152'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl
- stage: '153'
  role: paper_stage_tex
  path: paper/stages/stage_153.tex
- stage: '153'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage153_full_mouth_correction_status.md
- stage: '154'
  role: paper_stage_tex
  path: paper/stages/stage_154.tex
- stage: '154'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage154_coevolving_core_mouth_map.md
- stage: '154'
  role: sympy_script
  path: scripts/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.py
- stage: '154'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl
- stage: '155'
  role: paper_stage_tex
  path: paper/stages/stage_155.tex
- stage: '155'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage155_frozen_traction_fixedpoint.md
- stage: '155'
  role: sympy_script
  path: scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py
- stage: '155'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl
- stage: '156'
  role: paper_stage_tex
  path: paper/stages/stage_156.tex
- stage: '156'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
- stage: '156'
  role: sympy_script
  path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
- stage: '156'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage156_renormalized_canonical_branch_mathematica_audit.wl
- stage: '157'
  role: paper_stage_tex
  path: paper/stages/stage_157.tex
- stage: '157'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md
- stage: '157'
  role: sympy_script
  path: scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py
- stage: '157'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl
- stage: '158'
  role: paper_stage_tex
  path: paper/stages/stage_158.tex
- stage: '158'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage158_linear_defect_transport.md
- stage: '158'
  role: sympy_script
  path: scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py
- stage: '158'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl
- stage: '159'
  role: paper_stage_tex
  path: paper/stages/stage_159.tex
- stage: '159'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage159_hybrid_outlet_projection.md
- stage: '159'
  role: sympy_script
  path: scripts/moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.py
- stage: '159'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage159_hybrid_outlet_projection_mathematica_audit.wl
- stage: '160'
  role: paper_stage_tex
  path: paper/stages/stage_160.tex
- stage: '160'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage160_bare_mixed_port_slippage.md
- stage: '160'
  role: sympy_script
  path: scripts/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.py
- stage: '160'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl
- stage: '161'
  role: paper_stage_tex
  path: paper/stages/stage_161.tex
- stage: '161'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage161_dn_similarity_slippage.md
- stage: '161'
  role: sympy_script
  path: scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py
- stage: '161'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl
- stage: '162'
  role: paper_stage_tex
  path: paper/stages/stage_162.tex
- stage: '162'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage162_parent_compensation_rigidity.md
- stage: '162'
  role: sympy_script
  path: scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py
- stage: '162'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage162_parent_compensation_rigidity_mathematica_audit.wl
- stage: '163'
  role: paper_stage_tex
  path: paper/stages/stage_163.tex
- stage: '163'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage163_off_family_normal_coordinate.md
- stage: '163'
  role: sympy_script
  path: scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py
- stage: '163'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl
- stage: '164'
  role: paper_stage_tex
  path: paper/stages/stage_164.tex
- stage: '164'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage164_microscopic_log_channels.md
- stage: '164'
  role: sympy_script
  path: scripts/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.py
- stage: '164'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl
- stage: '165'
  role: paper_stage_tex
  path: paper/stages/stage_165.tex
- stage: '165'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage165_exact_branch_drifts.md
- stage: '165'
  role: sympy_script
  path: scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py
- stage: '165'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl
- stage: '166'
  role: paper_stage_tex
  path: paper/stages/stage_166.tex
- stage: '166'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage166_bundle_inversion_four_drifts.md
- stage: '166'
  role: sympy_script
  path: scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py
- stage: '166'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl
- stage: '167'
  role: paper_stage_tex
  path: paper/stages/stage_167.tex
- stage: '167'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage167_bundle_transport_tangent_compensation.md
- stage: '167'
  role: sympy_script
  path: scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py
- stage: '167'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl
- stage: '168'
  role: paper_stage_tex
  path: paper/stages/stage_168.tex
- stage: '168'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage168_off_bundle_slippage.md
- stage: '168'
  role: sympy_script
  path: scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py
- stage: '168'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl
- stage: '169'
  role: paper_stage_tex
  path: paper/stages/stage_169.tex
- stage: '169'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage169_no_linear_p2_scalar_slippage.md
- stage: '169'
  role: sympy_script
  path: scripts/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py
- stage: '169'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl
- stage: '170'
  role: paper_stage_tex
  path: paper/stages/stage_170.tex
- stage: '170'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md
- stage: '170'
  role: sympy_script
  path: scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py
- stage: '170'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl
- stage: '171'
  role: paper_stage_tex
  path: paper/stages/stage_171.tex
- stage: '171'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage171_microscopic_grouped_obstructions.md
- stage: '171'
  role: sympy_script
  path: scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py
- stage: '171'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl
- stage: '172'
  role: paper_stage_tex
  path: paper/stages/stage_172.tex
- stage: '172'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage172_physical_slope_collapse.md
- stage: '172'
  role: sympy_script
  path: scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py
- stage: '172'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl
- stage: '173'
  role: paper_stage_tex
  path: paper/stages/stage_173.tex
- stage: '173'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage173_axisymmetric_loading_mismatch.md
- stage: '173'
  role: sympy_script
  path: scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py
- stage: '173'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl
- stage: '174'
  role: paper_stage_tex
  path: paper/stages/stage_174.tex
- stage: '174'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage174_static_self_similarity.md
- stage: '174'
  role: sympy_script
  path: scripts/moving_throat_pde_stage174_static_self_similarity_sympy_audit.py
- stage: '174'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl
- stage: '175'
  role: paper_stage_tex
  path: paper/stages/stage_175.tex
- stage: '175'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage175_wall_normalized_load_shape.md
- stage: '175'
  role: sympy_script
  path: scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py
- stage: '175'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl
- stage: '176'
  role: paper_stage_tex
  path: paper/stages/stage_176.tex
- stage: '176'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage176_outgoing_load_factorization.md
- stage: '176'
  role: sympy_script
  path: scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py
- stage: '176'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl
- stage: '177'
  role: paper_stage_tex
  path: paper/stages/stage_177.tex
- stage: '177'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md
- stage: '177'
  role: sympy_script
  path: scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py
- stage: '177'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl
- stage: '178'
  role: paper_stage_tex
  path: paper/stages/stage_178.tex
- stage: '178'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage178_outgoing_port_coloading_theorem.md
- stage: '178'
  role: sympy_script
  path: scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py
- stage: '178'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl
- stage: '179'
  role: paper_stage_tex
  path: paper/stages/stage_179.tex
- stage: '179'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage179_transfer_shape_theorem.md
- stage: '179'
  role: sympy_script
  path: scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py
- stage: '179'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl
- stage: '180'
  role: paper_stage_tex
  path: paper/stages/stage_180.tex
- stage: '180'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage180_effective_transfer_shape_collapse.md
- stage: '180'
  role: sympy_script
  path: scripts/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py
- stage: '180'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl
- stage: '181'
  role: paper_stage_tex
  path: paper/stages/stage_181.tex
- stage: '181'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage181_coherent_tracking_defect.md
- stage: '181'
  role: sympy_script
  path: scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py
- stage: '181'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl
- stage: '182'
  role: paper_stage_tex
  path: paper/stages/stage_182.tex
- stage: '182'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage182_microscopic_coherent_slippage.md
- stage: '182'
  role: sympy_script
  path: scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py
- stage: '182'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl
- stage: '183'
  role: paper_stage_tex
  path: paper/stages/stage_183.tex
- stage: '183'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage183_triangular_normal_form.md
- stage: '183'
  role: sympy_script
  path: scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py
- stage: '183'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl
- stage: '184'
  role: paper_stage_tex
  path: paper/stages/stage_184.tex
- stage: '184'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage184_branch_invariant_coordinates.md
- stage: '184'
  role: sympy_script
  path: scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py
- stage: '184'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl
- stage: '185'
  role: paper_stage_tex
  path: paper/stages/stage_185.tex
- stage: '185'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage185_microscopic_monomials.md
- stage: '185'
  role: sympy_script
  path: scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py
- stage: '185'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl
- stage: '186'
  role: paper_stage_tex
  path: paper/stages/stage_186.tex
- stage: '186'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage186_similarity_orbit_closure.md
- stage: '186'
  role: sympy_script
  path: scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py
- stage: '186'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl
- stage: '187'
  role: paper_stage_tex
  path: paper/stages/stage_187.tex
- stage: '187'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage187_orbit_quotient_closure.md
- stage: '187'
  role: sympy_script
  path: scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py
- stage: '187'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl
- stage: '188'
  role: paper_stage_tex
  path: paper/stages/stage_188.tex
- stage: '188'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage188_branch_observables_completion.md
- stage: '188'
  role: sympy_script
  path: scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py
- stage: '188'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.wl
- stage: '189'
  role: paper_stage_tex
  path: paper/stages/stage_189.tex
- stage: '189'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md
- stage: '189'
  role: sympy_script
  path: scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py
- stage: '189'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.wl
- stage: '190'
  role: paper_stage_tex
  path: paper/stages/stage_190.tex
- stage: '190'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage190_direct_defect_vs_dressing_split.md
- stage: '190'
  role: sympy_script
  path: scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py
- stage: '190'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl
- stage: '191'
  role: paper_stage_tex
  path: paper/stages/stage_191.tex
- stage: '191'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md
- stage: '191'
  role: sympy_script
  path: scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py
- stage: '191'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.wl
- stage: '192'
  role: paper_stage_tex
  path: paper/stages/stage_192.tex
- stage: '192'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
- stage: '192'
  role: sympy_script
  path: scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py
- stage: '192'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage192_orbit_quotient_projectors_mathematica_audit.wl
- stage: '193'
  role: paper_stage_tex
  path: paper/stages/stage_193.tex
- stage: '193'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md
- stage: '193'
  role: sympy_script
  path: scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py
- stage: '193'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_mathematica_audit.wl
- stage: '194'
  role: paper_stage_tex
  path: paper/stages/stage_194.tex
- stage: '194'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra.md
- stage: '194'
  role: sympy_script
  path: scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py
- stage: '194'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.wl
- stage: '195'
  role: paper_stage_tex
  path: paper/stages/stage_195.tex
- stage: '195'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md
- stage: '195'
  role: sympy_script
  path: scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py
- stage: '195'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.wl
- stage: '196'
  role: paper_stage_tex
  path: paper/stages/stage_196.tex
- stage: '196'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md
- stage: '196'
  role: sympy_script
  path: scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py
- stage: '196'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.wl
- stage: '197'
  role: paper_stage_tex
  path: paper/stages/stage_197.tex
- stage: '197'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage197_conditional_packetA_closure_theorem.md
- stage: '197'
  role: sympy_script
  path: scripts/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.py
- stage: '197'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage197_conditional_packetA_closure_theorem_mathematica_audit.wl
- stage: '198'
  role: paper_stage_tex
  path: paper/stages/stage_198.tex
- stage: '198'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage198_exact_finite_orbit_law.md
- stage: '198'
  role: sympy_script
  path: scripts/moving_throat_pde_stage198_exact_finite_orbit_law_sympy_audit.py
- stage: '198'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage198_exact_finite_orbit_law_mathematica_audit.wl
- stage: '199'
  role: paper_stage_tex
  path: paper/stages/stage_199.tex
- stage: '199'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage199_pairwise_orbit_transport_law.md
- stage: '199'
  role: sympy_script
  path: scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py
- stage: '199'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage199_pairwise_orbit_transport_law_mathematica_audit.wl
- stage: '200'
  role: paper_stage_tex
  path: paper/stages/stage_200.tex
- stage: '200'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage200_reference_free_home_stretch_theorem.md
- stage: '200'
  role: sympy_script
  path: scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py
- stage: '200'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl
- stage: '201'
  role: paper_stage_tex
  path: paper/stages/stage_201.tex
- stage: '201'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection.md
- stage: '201'
  role: sympy_script
  path: scripts/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.py
- stage: '201'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.wl
- stage: '202'
  role: paper_stage_tex
  path: paper/stages/stage_202.tex
- stage: '202'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md
- stage: '202'
  role: sympy_script
  path: scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py
- stage: '202'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl
- stage: '203'
  role: paper_stage_tex
  path: paper/stages/stage_203.tex
- stage: '203'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem.md
- stage: '203'
  role: sympy_script
  path: scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py
- stage: '203'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl
- stage: '204'
  role: paper_stage_tex
  path: paper/stages/stage_204.tex
- stage: '204'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor.md
- stage: '204'
  role: sympy_script
  path: scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py
- stage: '204'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.wl
- stage: '205'
  role: paper_stage_tex
  path: paper/stages/stage_205.tex
- stage: '205'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md
- stage: '205'
  role: sympy_script
  path: scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py
- stage: '205'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl
- stage: '206'
  role: paper_stage_tex
  path: paper/stages/stage_206.tex
- stage: '206'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md
- stage: '206'
  role: sympy_script
  path: scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py
- stage: '206'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl
- stage: '207'
  role: paper_stage_tex
  path: paper/stages/stage_207.tex
- stage: '207'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table.md
- stage: '207'
  role: sympy_script
  path: scripts/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.py
- stage: '207'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl
- stage: '208'
  role: paper_stage_tex
  path: paper/stages/stage_208.tex
- stage: '208'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md
- stage: '208'
  role: sympy_script
  path: scripts/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.py
- stage: '208'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl
- stage: '209'
  role: paper_stage_tex
  path: paper/stages/stage_209.tex
- stage: '209'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md
- stage: '209'
  role: sympy_script
  path: scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py
- stage: '209'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl
- stage: '210'
  role: paper_stage_tex
  path: paper/stages/stage_210.tex
- stage: '210'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md
- stage: '210'
  role: sympy_script
  path: scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py
- stage: '210'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl
- stage: '211'
  role: paper_stage_tex
  path: paper/stages/stage_211.tex
- stage: '211'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction.md
- stage: '211'
  role: sympy_script
  path: scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py
- stage: '211'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl
- stage: '212'
  role: paper_stage_tex
  path: paper/stages/stage_212.tex
- stage: '212'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md
- stage: '212'
  role: sympy_script
  path: scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py
- stage: '212'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl
- stage: '213'
  role: paper_stage_tex
  path: paper/stages/stage_213.tex
- stage: '213'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md
- stage: '213'
  role: sympy_script
  path: scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py
- stage: '213'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl
- stage: '214'
  role: paper_stage_tex
  path: paper/stages/stage_214.tex
- stage: '214'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md
- stage: '214'
  role: sympy_script
  path: scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py
- stage: '214'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl
- stage: '215'
  role: paper_stage_tex
  path: paper/stages/stage_215.tex
- stage: '215'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem.md
- stage: '215'
  role: sympy_script
  path: scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py
- stage: '215'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl
- stage: '216'
  role: paper_stage_tex
  path: paper/stages/stage_216.tex
- stage: '216'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md
- stage: '216'
  role: sympy_script
  path: scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py
- stage: '216'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl
- stage: '217'
  role: paper_stage_tex
  path: paper/stages/stage_217.tex
- stage: '217'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md
- stage: '217'
  role: sympy_script
  path: scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py
- stage: '217'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl
- stage: '218'
  role: paper_stage_tex
  path: paper/stages/stage_218.tex
- stage: '218'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md
- stage: '218'
  role: sympy_script
  path: scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py
- stage: '218'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl
- stage: '219'
  role: paper_stage_tex
  path: paper/stages/stage_219.tex
- stage: '219'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md
- stage: '219'
  role: sympy_script
  path: scripts/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.py
- stage: '219'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl
- stage: '220'
  role: paper_stage_tex
  path: paper/stages/stage_220.tex
- stage: '220'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md
- stage: '220'
  role: sympy_script
  path: scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py
- stage: '220'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl
- stage: '221'
  role: paper_stage_tex
  path: paper/stages/stage_221.tex
- stage: '221'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md
- stage: '221'
  role: sympy_script
  path: scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py
- stage: '221'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl
- stage: '222'
  role: paper_stage_tex
  path: paper/stages/stage_222.tex
- stage: '222'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
- stage: '222'
  role: sympy_script
  path: scripts/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py
- stage: '222'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_mathematica_audit.wl
- stage: '223'
  role: paper_stage_tex
  path: paper/stages/stage_223.tex
- stage: '223'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
- stage: '223'
  role: sympy_script
  path: scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py
- stage: '223'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_mathematica_audit.wl
- stage: '224'
  role: paper_stage_tex
  path: paper/stages/stage_224.tex
- stage: '224'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
- stage: '224'
  role: sympy_script
  path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
- stage: '224'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl
- stage: '225'
  role: paper_stage_tex
  path: paper/stages/stage_225.tex
- stage: '225'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md
- stage: '225'
  role: sympy_script
  path: scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py
- stage: '225'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.wl
- stage: '226'
  role: paper_stage_tex
  path: paper/stages/stage_226.tex
- stage: '226'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md
- stage: '226'
  role: sympy_script
  path: scripts/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.py
- stage: '226'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl
- stage: '227'
  role: paper_stage_tex
  path: paper/stages/stage_227.tex
- stage: '227'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md
- stage: '227'
  role: sympy_script
  path: scripts/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py
- stage: '227'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.wl
- stage: '228'
  role: paper_stage_tex
  path: paper/stages/stage_228.tex
- stage: '228'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md
- stage: '228'
  role: sympy_script
  path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
- stage: '228'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.wl
- stage: '229'
  role: paper_stage_tex
  path: paper/stages/stage_229.tex
- stage: '229'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.md
- stage: '229'
  role: sympy_script
  path: scripts/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py
- stage: '229'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage229_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_mathematica_audit.wl
- stage: '230'
  role: paper_stage_tex
  path: paper/stages/stage_230.tex
- stage: '230'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md
- stage: '230'
  role: sympy_script
  path: scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py
- stage: '230'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl
- stage: '231'
  role: paper_stage_tex
  path: paper/stages/stage_231.tex
- stage: '231'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md
- stage: '231'
  role: sympy_script
  path: scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py
- stage: '231'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl
- stage: '232'
  role: paper_stage_tex
  path: paper/stages/stage_232.tex
- stage: '232'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md
- stage: '232'
  role: sympy_script
  path: scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py
- stage: '232'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl
- stage: '233'
  role: paper_stage_tex
  path: paper/stages/stage_233.tex
- stage: '233'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.md
- stage: '233'
  role: sympy_script
  path: scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py
- stage: '233'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_mathematica_audit.wl
- stage: '234'
  role: paper_stage_tex
  path: paper/stages/stage_234.tex
- stage: '234'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md
- stage: '234'
  role: sympy_script
  path: scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py
- stage: '234'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.wl
- stage: '235'
  role: paper_stage_tex
  path: paper/stages/stage_235.tex
- stage: '235'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md
- stage: '235'
  role: sympy_script
  path: scripts/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.py
- stage: '235'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl
- stage: '236'
  role: paper_stage_tex
  path: paper/stages/stage_236.tex
- stage: '236'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md
- stage: '236'
  role: sympy_script
  path: scripts/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.py
- stage: '236'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl
- stage: '237'
  role: paper_stage_tex
  path: paper/stages/stage_237.tex
- stage: '237'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.md
- stage: '237'
  role: sympy_script
  path: scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py
- stage: '237'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl
- stage: '238'
  role: paper_stage_tex
  path: paper/stages/stage_238.tex
- stage: '238'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.md
- stage: '238'
  role: sympy_script
  path: scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py
- stage: '238'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.wl
- stage: '239'
  role: paper_stage_tex
  path: paper/stages/stage_239.tex
- stage: '239'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md
- stage: '239'
  role: sympy_script
  path: scripts/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py
- stage: '239'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl
- stage: '240'
  role: paper_stage_tex
  path: paper/stages/stage_240.tex
- stage: '240'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.md
- stage: '240'
  role: sympy_script
  path: scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py
- stage: '240'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl
- stage: '241'
  role: paper_stage_tex
  path: paper/stages/stage_241.tex
- stage: '241'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md
- stage: '241'
  role: sympy_script
  path: scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py
- stage: '241'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl
- stage: '242'
  role: paper_stage_tex
  path: paper/stages/stage_242.tex
- stage: '242'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.md
- stage: '242'
  role: sympy_script
  path: scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py
- stage: '242'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl
- stage: '243'
  role: paper_stage_tex
  path: paper/stages/stage_243.tex
- stage: '243'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md
- stage: '243'
  role: sympy_script
  path: scripts/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py
- stage: '243'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl
- stage: '244'
  role: paper_stage_tex
  path: paper/stages/stage_244.tex
- stage: '244'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md
- stage: '244'
  role: sympy_script
  path: scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py
- stage: '244'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl
- stage: '245'
  role: paper_stage_tex
  path: paper/stages/stage_245.tex
- stage: '245'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md
- stage: '245'
  role: sympy_script
  path: scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py
- stage: '245'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl
- stage: '246'
  role: paper_stage_tex
  path: paper/stages/stage_246.tex
- stage: '246'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.md
- stage: '246'
  role: sympy_script
  path: scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py
- stage: '246'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl
- stage: '247'
  role: paper_stage_tex
  path: paper/stages/stage_247.tex
- stage: '247'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
- stage: '247'
  role: sympy_script
  path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
- stage: '247'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.wl
- stage: '248'
  role: paper_stage_tex
  path: paper/stages/stage_248.tex
- stage: '248'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
- stage: '248'
  role: sympy_script
  path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
- stage: '248'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl
- stage: '249'
  role: paper_stage_tex
  path: paper/stages/stage_249.tex
- stage: '249'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.md
- stage: '249'
  role: sympy_script
  path: scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py
- stage: '249'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl
- stage: '250'
  role: paper_stage_tex
  path: paper/stages/stage_250.tex
- stage: '250'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
- stage: '250'
  role: sympy_script
  path: scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
- stage: '250'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl
- stage: '251'
  role: paper_stage_tex
  path: paper/stages/stage_251.tex
- stage: '251'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md
- stage: '251'
  role: sympy_script
  path: scripts/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.py
- stage: '251'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl
- stage: '252'
  role: paper_stage_tex
  path: paper/stages/stage_252.tex
- stage: '252'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.md
- stage: '252'
  role: sympy_script
  path: scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py
- stage: '252'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl
- stage: '253'
  role: paper_stage_tex
  path: paper/stages/stage_253.tex
- stage: '253'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md
- stage: '253'
  role: sympy_script
  path: scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py
- stage: '253'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`

Task:
Find every claim label or prose claim that says or implies a value is derived, forced, exact, non-tunable, matched, canonical, not free, or fixed. The target is the claim, not whether the claim is true.

Emit only YAML:

```yaml
modality: claim_label
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
```
