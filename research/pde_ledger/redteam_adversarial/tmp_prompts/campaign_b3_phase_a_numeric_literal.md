# Phase A Numeric-Literal Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 082, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 093, 094, 095, 096`
- Source files: `- stage: '065'
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
  path: mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl`

Task:
Find numerical literals, rational coefficients, closed-form constants, or parameter assignments that could be fit-insertion-point candidates. Exclude pure identities such as harmonic normalizations, dimension labels, residual tolerances, line numbers, and pass/fail counters unless the source itself claims the value is matched, fixed, canonical, forced, or derived.

Emit only YAML:

```yaml
modality: numeric_literal
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
