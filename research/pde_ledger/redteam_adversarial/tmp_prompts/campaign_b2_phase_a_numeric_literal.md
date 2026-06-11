# Phase A Numeric-Literal Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `033, 034, 035, 036, 037, 038, 039, 040, 041, 042, 043, 044, 045, 046, 047, 048, 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064`
- Source files: `- stage: '033'
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
  path: mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`

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
