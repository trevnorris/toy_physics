# Phase A Numeric-Literal Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `001, 002, 003, 004, 005, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 023, 024, 025, 026, 027, 028, 029, 030, 031, 032`
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
  path: mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`

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
