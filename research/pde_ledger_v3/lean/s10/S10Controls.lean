import S10Controls.Checks
import S10Controls.PhaseAverage
import S10Controls.Scalar

/-! Full-gradient and divergence-only controls, from supplied actions through
integrated variation to complete real plane-wave amplitude classifications.
See CONTROLS_RESULT.md for scope and interpretation. -/

#print axioms S10Controls.lagrangian_variation
#print axioms S10Controls.fullGradient_mode_pair
#print axioms S10Controls.divergence_modeJet
#print axioms S10Controls.modalAction_eq
#print axioms S10Controls.modalAction_variation
#print axioms S10Controls.momentum_eq
#print axioms S10Controls.momentum_contraction
#print axioms S10Controls.eulerLagrange_planeWave
#print axioms S10Controls.relative_density_integrable
#print axioms S10Controls.relativeAction_hasDerivAt
#print axioms S10Controls.actionStationary_iff_eulerLagrange
#print axioms S10Controls.finiteAction_hasDerivAt
#print axioms S10Controls.phaseAverage_eq
#print axioms S10Controls.phaseAverage_variation
#print axioms S10Controls.phaseAverage_stationary_iff
#print axioms S10Controls.full_stationary_iff
#print axioms S10Controls.div_propagating_iff
#print axioms S10Controls.div_zero_frequency_iff
#print axioms S10Controls.actionStationary_iff_mem_amplitudeSpace
#print axioms S10Controls.full_variational_iff
#print axioms S10Controls.div_propagating_variational_iff
#print axioms S10Controls.full_on_cone_space
#print axioms S10Controls.full_zero_space
#print axioms S10Controls.div_on_cone_space
#print axioms S10Controls.div_zero_space
#print axioms S10Controls.full_mode_counts
#print axioms S10Controls.div_mode_counts
#print axioms S10Controls.stiffness_comparison_certificate
#print axioms S10Controls.longitudinal_control_witness
#print axioms S10Controls.transverse_control_witness
#print axioms S10Controls.dimension_one_stiffness_eq
#print axioms S10Controls.concrete_distinct_densities

#print axioms S10ScalarControls.lagrangian_variation
#print axioms S10ScalarControls.relativeAction_hasDerivAt
#print axioms S10ScalarControls.modalAction_eq
#print axioms S10ScalarControls.phaseAverage_eq
#print axioms S10ScalarControls.spectralAction_variation
#print axioms S10ScalarControls.variational_planeWave_iff
#print axioms S10ScalarControls.nonzero_root_iff
#print axioms S10ScalarControls.cone_modeSpace
#print axioms S10ScalarControls.zero_modeSpace
#print axioms S10ScalarControls.positive_variational_certificate
#print axioms S10ScalarControls.negative_control_no_real_wave
#print axioms S10ScalarControls.coneValue_scaling_ratio
#print axioms S10ScalarControls.signflip_counts
#print axioms S10ScalarControls.zero_counts
#print axioms S10ScalarControls.negative_root_exists
#print axioms S10ScalarControls.coefficient_changes_frequency
