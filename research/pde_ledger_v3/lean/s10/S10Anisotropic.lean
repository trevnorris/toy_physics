import S10Anisotropic.Checks
import S10Anisotropic.Scaling

/-! One-axis anisotropic inertia: supplied action, integrated variation,
complete real plane-wave spectrum, and directional strata.
See ANISOTROPIC_RESULT.md for scope and the additional perpendicular stratum. -/

#print axioms S10Anisotropic.kinetic_eq
#print axioms S10Anisotropic.lagrangian_variation
#print axioms S10Anisotropic.modalAction_eq
#print axioms S10Anisotropic.modalAction_variation
#print axioms S10Anisotropic.momentum_eq
#print axioms S10Anisotropic.momentum_contraction
#print axioms S10Anisotropic.eulerLagrange_planeWave
#print axioms S10Anisotropic.relative_density_integrable
#print axioms S10Anisotropic.relativeAction_hasDerivAt
#print axioms S10Anisotropic.actionStationary_iff_eulerLagrange
#print axioms S10Anisotropic.finiteAction_hasDerivAt
#print axioms S10Anisotropic.phaseAverage_eq
#print axioms S10Anisotropic.phaseAverage_variation
#print axioms S10Anisotropic.phaseAverage_stationary_iff
#print axioms S10Anisotropic.perpSq_zero_iff
#print axioms S10Anisotropic.ordinarySpace_finrank
#print axioms S10Anisotropic.extraValue_pos
#print axioms S10Anisotropic.extraValue_ne_ordinary
#print axioms S10Anisotropic.nonzero_root_positive
#print axioms S10Anisotropic.modalOperator_normalized
#print axioms S10Anisotropic.weighted_transverse
#print axioms S10Anisotropic.split_propagating_iff
#print axioms S10Anisotropic.ordinary_kernel_iff
#print axioms S10Anisotropic.extra_kernel_iff
#print axioms S10Anisotropic.parallel_kernel_iff
#print axioms S10Anisotropic.zero_frequency_iff
#print axioms S10Anisotropic.oblique_counts
#print axioms S10Anisotropic.perpendicular_counts
#print axioms S10Anisotropic.parallel_counts
#print axioms S10Anisotropic.zero_counts
#print axioms S10Anisotropic.split_total_count
#print axioms S10Anisotropic.actionStationary_iff_mem_modeSpace
#print axioms S10Anisotropic.split_variational_iff
#print axioms S10Anisotropic.parallel_variational_iff
#print axioms S10Anisotropic.zero_variational_iff
#print axioms S10Anisotropic.split_variational_certificate
#print axioms S10Anisotropic.frequency_coincidence_iff
#print axioms S10Anisotropic.extra_exactly_transverse_iff
#print axioms S10Anisotropic.isotropic_variational_iff
#print axioms S10Anisotropic.extra_wave_exists
#print axioms S10Anisotropic.ordinary_wave_exists
#print axioms S10Anisotropic.concrete_oblique_wave
#print axioms S10Anisotropic.concrete_perpendicular_wave
#print axioms S10Anisotropic.concrete_perpendicular_transverse_dimension
#print axioms S10Anisotropic.concrete_one_axis_kinetic
#print axioms S10Anisotropic.dimension_one_no_propagating_mode
#print axioms S10Anisotropic.extraConeValue_scaling
#print axioms S10Anisotropic.extraConeValue_scaling_ratio
