import S10Pilot.Specialization
import S10Pilot.PhaseAverage
import S10Pilot.EdgeCases

/-! Arbitrary-dimensional baseline S10 certificate, with explicit D=3 agreement
with S9. See RESULT.md for mathematical scope and remaining controls. -/

#print axioms S10Pilot.lagrangian_variation
#print axioms S10Pilot.antisym_mode_pair
#print axioms S10Pilot.modalAction_eq
#print axioms S10Pilot.modalAction_variation
#print axioms S10Pilot.modal_stationary_iff
#print axioms S10Pilot.momentum_eq
#print axioms S10Pilot.eulerLagrange_planeWave
#print axioms S10Pilot.phaseAverage_eq
#print axioms S10Pilot.phaseAverage_variation
#print axioms S10Pilot.phaseAverage_stationary_iff
#print axioms S10Pilot.coord_integration_by_parts
#print axioms S10Pilot.relative_density_integrable
#print axioms S10Pilot.relativeAction_hasDerivAt
#print axioms S10Pilot.actionStationary_iff_eulerLagrange
#print axioms S10Pilot.finiteAction_hasDerivAt
#print axioms S10Pilot.s10_plane_wave_certificate
#print axioms S10Pilot.s10_variational_certificate
#print axioms S10Pilot.propagating_variational_mode_iff
#print axioms S10Pilot.coneValue_scaling
#print axioms S10Pilot.stiffness_three
#print axioms S10Pilot.lagrangian_three
#print axioms S10Pilot.modalOperator_three
#print axioms S10Pilot.coneValue_three
#print axioms S10Pilot.transverseSpace_three
#print axioms S10Pilot.longitudinalSpace_three
#print axioms S10Pilot.planeWave_three
#print axioms S10Pilot.eulerLagrange_three
#print axioms S10Pilot.relativeAction_three
#print axioms S10Pilot.actionStationary_three
#print axioms S10Pilot.exists_propagating_mode
#print axioms S10Pilot.dimension_one_stiffness
#print axioms S10Pilot.dimension_one_no_propagating_mode
#print axioms S10Pilot.concrete_five_transverse_stationary
#print axioms S10Pilot.concrete_five_longitudinal_detected
