import Physlib.ClassicalMechanics.WaveEquation.Basic
import S9Pilot.Certificate
import S9Pilot.VariationalCertificate
import S9Pilot.Madelung

/-!
S9 formalization pilot: integrated variation, local PDE, and plane-wave spectrum.

The original elementary installation check is retained below for provenance.
The derivation and integrated variational principle are in the imported modules
under S9Pilot/. See RESULT.md for assumptions and the remaining proof boundary.
-/

namespace S9Pilot

/-- Original scalar dispersion expression retained as an installation check. -/
noncomputable def candidateOmegaSq (rho mu k : ℝ) : ℝ := (mu / rho) * k ^ 2

/-- The candidate expression is homogeneous of degree two in wave number. -/
theorem candidateOmegaSq_scaling (rho mu k scale : ℝ) :
    candidateOmegaSq rho mu (scale * k) = scale ^ 2 * candidateOmegaSq rho mu k := by
  unfold candidateOmegaSq
  ring

end S9Pilot

-- Verify that the physics library is available and display proof dependencies.
#check ClassicalMechanics.planeWave_waveEquation
#print axioms S9Pilot.candidateOmegaSq_scaling
#print axioms ClassicalMechanics.planeWave_waveEquation
#print axioms S9Pilot.lagrangian_variation
#print axioms S9Pilot.modal_stationary_iff
#print axioms S9Pilot.planeWave_lagrangian
#print axioms S9Pilot.eulerLagrange_planeWave
#print axioms S9Pilot.s9_plane_wave_certificate
#print axioms S9Pilot.propagating_planeWave_iff
#print axioms S9Pilot.coneValue_scaling
#print axioms S9Pilot.concrete_transverse_wave
#print axioms S9Pilot.concrete_longitudinal_not_wave
#print axioms S9Pilot.coord_integration_by_parts
#print axioms S9Pilot.relative_density_integrable
#print axioms S9Pilot.relativeAction_hasDerivAt
#print axioms S9Pilot.relativeAction_deriv_eq_eulerLagrange
#print axioms S9Pilot.actionStationary_iff_eulerLagrange
#print axioms S9Pilot.finiteAction_hasDerivAt
#print axioms S9Pilot.s9_variational_certificate
#print axioms S9Pilot.propagating_variational_mode_iff
#print axioms S9Pilot.concrete_transverse_stationary
#print axioms S9Pilot.concrete_longitudinal_detected_by_test
#print axioms S9Pilot.Madelung.phasePerturbation_coordDeriv
#print axioms S9Pilot.Madelung.velocity_phasePerturbation
#print axioms S9Pilot.Madelung.linearVelocity_eq
#print axioms S9Pilot.Madelung.cosAmplitude_mem
#print axioms S9Pilot.Madelung.sinAmplitude_mem
#print axioms S9Pilot.Madelung.cosAmplitude_range
#print axioms S9Pilot.Madelung.longitudinal_transverse_eq_zero
#print axioms S9Pilot.Madelung.no_transverse_velocity
#print axioms S9Pilot.Madelung.transverse_linearVelocity_zero
#print axioms S9Pilot.Madelung.zero_wavevector
#print axioms S9Pilot.Madelung.concrete_longitudinal_velocity
