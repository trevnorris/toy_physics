import S10Pilot.PlaneWave
import S10Pilot.Spectrum

/-! A dimension-general certificate for the baseline S10 action. -/

namespace S10Pilot
noncomputable section

variable {D : ℕ}

/-- For every nonzero wavevector and positive material coefficients, the cone
frequency is positive; its entire stationary amplitude space is transverse with
dimension D-1. The zero-frequency space is the one-dimensional longitudinal span.
At D=1 the positive cone parameter has no nonzero stationary amplitude. -/
theorem s10_plane_wave_certificate {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue rho mu k ∧
      (∀ a : Vec D, (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
        a ∈ transverseSpace k) ∧
      Module.finrank ℝ (transverseSpace k) = D - 1 ∧
      (∀ a : Vec D, (∀ x, eulerLagrange rho mu (planeWave 0 k a) x = 0) ↔
        a ∈ longitudinalSpace k) ∧
      Module.finrank ℝ (longitudinalSpace k) = 1 := by
  have hp := coneValue_pos hrho hmu hk
  have hw : 0 < Real.sqrt (coneValue rho mu k) := Real.sqrt_pos.mpr hp
  have hf := Real.sq_sqrt (le_of_lt hp)
  refine ⟨Real.sqrt (coneValue rho mu k), hw, hf, ?_, transverseSpace_finrank hk,
    ?_, longitudinalSpace_finrank hk⟩
  · intro a
    rw [planeWave_solves_iff]
    exact stationary_on_cone_iff (ne_of_gt hrho) (ne_of_gt hw) hf
  · intro a
    rw [planeWave_solves_iff]
    exact zero_frequency_iff (ne_of_gt hmu) hk

/-- Completeness for arbitrary nonzero real frequency, rather than only the root
chosen in the certificate. -/
theorem propagating_planeWave_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (homega : omega ≠ 0) (ha : a ≠ 0) :
    (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
      dot k a = 0 ∧ omega ^ 2 = coneValue rho mu k := by
  rw [planeWave_solves_iff]
  exact propagating_mode_iff hrho homega ha

end
end S10Pilot
