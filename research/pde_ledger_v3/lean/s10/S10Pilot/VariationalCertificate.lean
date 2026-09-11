import S10Pilot.FiniteAction
import S10Pilot.Certificate

/-! S10's mode census in arbitrary dimension, now characterized by the integrated variational principle. -/

namespace S10Pilot
noncomputable section

variable {D : ℕ}

/-- The S10 census expressed through stationarity against all smooth compact test fields. -/
theorem s10_variational_certificate {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue rho mu k ∧
      (∀ a : Vec D, ActionStationary rho mu (planeWave omega k a) ↔ a ∈ transverseSpace k) ∧
      Module.finrank ℝ (transverseSpace k) = D - 1 ∧
      (∀ a : Vec D, ActionStationary rho mu (planeWave 0 k a) ↔ a ∈ longitudinalSpace k) ∧
      Module.finrank ℝ (longitudinalSpace k) = 1 := by
  obtain ⟨omega, hw, hf, ht, htdim, hl, hldim⟩ := s10_plane_wave_certificate hrho hmu hk
  refine ⟨omega, hw, hf, ?_, htdim, ?_, hldim⟩
  · intro a
    rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a)]
    exact ht a
  · intro a
    rw [actionStationary_iff_eulerLagrange (smooth_planeWave 0 k a)]
    exact hl a

theorem propagating_variational_mode_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (homega : omega ≠ 0) (ha : a ≠ 0) :
    ActionStationary rho mu (planeWave omega k a) ↔
      dot k a = 0 ∧ omega ^ 2 = coneValue rho mu k := by
  rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a)]
  exact propagating_planeWave_iff hrho homega ha

end
end S10Pilot
