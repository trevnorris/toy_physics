import S9Pilot.FiniteAction
import S9Pilot.Certificate

/-! S9's mode census, now characterized by the integrated variational principle. -/

namespace S9Pilot
noncomputable section

/-- The S9 census expressed through stationarity against all smooth compact test fields. -/
theorem s9_variational_certificate {rho mu : ℝ} {k : Vec}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue rho mu k ∧
      (∀ a : Vec, ActionStationary rho mu (planeWave omega k a) ↔ a ∈ transverseSpace k) ∧
      Module.finrank ℝ (transverseSpace k) = 2 ∧
      (∀ a : Vec, ActionStationary rho mu (planeWave 0 k a) ↔ a ∈ longitudinalSpace k) ∧
      Module.finrank ℝ (longitudinalSpace k) = 1 := by
  obtain ⟨omega, hw, hf, ht, htdim, hl, hldim⟩ := s9_plane_wave_certificate hrho hmu hk
  refine ⟨omega, hw, hf, ?_, htdim, ?_, hldim⟩
  · intro a
    rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a)]
    exact ht a
  · intro a
    rw [actionStationary_iff_eulerLagrange (smooth_planeWave 0 k a)]
    exact hl a

theorem propagating_variational_mode_iff {rho mu omega : ℝ} {k a : Vec}
    (hrho : rho ≠ 0) (homega : omega ≠ 0) (ha : a ≠ 0) :
    ActionStationary rho mu (planeWave omega k a) ↔
      dot k a = 0 ∧ omega ^ 2 = coneValue rho mu k := by
  rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a)]
  exact propagating_planeWave_iff hrho homega ha

/-- The transverse witness is stationary against every compact test variation. -/
theorem concrete_transverse_stationary :
    ActionStationary 1 1 (planeWave 1 ![0, 0, 1] ![1, 0, 0]) := by
  rw [actionStationary_iff_eulerLagrange (smooth_planeWave _ _ _)]
  exact concrete_transverse_wave

/-- There exists a smooth compact test variation detecting the longitudinal nonsolution. -/
theorem concrete_longitudinal_detected_by_test :
    ∃ h : Point → Vec, TestField h ∧
      deriv (relativeAction 1 1 (planeWave 1 ![0, 0, 1] ![0, 0, 1]) h) 0 ≠ 0 := by
  have hn : ¬ ActionStationary 1 1 (planeWave 1 ![0, 0, 1] ![0, 0, 1]) := by
    rw [actionStationary_iff_eulerLagrange (smooth_planeWave _ _ _)]
    exact concrete_longitudinal_not_wave
  simpa [ActionStationary] using hn

end
end S9Pilot
