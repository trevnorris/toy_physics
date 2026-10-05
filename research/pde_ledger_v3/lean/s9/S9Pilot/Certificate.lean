import S9Pilot.PlaneWave
import S9Pilot.Spectrum

/-! Combined theorem and concrete nonvacuity checks for the first S9 pilot. -/

namespace S9Pilot
noncomputable section

/-- For every nonzero spatial wavevector and positive supplied material constants,
there is a positive cone frequency. Its full real plane-wave amplitude space has
dimension two and is transverse. The zero-frequency space is longitudinal and
has dimension one. Both spaces solve the differential expression obtained from
the supplied density. This theorem is for D=3. -/
theorem s9_plane_wave_certificate {rho mu : ℝ} {k : Vec}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue rho mu k ∧
      (∀ a : Vec, (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
        a ∈ transverseSpace k) ∧
      Module.finrank ℝ (transverseSpace k) = 2 ∧
      (∀ a : Vec, (∀ x, eulerLagrange rho mu (planeWave 0 k a) x = 0) ↔
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
theorem propagating_planeWave_iff {rho mu omega : ℝ} {k a : Vec}
    (hrho : rho ≠ 0) (homega : omega ≠ 0) (ha : a ≠ 0) :
    (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
      dot k a = 0 ∧ omega ^ 2 = coneValue rho mu k := by
  rw [planeWave_solves_iff]
  exact propagating_mode_iff hrho homega ha

/-- A concrete transverse solution with all supplied physical inequalities satisfiable. -/
theorem concrete_transverse_wave :
    ∀ x, eulerLagrange 1 1 (planeWave 1 ![0, 0, 1] ![1, 0, 0]) x = 0 := by
  intro x
  rw [eulerLagrange_planeWave]
  have h : modalOperator 1 1 1 ![0, 0, 1] ![1, 0, 0] = 0 := by
    ext i
    fin_cases i <;> simp [modalOperator, normSq, dot, Matrix.cons_val_two]
  rw [h, smul_zero]

/-- At the same nonzero frequency, a purely longitudinal amplitude fails the PDE. -/
theorem concrete_longitudinal_not_wave :
    ¬ (∀ x, eulerLagrange 1 1 (planeWave 1 ![0, 0, 1] ![0, 0, 1]) x = 0) := by
  intro h
  have hz := congrFun (h 0) 2
  simp [eulerLagrange_planeWave, phase, modalOperator, normSq, dot, Matrix.cons_val_two] at hz

end
end S9Pilot
