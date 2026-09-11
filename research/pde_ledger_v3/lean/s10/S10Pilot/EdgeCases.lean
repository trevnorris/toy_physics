import S10Pilot.VariationalCertificate

/-! Nonvacuity and edge cases, including dimension one and axis-aligned vectors. -/

namespace S10Pilot
noncomputable section

/-- The cone contains nonzero stationary amplitudes in every D ≥ 2. -/
theorem exists_propagating_mode {D : ℕ} {rho mu : ℝ} {k : Vec D}
    (hD : 2 ≤ D) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, ∃ a : Vec D, 0 < omega ∧ a ≠ 0 ∧
      ActionStationary rho mu (planeWave omega k a) := by
  obtain ⟨w, hw, _, ht, hdim, _, _⟩ := s10_variational_certificate hrho hmu hk
  have hp : 0 < Module.finrank ℝ (transverseSpace k) := by
    rw [hdim]
    omega
  obtain ⟨a, ha⟩ := Module.finrank_pos_iff_exists_ne_zero.mp hp
  refine ⟨w, a, hw, ?_, (ht a).mpr a.property⟩
  intro hz
  apply ha
  exact Subtype.ext hz

theorem dimension_one_stiffness (J : Jet 1) : stiffness J = 0 := by
  simp [stiffness, antisym]

theorem dimension_one_operator (rho mu omega : ℝ) (k a : Vec 1) :
    modalOperator rho mu omega k a = (rho * omega ^ 2) • a := by
  ext i
  have hi : i = 0 := Subsingleton.elim _ _
  subst i
  simp only [modalOperator, normSq, dot, Fin.sum_univ_one, Pi.smul_apply, smul_eq_mul]
  ring

/-- In D=1, a positive cone parameter is not a propagating mode: its amplitude
space is zero. This statement excludes every nonzero real frequency. -/
theorem dimension_one_no_propagating_mode {rho mu omega : ℝ} (k : Vec 1) {a : Vec 1}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (ha : a ≠ 0) :
    ¬ ActionStationary rho mu (planeWave omega k a) := by
  rw [actionStationary_planeWave_iff, modal_stationary_iff, dimension_one_operator]
  exact fun h => (mul_ne_zero hrho (pow_ne_zero 2 hw)) ((smul_eq_zero.mp h).resolve_right ha)

/-- An axis-aligned five-dimensional wave, with a non-unit speed, is stationary
against every admissible compact variation. -/
theorem concrete_five_transverse_stationary :
    ActionStationary 2 8 (planeWave 2 ![0, 0, 0, 0, 1] ![1, 0, 0, 0, 0]) := by
  rw [actionStationary_planeWave_iff, modal_stationary_iff]
  ext i
  fin_cases i <;> norm_num [modalOperator, normSq, dot, Fin.sum_univ_succ]

theorem concrete_five_longitudinal_detected :
    ∃ h : Point 5 → Vec 5, TestField h ∧
      deriv (relativeAction 2 8 (planeWave 2 ![0, 0, 0, 0, 1] ![0, 0, 0, 0, 1]) h) 0 ≠ 0 := by
  have hn : ¬ ActionStationary 2 8 (planeWave 2 ![0, 0, 0, 0, 1] ![0, 0, 0, 0, 1]) := by
    intro h
    have ht := stationary_transverse (by norm_num : (2 : ℝ) ≠ 0)
      (by norm_num : (2 : ℝ) ≠ 0) ((actionStationary_planeWave_iff _ _ _ _ _).mp h)
    norm_num [dot, Fin.sum_univ_succ] at ht
  simpa [ActionStationary] using hn

end
end S10Pilot
