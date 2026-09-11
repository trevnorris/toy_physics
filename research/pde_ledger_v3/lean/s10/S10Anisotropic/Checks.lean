import S10Anisotropic.Certificate

/-! Nonvacuity, explicit oblique/perpendicular examples, and the D=1 boundary. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section

theorem extra_wave_exists {D : ℕ} {e : Fin D} {sigma rho mu : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hrho : 0 < rho) (hmu : 0 < mu)
    (hk : k ≠ 0) (hq : perpSq e k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ extraVector e sigma k ≠ 0 ∧
      ActionStationary e sigma rho mu (planeWave omega k (extraVector e sigma k)) := by
  obtain ⟨_, w, _, hw, _, _, _, _, he, _, _⟩ :=
    split_variational_certificate hs hs1 hrho hmu hk hq
  exact ⟨w, hw, extraVector_ne_zero hq,
    (he _).mpr (Submodule.mem_span_singleton.mpr ⟨1, one_smul _ _⟩)⟩

theorem ordinary_wave_exists {D : ℕ} {e : Fin D} {sigma rho mu : ℝ} {k : Vec D}
    (hD : 3 ≤ D) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hrho : 0 < rho) (hmu : 0 < mu)
    (hk : k ≠ 0) (hq : perpSq e k ≠ 0) :
    ∃ omega : ℝ, ∃ a : Vec D, 0 < omega ∧ a ≠ 0 ∧
      a ∈ ordinarySpace e k ∧ ActionStationary e sigma rho mu (planeWave omega k a) := by
  obtain ⟨w, _, hw, _, _, _, _, ho, _, hd, _⟩ :=
    split_variational_certificate hs hs1 hrho hmu hk hq
  have hp : 0 < Module.finrank ℝ (ordinarySpace e k) := by rw [hd]; omega
  obtain ⟨a, ha⟩ := Module.finrank_pos_iff_exists_ne_zero.mp hp
  refine ⟨w, a, hw, ?_, a.property, (ho a).mpr a.property⟩
  intro hz
  exact ha (Subtype.ext hz)

theorem concrete_oblique_wave :
    ActionStationary (0 : Fin 3) 2 1 1
      (planeWave (Real.sqrt (3 / 2)) ![1, 1, 0] ![1, -2, 0]) ∧
    dot (![1, 1, 0] : Vec 3) ![1, -2, 0] = -1 := by
  constructor
  · rw [actionStationary_planeWave_iff, modal_stationary_normalized (by norm_num)]
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 3 / 2)]
    ext i
    fin_cases i <;>
      norm_num [normalizedOperator, unit, normSq, dot, Fin.sum_univ_succ]
  · norm_num [dot, Fin.sum_univ_succ]

/-- The extra polarization in the perpendicular stratum is a nonzero,
exactly transverse stationary field at the distinct extra frequency. -/
theorem concrete_perpendicular_wave :
    ActionStationary (0 : Fin 3) 2 1 1
      (planeWave (Real.sqrt (1 / 2)) ![0, 1, 0] ![1, 0, 0]) ∧
    dot (![0, 1, 0] : Vec 3) ![1, 0, 0] = 0 ∧
    extraValue (0 : Fin 3) 2 ![0, 1, 0] = 1 / 2 ∧
    normSq (![0, 1, 0] : Vec 3) = 1 := by
  constructor
  · rw [actionStationary_planeWave_iff, modal_stationary_normalized (by norm_num)]
    ext i
    fin_cases i <;>
      norm_num [normalizedOperator, unit, normSq, dot, Fin.sum_univ_succ]
  · norm_num [dot, extraValue, extraNumerator, perpSq, normSq, Fin.sum_univ_succ]

theorem concrete_perpendicular_transverse_dimension :
    Module.finrank ℝ (modeSpace (0 : Fin 3) 2 (1 / 2) ![0, 1, 0] ⊓
      transverseSpace ![0, 1, 0] : Submodule ℝ (Vec 3)) = 1 := by
  have hk : (![0, 1, 0] : Vec 3) ≠ 0 := by
    intro h
    have hi := congrFun h 1
    norm_num at hi
  have h := perpendicular_counts (e := (0 : Fin 3)) (sigma := 2)
    (by norm_num) (by norm_num) hk (by norm_num)
  rw [concrete_perpendicular_wave.2.2.1] at h
  exact h.2.2.2

theorem concrete_one_axis_kinetic :
    kinetic (0 : Fin 2) 2 ![0, 1] = 1 ∧ kinetic (0 : Fin 2) 2 ![1, 0] = 2 := by
  norm_num [kinetic, Fin.sum_univ_succ]

/-- There is no curl-restored propagation in one dimension even after changing
the only inertia coefficient. -/
theorem dimension_one_no_propagating_mode {sigma rho mu omega : ℝ} {k a : Vec 1}
    (hs : sigma ≠ 0) (hrho : rho ≠ 0) (hmu : mu ≠ 0)
    (hw : omega ≠ 0) (hk : k ≠ 0) (ha : a ≠ 0) :
    ¬ ActionStationary (0 : Fin 1) sigma rho mu (planeWave omega k a) := by
  have hq : perpSq (0 : Fin 1) k = 0 := by simp [perpSq, normSq, dot, sq]
  have hz : rho * omega ^ 2 / mu ≠ 0 := div_ne_zero (mul_ne_zero hrho (pow_ne_zero 2 hw)) hmu
  intro h
  rw [actionStationary_planeWave_iff, modal_stationary_normalized hmu,
    parallel_kernel_iff hs hz hk hq] at h
  rcases h with h | ⟨_, ht⟩
  · exact ha h
  · have he := ((mem_ordinarySpace (0 : Fin 1) k a).mp ht).1
    apply ha
    ext i
    have hi : i = 0 := Subsingleton.elim _ _
    subst i
    exact he

end
end S10Anisotropic
