import S10Anisotropic.Census
import S10Anisotropic.FiniteAction
import S10Anisotropic.PhaseAverage

/-! The spectral census is tied to the integrated action, including arbitrary
smooth compact variations outside the plane-wave ansatz. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

def extraConeValue (e : Fin D) (sigma rho mu : ℝ) (k : Vec D) : ℝ :=
  (mu / rho) * extraValue e sigma k

theorem normalized_frequency_iff {rho mu omega z : ℝ} (hrho : rho ≠ 0) (hmu : mu ≠ 0) :
    rho * omega ^ 2 / mu = z ↔ omega ^ 2 = (mu / rho) * z := by
  constructor <;> intro h <;> field_simp at h ⊢ <;> nlinarith

theorem extra_frequency_iff {e : Fin D} {sigma z : ℝ} {k : Vec D} (hs : sigma ≠ 0) :
    sigma * z = extraNumerator e sigma k ↔ z = extraValue e sigma k := by
  unfold extraValue
  rw [eq_div_iff hs, mul_comm z sigma]

/-- This treats z as an arbitrary real squared-frequency parameter, not as an
already nonnegative square. Every nonzero real root is positive. -/
theorem nonzero_root_positive {e : Fin D} {sigma z : ℝ} {k a : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hz : z ≠ 0) (hk : k ≠ 0) (ha : a ≠ 0)
    (h : normalizedOperator e sigma z k a = 0) : 0 < z := by
  by_cases hq : perpSq e k = 0
  · have hf := ((parallel_kernel_iff (ne_of_gt hs) hz hk hq).mp h).resolve_left ha
    rw [hf.1]
    exact dot_self_pos hk
  · rcases (split_propagating_iff hs1 hq hz ha).mp h with ho | he
    · rw [ho.1]
      exact dot_self_pos hk
    · rw [(extra_frequency_iff (ne_of_gt hs)).mp he.1]
      exact extraValue_pos hs hk

theorem actionStationary_iff_mem_modeSpace {e : Fin D} {sigma rho mu omega : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) :
    ActionStationary e sigma rho mu (planeWave omega k a) ↔
      a ∈ modeSpace e sigma (rho * omega ^ 2 / mu) k := by
  rw [actionStationary_planeWave_iff, modal_stationary_normalized hmu, mem_modeSpace]

/-- Complete nonzero-frequency classification in physical frequency variables. -/
theorem split_variational_iff {e : Fin D} {sigma rho mu omega : ℝ} {k a : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hrho : 0 < rho) (hmu : 0 < mu)
    (hq : perpSq e k ≠ 0) (hw : omega ≠ 0) (ha : a ≠ 0) :
    ActionStationary e sigma rho mu (planeWave omega k a) ↔
      (omega ^ 2 = coneValue rho mu k ∧ a ∈ ordinarySpace e k) ∨
      (omega ^ 2 = extraConeValue e sigma rho mu k ∧
        a ∈ longitudinalSpace (extraVector e sigma k)) := by
  have hz : rho * omega ^ 2 / mu ≠ 0 :=
    div_ne_zero (mul_ne_zero (ne_of_gt hrho) (pow_ne_zero 2 hw)) (ne_of_gt hmu)
  rw [actionStationary_planeWave_iff, modal_stationary_normalized (ne_of_gt hmu),
    split_propagating_iff hs1 hq hz ha, extra_frequency_iff (ne_of_gt hs),
    normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu),
    normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu)]
  rfl

theorem parallel_variational_iff {e : Fin D} {sigma rho mu omega : ℝ} {k a : Vec D}
    (hs : 0 < sigma) (hrho : 0 < rho) (hmu : 0 < mu)
    (hk : k ≠ 0) (hq : perpSq e k = 0) (hw : omega ≠ 0) :
    ActionStationary e sigma rho mu (planeWave omega k a) ↔
      a = 0 ∨ (omega ^ 2 = coneValue rho mu k ∧ a ∈ ordinarySpace e k) := by
  have hz : rho * omega ^ 2 / mu ≠ 0 :=
    div_ne_zero (mul_ne_zero (ne_of_gt hrho) (pow_ne_zero 2 hw)) (ne_of_gt hmu)
  rw [actionStationary_planeWave_iff, modal_stationary_normalized (ne_of_gt hmu),
    parallel_kernel_iff (ne_of_gt hs) hz hk hq,
    normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu)]
  rfl

theorem zero_variational_iff {e : Fin D} {sigma rho mu : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) (hk : k ≠ 0) :
    ActionStationary e sigma rho mu (planeWave 0 k a) ↔ a ∈ longitudinalSpace k := by
  rw [actionStationary_planeWave_iff, zero_frequency_iff hmu hk]

theorem positive_frequency_exists {rho mu z : ℝ} (hrho : 0 < rho) (hmu : 0 < mu) (hz : 0 < z) :
    ∃ omega : ℝ, 0 < omega ∧ rho * omega ^ 2 / mu = z := by
  refine ⟨Real.sqrt ((mu / rho) * z), Real.sqrt_pos.mpr (mul_pos (div_pos hmu hrho) hz), ?_⟩
  rw [normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu)]
  exact Real.sq_sqrt (le_of_lt (mul_pos (div_pos hmu hrho) hz))

/-- Positive, distinct squared frequencies, full amplitude spaces, and their
dimensions, all obtained from the integrated variational principle. -/
theorem split_variational_certificate {e : Fin D} {sigma rho mu : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hrho : 0 < rho) (hmu : 0 < mu)
    (hk : k ≠ 0) (hq : perpSq e k ≠ 0) :
    ∃ wo we : ℝ, 0 < wo ∧ 0 < we ∧ wo ^ 2 ≠ we ^ 2 ∧
      wo ^ 2 = coneValue rho mu k ∧ we ^ 2 = extraConeValue e sigma rho mu k ∧
      (∀ a, ActionStationary e sigma rho mu (planeWave wo k a) ↔ a ∈ ordinarySpace e k) ∧
      (∀ a, ActionStationary e sigma rho mu (planeWave we k a) ↔
        a ∈ longitudinalSpace (extraVector e sigma k)) ∧
      Module.finrank ℝ (ordinarySpace e k) = D - 2 ∧
      Module.finrank ℝ (longitudinalSpace (extraVector e sigma k)) = 1 := by
  obtain ⟨wo, hwo, hfo⟩ := positive_frequency_exists hrho hmu (dot_self_pos hk)
  obtain ⟨we, hwe, hfe⟩ := positive_frequency_exists hrho hmu (extraValue_pos hs hk)
  have hne : wo ^ 2 ≠ we ^ 2 := by
    intro h
    rw [h] at hfo
    rw [hfe] at hfo
    exact extraValue_ne_ordinary (ne_of_gt hs) hs1 hq hfo
  refine ⟨wo, we, hwo, hwe, hne,
    (normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu)).mp hfo,
    (normalized_frequency_iff (ne_of_gt hrho) (ne_of_gt hmu)).mp hfe, ?_, ?_,
    ordinarySpace_finrank hq, longitudinalSpace_finrank (extraVector_ne_zero hq)⟩
  · intro a
    rw [actionStationary_iff_mem_modeSpace (ne_of_gt hmu), hfo, ordinary_modeSpace hs1 hq hk]
  · intro a
    rw [actionStationary_iff_mem_modeSpace (ne_of_gt hmu), hfe, extra_modeSpace hs hs1 hq hk]

theorem frequency_coincidence_iff {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    extraValue e sigma k = normSq k ↔ perpSq e k = 0 := by
  constructor
  · intro h
    by_contra hq
    exact extraValue_ne_ordinary hs hs1 hq h
  · intro hq
    have hn : normSq k = k e ^ 2 := sub_eq_zero.mp hq
    unfold extraValue extraNumerator
    rw [hq, zero_add, hn, mul_div_cancel_left₀ _ hs]

/-- Exact transverse membership of the extra polarization: the perpendicular
case is an allowed exceptional stratum even though no frequencies coincide. -/
theorem extra_exactly_transverse_iff {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) :
    extraVector e sigma k ∈ transverseSpace k ↔ k e = 0 := by
  change dot k (extraVector e sigma k) = 0 ↔ _
  rw [dot_extraVector]
  simp [mul_eq_zero, sub_ne_zero.mpr hs1.symm, hq]

theorem isotropic_density (e : Fin D) (rho mu : ℝ) (J : Jet D) :
    lagrangian e 1 rho mu J = S10Pilot.lagrangian rho mu J := by
  simp [lagrangian_eq]

theorem isotropic_variational_iff (e : Fin D) (rho mu : ℝ) (u : Point D → Vec D) :
    ActionStationary e 1 rho mu u ↔ S10Pilot.ActionStationary rho mu u := by
  have heq : ∀ h, relativeAction e 1 rho mu u h = S10Pilot.relativeAction rho mu u h := by
    intro h
    funext s
    simp only [relativeAction, S10Pilot.relativeAction, isotropic_density]
  simp only [ActionStationary, S10Pilot.ActionStationary, heq]

end
end S10Anisotropic
