import S11Homogeneous.Action
import S10Pilot.Spectrum

/-! Entire amplitude spaces of the combined operator. The formulas include
coincident roots; geometrical polarizations are not inferred from an arbitrary
basis of a degenerate eigenspace. -/

namespace S11Homogeneous
open S10Pilot
noncomputable section
variable {D : ℕ}

theorem coefficient_zero_iff {rho c omega : ℝ} (hrho : rho ≠ 0) (k : Vec D) :
    rho * omega ^ 2 - c * normSq k = 0 ↔ omega ^ 2 = coneValue rho c k := by
  unfold coneValue
  constructor <;> intro h <;> field_simp at * <;> nlinarith

theorem dot_modalOperator (rho mu B omega : ℝ) (k a : Vec D) :
    dot k (modalOperator rho mu B omega k a) =
      (rho * omega ^ 2 - B * normSq k) * dot k a := by
  rw [modalOperator, dot_add_right, dot_smul_right, dot_smul_right]
  unfold normSq
  ring

theorem transverse_operator (rho mu B omega : ℝ) (k a : Vec D) (ha : dot k a = 0) :
    modalOperator rho mu B omega k a = (rho * omega ^ 2 - mu * normSq k) • a := by
  simp [modalOperator, ha]

theorem longitudinal_operator (rho mu B omega s : ℝ) (k : Vec D) :
    modalOperator rho mu B omega k (s • k) =
      ((rho * omega ^ 2 - B * normSq k) * s) • k := by
  ext i
  simp [modalOperator, dot_smul_right, normSq]
  ring

theorem frequency_coincidence_iff {rho mu B : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hk : k ≠ 0) :
    coneValue rho mu k = coneValue rho B k ↔ B = mu := by
  have hn := ne_of_gt (dot_self_pos hk)
  unfold coneValue
  constructor
  · intro h
    have hdiv := mul_right_cancel₀ hn h
    field_simp at hdiv
    nlinarith
  · rintro rfl
    rfl

theorem transverse_kernel_iff {rho mu B omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hk : k ≠ 0) (hB : B ≠ mu)
    (hf : omega ^ 2 = coneValue rho mu k) :
    ModalStationary rho mu B omega k a ↔ a ∈ transverseSpace k := by
  rw [modal_stationary_iff]
  have hc := (coefficient_zero_iff hrho k).mpr hf
  change modalOperator rho mu B omega k a = 0 ↔ dot k a = 0
  simp only [modalOperator, hc, zero_smul, zero_add, smul_eq_zero, hk, or_false,
    mul_eq_zero, sub_eq_zero, Ne.symm hB, false_or]

theorem operator_on_longitudinal_cone {rho mu B omega : ℝ} (hrho : rho ≠ 0)
    (k a : Vec D) (hf : omega ^ 2 = coneValue rho B k) :
    modalOperator rho mu B omega k a = S10Pilot.modalOperator rho (mu - B) 0 k a := by
  have hc := (coefficient_zero_iff hrho k).mpr hf
  have hr : rho * omega ^ 2 = B * normSq k := sub_eq_zero.mp hc
  ext i
  simp only [modalOperator, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    S10Pilot.modalOperator, zero_pow (by decide : 2 ≠ 0), mul_zero]
  rw [hr]
  ring

theorem longitudinal_kernel_iff {rho mu B omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hk : k ≠ 0) (hB : B ≠ mu)
    (hf : omega ^ 2 = coneValue rho B k) :
    ModalStationary rho mu B omega k a ↔ a ∈ longitudinalSpace k := by
  rw [modal_stationary_iff, operator_on_longitudinal_cone hrho k a hf,
    ← S10Pilot.modal_stationary_iff]
  exact S10Pilot.zero_frequency_iff (sub_ne_zero.mpr (Ne.symm hB)) hk

theorem off_roots_iff {rho mu B omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (ht : omega ^ 2 ≠ coneValue rho mu k)
    (hl : omega ^ 2 ≠ coneValue rho B k) :
    ModalStationary rho mu B omega k a ↔ a = 0 := by
  rw [modal_stationary_iff]
  constructor
  · intro h
    have hd := dot_modalOperator rho mu B omega k a
    rw [h, dot_zero] at hd
    have hcoeff : rho * omega ^ 2 - B * normSq k ≠ 0 :=
      fun he => hl ((coefficient_zero_iff hrho k).mp he)
    have ha := (mul_eq_zero.mp hd.symm).resolve_left hcoeff
    rw [transverse_operator _ _ _ _ _ _ ha] at h
    exact (smul_eq_zero.mp h).resolve_left
      (fun he => ht ((coefficient_zero_iff hrho k).mp he))
  · rintro rfl
    simp [modalOperator]

def modalMap (rho mu B omega : ℝ) (k : Vec D) : Vec D →ₗ[ℝ] Vec D where
  toFun := modalOperator rho mu B omega k
  map_add' a b := by
    simp [modalOperator, dot_add_right, mul_add, add_smul, smul_add]
    abel
  map_smul' s a := by
    ext i
    simp [modalOperator, dot_smul_right]
    ring

def modeSpace (rho mu B omega : ℝ) (k : Vec D) : Submodule ℝ (Vec D) :=
  (modalMap rho mu B omega k).ker

theorem mem_modeSpace (rho mu B omega : ℝ) (k a : Vec D) :
    a ∈ modeSpace rho mu B omega k ↔ ModalStationary rho mu B omega k a := by
  rw [modal_stationary_iff]
  rfl

theorem coincidence_space {rho mu omega : ℝ} (hrho : rho ≠ 0) (k : Vec D)
    (hf : omega ^ 2 = coneValue rho mu k) :
    modeSpace rho mu mu omega k = ⊤ := by
  ext a
  simp only [mem_modeSpace, modal_stationary_iff, Submodule.mem_top, iff_true]
  simp [modalOperator, (coefficient_zero_iff hrho k).mpr hf]

/-- Exhaustive classification, including the merged root and all off-root frequencies. -/
theorem modeSpace_classification {rho mu B omega : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hk : k ≠ 0) :
    modeSpace rho mu B omega k =
      if omega ^ 2 = coneValue rho mu k then
        if B = mu then ⊤ else transverseSpace k
      else if omega ^ 2 = coneValue rho B k then longitudinalSpace k else ⊥ := by
  split_ifs with ht hb hl
  · subst B
    exact coincidence_space hrho k ht
  · ext a
    rw [mem_modeSpace, transverse_kernel_iff hrho hk hb ht]
  · have hb : B ≠ mu := by
      rintro rfl
      exact ht hl
    ext a
    rw [mem_modeSpace, longitudinal_kernel_iff hrho hk hb hl]
  · ext a
    rw [mem_modeSpace, off_roots_iff hrho ht hl]
    rfl

theorem kernel_census_three {rho mu B omega : ℝ} {k : Vec 3}
    (hrho : rho ≠ 0) (hk : k ≠ 0) :
    Module.finrank ℝ (modeSpace rho mu B omega k) =
      if omega ^ 2 = coneValue rho mu k then (if B = mu then 3 else 2)
      else if omega ^ 2 = coneValue rho B k then 1 else 0 := by
  have hs := modeSpace_classification (mu := mu) (B := B) (omega := omega) hrho hk
  by_cases ht : omega ^ 2 = coneValue rho mu k
  · by_cases hb : B = mu
    · rw [if_pos ht, if_pos hb] at hs
      rw [hs, if_pos ht, if_pos hb]
      simp
    · rw [if_pos ht, if_neg hb] at hs
      rw [hs, if_pos ht, if_neg hb]
      exact transverseSpace_finrank hk
  · by_cases hl : omega ^ 2 = coneValue rho B k
    · rw [if_neg ht, if_pos hl] at hs
      rw [hs, if_neg ht, if_pos hl]
      exact longitudinalSpace_finrank hk
    · rw [if_neg ht, if_neg hl] at hs
      rw [hs, if_neg ht, if_neg hl]
      simp

theorem positive_frequencies {rho mu B : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hB : 0 < B) (hk : k ≠ 0) :
    ∃ wt wl : ℝ, 0 < wt ∧ 0 < wl ∧ wt ^ 2 = coneValue rho mu k ∧
      wl ^ 2 = coneValue rho B k := by
  have ht := coneValue_pos hrho hmu hk
  have hl := coneValue_pos hrho hB hk
  exact ⟨Real.sqrt (coneValue rho mu k), Real.sqrt (coneValue rho B k),
    Real.sqrt_pos.mpr ht, Real.sqrt_pos.mpr hl,
    Real.sq_sqrt (le_of_lt ht), Real.sq_sqrt (le_of_lt hl)⟩

theorem zero_compression_stationarity (rho mu : ℝ) (u : Point D → Vec D) :
    ActionStationary rho mu 0 u ↔ S10Pilot.ActionStationary rho mu u := by
  have he : ∀ h, relativeAction rho mu 0 u h = S10Pilot.relativeAction rho mu u h := by
    intro h
    funext s
    unfold relativeAction S10Pilot.relativeAction
    simp only [zero_compression_action]
  simp only [ActionStationary, S10Pilot.ActionStationary, he]

theorem zero_wavevector_iff {rho mu B omega : ℝ} (hrho : rho ≠ 0) (a : Vec D) :
    ModalStationary rho mu B omega 0 a ↔ omega = 0 ∨ a = 0 := by
  rw [modal_stationary_iff]
  simp [modalOperator, normSq, smul_eq_zero, mul_eq_zero, hrho]

end
end S11Homogeneous
