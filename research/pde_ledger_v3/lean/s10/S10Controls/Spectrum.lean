import S10Controls.Action
import S10Pilot.Spectrum

/-! Complete real plane-wave amplitude spaces for the two stiffness controls. -/

namespace S10Controls
open S10Pilot
noncomputable section
variable {D : ℕ}

theorem full_stationary_iff {rho mu omega : ℝ} {k a : Vec D} (hrho : rho ≠ 0) :
    ModalStationary .fullGradient rho mu omega k a ↔
      a = 0 ∨ omega ^ 2 = coneValue rho mu k := by
  rw [modal_stationary_iff, modalOperator, smul_eq_zero]
  constructor
  · rintro (hc | ha)
    · right
      unfold coneValue
      field_simp
      nlinarith
    · exact Or.inl ha
  · rintro (ha | hf)
    · exact Or.inr ha
    · left
      rw [hf]
      unfold coneValue
      field_simp
      ring

theorem div_longitudinal_operator (rho mu omega s : ℝ) (k : Vec D) :
    modalOperator .divergenceOnly rho mu omega k (s • k) =
      ((rho * omega ^ 2 - mu * normSq k) * s) • k := by
  ext i
  simp [modalOperator, dot_smul_right, normSq]
  ring

theorem div_stationary_longitudinal {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0)
    (h : ModalStationary .divergenceOnly rho mu omega k a) : a ∈ longitudinalSpace k := by
  have hm := (modal_stationary_iff _ _ _ _ _ _).mp h
  have hc : rho * omega ^ 2 ≠ 0 := mul_ne_zero hrho (pow_ne_zero 2 hw)
  have heq : ((mu * dot k a) / (rho * omega ^ 2)) • k = a := by
    ext i
    have hi := congrFun hm i
    simp only [modalOperator, Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at hi
    simp only [Pi.smul_apply, smul_eq_mul]
    field_simp
    nlinarith
  exact Submodule.mem_span_singleton.mpr ⟨_, heq⟩

theorem div_propagating_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0) (ha : a ≠ 0) :
    ModalStationary .divergenceOnly rho mu omega k a ↔
      a ∈ longitudinalSpace k ∧ omega ^ 2 = coneValue rho mu k := by
  constructor
  · intro h
    have hl := div_stationary_longitudinal hrho hw h
    refine ⟨hl, ?_⟩
    obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp hl
    have hs : s ≠ 0 := fun hzero => ha (by simp [hzero])
    have hm := (modal_stationary_iff _ _ _ _ _ _).mp h
    rw [div_longitudinal_operator] at hm
    have hc := (mul_eq_zero.mp ((smul_eq_zero.mp hm).resolve_right hk)).resolve_right hs
    unfold coneValue
    field_simp
    nlinarith
  · rintro ⟨hl, hf⟩
    obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp hl
    rw [modal_stationary_iff, div_longitudinal_operator, hf]
    have hc : rho * coneValue rho mu k - mu * normSq k = 0 := by
      unfold coneValue
      field_simp
      ring
    rw [hc, zero_mul, zero_smul]

theorem div_on_cone_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0)
    (hf : omega ^ 2 = coneValue rho mu k) :
    ModalStationary .divergenceOnly rho mu omega k a ↔ a ∈ longitudinalSpace k := by
  by_cases ha : a = 0
  · simp [ha, modal_stationary_iff, modalOperator]
  · rw [div_propagating_iff hrho hw hk ha]
    simp [hf]

theorem div_zero_frequency_iff {rho mu : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) (hk : k ≠ 0) :
    ModalStationary .divergenceOnly rho mu 0 k a ↔ a ∈ transverseSpace k := by
  change ModalStationary .divergenceOnly rho mu 0 k a ↔ dot k a = 0
  simp [modal_stationary_iff, modalOperator, hk, hmu]

def modalMap (form : Form) (rho mu omega : ℝ) (k : Vec D) : Vec D →ₗ[ℝ] Vec D where
  toFun := modalOperator form rho mu omega k
  map_add' a b := by
    cases form <;> ext i <;> simp [modalOperator, dot_add_right] <;> ring
  map_smul' s a := by
    cases form <;> ext i <;> simp [modalOperator, dot_smul_right] <;> ring

def amplitudeSpace (form : Form) (rho mu omega : ℝ) (k : Vec D) : Submodule ℝ (Vec D) :=
  (modalMap form rho mu omega k).ker

theorem mem_amplitudeSpace (form : Form) (rho mu omega : ℝ) (k a : Vec D) :
    a ∈ amplitudeSpace form rho mu omega k ↔ ModalStationary form rho mu omega k a := by
  change modalOperator form rho mu omega k a = 0 ↔ _
  exact (modal_stationary_iff _ _ _ _ _ _).symm

theorem full_on_cone_space {rho mu omega : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hf : omega ^ 2 = coneValue rho mu k) :
    amplitudeSpace .fullGradient rho mu omega k = ⊤ := by
  ext a
  rw [mem_amplitudeSpace, full_stationary_iff hrho]
  simp [hf]

theorem full_zero_space {rho mu : ℝ} {k : Vec D} (hmu : mu ≠ 0) (hk : k ≠ 0) :
    amplitudeSpace .fullGradient rho mu 0 k = ⊥ := by
  have hn : normSq k ≠ 0 := ne_of_gt (dot_self_pos hk)
  ext a
  rw [mem_amplitudeSpace, modal_stationary_iff]
  simp [modalOperator, hn, hmu]

theorem div_on_cone_space {rho mu omega : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0)
    (hf : omega ^ 2 = coneValue rho mu k) :
    amplitudeSpace .divergenceOnly rho mu omega k = longitudinalSpace k := by
  ext a
  rw [mem_amplitudeSpace, div_on_cone_iff hrho hw hk hf]

theorem div_zero_space {rho mu : ℝ} {k : Vec D} (hmu : mu ≠ 0) (hk : k ≠ 0) :
    amplitudeSpace .divergenceOnly rho mu 0 k = transverseSpace k := by
  ext a
  rw [mem_amplitudeSpace, div_zero_frequency_iff hmu hk]

theorem longitudinal_inf_transverse {k : Vec D} (hk : k ≠ 0) :
    longitudinalSpace k ⊓ transverseSpace k = ⊥ := by
  have hn : normSq k ≠ 0 := ne_of_gt (dot_self_pos hk)
  ext a
  change (a ∈ longitudinalSpace k ∧ a ∈ transverseSpace k) ↔ a = 0
  constructor
  · rintro ⟨hl, ht⟩
    obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp hl
    change dot k (s • k) = 0 at ht
    rw [dot_smul_right] at ht
    have hs : s = 0 := (mul_eq_zero.mp ht).resolve_right hn
    simp [hs]
  · intro h
    simp [h]

end
end S10Controls
