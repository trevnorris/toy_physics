import S10Pilot.Action

/-! Spectrum of the action-derived S10 operator for every finite spatial dimension. -/

namespace S10Pilot
noncomputable section

variable {D : ℕ}

theorem dot_modalOperator (rho mu omega : ℝ) (k a : Vec D) :
    dot k (modalOperator rho mu omega k a) = rho * omega ^ 2 * dot k a := by
  rw [modalOperator_eq, dot_add_right, dot_smul_right, dot_smul_right]
  unfold normSq
  ring

theorem transverse_operator (rho mu omega : ℝ) (k a : Vec D) (ha : dot k a = 0) :
    modalOperator rho mu omega k a = (rho * omega ^ 2 - mu * normSq k) • a := by
  ext i
  simp [modalOperator, ha]
  ring

theorem longitudinal_operator (rho mu omega s : ℝ) (k : Vec D) :
    modalOperator rho mu omega k (s • k) = (rho * omega ^ 2 * s) • k := by
  ext i
  simp [modalOperator, dot_smul_right, normSq]
  ring

theorem stationary_transverse {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (homega : omega ≠ 0)
    (h : ModalStationary rho mu omega k a) : dot k a = 0 := by
  have hm := (modal_stationary_iff _ _ _ _ _).mp h
  have hd := dot_modalOperator rho mu omega k a
  rw [hm] at hd
  rw [dot_zero] at hd
  exact (mul_eq_zero.mp hd.symm).resolve_left (mul_ne_zero hrho (pow_ne_zero 2 homega))

/-- Completeness: every nonzero-frequency, nonzero-amplitude mode is transverse
and lies on the derived cone. No nonzero-frequency longitudinal solution remains. -/
theorem propagating_mode_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (homega : omega ≠ 0) (ha : a ≠ 0) :
    ModalStationary rho mu omega k a ↔
      dot k a = 0 ∧ omega ^ 2 = (mu / rho) * normSq k := by
  constructor
  · intro h
    have ht := stationary_transverse hrho homega h
    refine ⟨ht, ?_⟩
    have hm := (modal_stationary_iff _ _ _ _ _).mp h
    rw [transverse_operator _ _ _ _ _ ht] at hm
    have hc : rho * omega ^ 2 - mu * normSq k = 0 :=
      (smul_eq_zero.mp hm).resolve_right ha
    field_simp
    nlinarith
  · rintro ⟨ht, hf⟩
    rw [modal_stationary_iff, transverse_operator _ _ _ _ _ ht, hf]
    have hc : rho * (mu / rho * normSq k) - mu * normSq k = 0 := by
      field_simp
      ring
    rw [hc, zero_smul]

def dotLinear (k : Vec D) : Vec D →ₗ[ℝ] ℝ where
  toFun := dot k
  map_add' a b := dot_add_right k a b
  map_smul' s a := dot_smul_right k a s

def transverseSpace (k : Vec D) : Submodule ℝ (Vec D) := (dotLinear k).ker
def longitudinalSpace (k : Vec D) : Submodule ℝ (Vec D) := Submodule.span ℝ {k}

theorem transverseSpace_finrank {k : Vec D} (hk : k ≠ 0) :
    Module.finrank ℝ (transverseSpace k) = D - 1 := by
  have hn : normSq k ≠ 0 := ne_of_gt (dot_self_pos hk)
  have hs : Function.Surjective (dotLinear k) := by
    intro r
    refine ⟨(r / normSq k) • k, ?_⟩
    change dot k ((r / normSq k) • k) = r
    have hd : dot k ((r / normSq k) • k) = (r / normSq k) * normSq k := by
      exact dot_smul_right k k _
    rw [hd, div_mul_cancel₀ _ hn]
  have h := (dotLinear k).finrank_range_add_finrank_ker
  rw [LinearMap.range_eq_top.mpr hs] at h
  simp at h
  unfold transverseSpace
  omega

theorem longitudinalSpace_finrank {k : Vec D} (hk : k ≠ 0) :
    Module.finrank ℝ (longitudinalSpace k) = 1 :=
  finrank_span_singleton hk

theorem stationary_on_cone_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (homega : omega ≠ 0)
    (hf : omega ^ 2 = (mu / rho) * normSq k) :
    ModalStationary rho mu omega k a ↔ a ∈ transverseSpace k := by
  change ModalStationary rho mu omega k a ↔ dot k a = 0
  constructor
  · exact stationary_transverse hrho homega
  · intro ht
    by_cases ha : a = 0
    · rw [ha, modal_stationary_iff]
      ext i
      simp [modalOperator, dot]
    · exact (propagating_mode_iff hrho homega ha).mpr ⟨ht, hf⟩

/-- At zero frequency, the entire nullspace is precisely the longitudinal span. -/
theorem zero_frequency_iff {rho mu : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) (hk : k ≠ 0) :
    ModalStationary rho mu 0 k a ↔ a ∈ longitudinalSpace k := by
  rw [modal_stationary_iff]
  constructor
  · intro h
    have hn : normSq k ≠ 0 := ne_of_gt (dot_self_pos hk)
    have heq : (dot k a / normSq k) • k = a := by
      ext i
      have hi := congrFun h i
      simp only [modalOperator, zero_pow (by decide : 2 ≠ 0), mul_zero, zero_mul,
        zero_sub, Pi.zero_apply, neg_eq_zero, mul_eq_zero] at hi
      have hz := hi.resolve_left hmu
      simp only [Pi.smul_apply, smul_eq_mul]
      field_simp
      nlinarith
    exact Submodule.mem_span_singleton.mpr ⟨dot k a / normSq k, heq⟩
  · intro h
    obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp h
    simp [longitudinal_operator]

def coneValue (rho mu : ℝ) (k : Vec D) : ℝ := (mu / rho) * normSq k

theorem coneValue_pos {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) : 0 < coneValue rho mu k :=
  mul_pos (div_pos hmu hrho) (dot_self_pos hk)

theorem coneValue_scaling (rho mu s : ℝ) (k : Vec D) :
    coneValue rho mu (s • k) = s ^ 2 * coneValue rho mu k := by
  simp [coneValue, normSq, dot_smul_left, dot_smul_right]
  ring

end
end S10Pilot
