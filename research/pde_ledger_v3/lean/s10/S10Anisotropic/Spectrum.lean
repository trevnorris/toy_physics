import S10Anisotropic.Geometry

/-! Exhaustive real amplitude classification. z denotes rho*omega^2/mu.
No division by a wavevector component is made in the action or operator. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

def normalizedOperator (e : Fin D) (sigma z : ℝ) (k a : Vec D) : Vec D :=
  (z - normSq k) • a + dot k a • k + (z * (sigma - 1) * a e) • unit e

theorem modalOperator_normalized {e : Fin D} {sigma rho mu omega : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) :
    modalOperator e sigma rho mu omega k a =
      mu • normalizedOperator e sigma (rho * omega ^ 2 / mu) k a := by
  ext i
  simp only [modalOperator, S10Pilot.modalOperator, normalizedOperator,
    Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  field_simp
  ring

theorem modal_stationary_normalized {e : Fin D} {sigma rho mu omega : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) :
    ModalStationary e sigma rho mu omega k a ↔
      normalizedOperator e sigma (rho * omega ^ 2 / mu) k a = 0 := by
  rw [modal_stationary_iff, modalOperator_normalized hmu, smul_eq_zero]
  simp [hmu]

theorem normalizedOperator_smul (e : Fin D) (sigma z c : ℝ) (k a : Vec D) :
    normalizedOperator e sigma z k (c • a) = c • normalizedOperator e sigma z k a := by
  ext i
  simp [normalizedOperator, dot_smul_right]
  ring

theorem weighted_transverse {e : Fin D} {sigma z : ℝ} {k a : Vec D} (hz : z ≠ 0)
    (h : normalizedOperator e sigma z k a = 0) :
    dot k a + (sigma - 1) * k e * a e = 0 := by
  have hd := congrArg (dot k) h
  simp only [normalizedOperator, dot_add_right, dot_smul_right, dot_unit_right, dot_zero] at hd
  have he : z * (dot k a + (sigma - 1) * k e * a e) = 0 := by
    change _ + _ * normSq k + _ = 0 at hd
    nlinarith
  exact (mul_eq_zero.mp he).resolve_left hz

theorem axis_factor {e : Fin D} {sigma z : ℝ} {k a : Vec D} (hz : z ≠ 0)
    (h : normalizedOperator e sigma z k a = 0) :
    (sigma * z - extraNumerator e sigma k) * a e = 0 := by
  have hd := weighted_transverse hz h
  have hi := congrFun h e
  simp only [normalizedOperator, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    unit_self, mul_one, Pi.zero_apply] at hi
  unfold extraNumerator perpSq
  linear_combination hi - k e * hd

theorem ordinary_operator {e : Fin D} {sigma z : ℝ} {k a : Vec D}
    (ha : a ∈ ordinarySpace e k) :
    normalizedOperator e sigma z k a = (z - normSq k) • a := by
  obtain ⟨he, hd⟩ := (mem_ordinarySpace e k a).mp ha
  simp [normalizedOperator, he, hd]

theorem extra_operator (e : Fin D) (sigma z : ℝ) (k : Vec D) :
    normalizedOperator e sigma z k (extraVector e sigma k) =
      (sigma * z - extraNumerator e sigma k) • (normSq k • unit e - k e • k) := by
  unfold normalizedOperator
  rw [dot_extraVector, extraVector_axis]
  ext i
  simp only [extraVector, extraNumerator, perpSq, Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem off_ordinary_reconstruction {e : Fin D} {sigma z : ℝ} {k a : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hz : z ≠ 0)
    (hf : sigma * z = extraNumerator e sigma k)
    (h : normalizedOperator e sigma z k a = 0) :
    a = (a e / perpSq e k) • extraVector e sigma k := by
  have hd := weighted_transverse hz h
  have hf0 : sigma * z - extraNumerator e sigma k = 0 := sub_eq_zero.mpr hf
  ext i
  have hi := congrFun h i
  simp only [normalizedOperator, Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at hi
  have hp : (sigma - 1) * (perpSq e k * a i - a e *
      (extraNumerator e sigma k * unit e i - sigma * k e * k i)) = 0 := by
    unfold extraNumerator perpSq at hf0 ⊢
    linear_combination -sigma * hi +
      (a i + (sigma - 1) * a e * unit e i) * hf0 + sigma * k i * hd
  have he := (mul_eq_zero.mp hp).resolve_left (sub_ne_zero.mpr hs1)
  simp only [extraVector, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  field_simp
  nlinarith

/-- For a nonzero amplitude at nonzero frequency, these are all possibilities
whenever the wavevector has a nonzero component perpendicular to the axis. -/
theorem split_propagating_iff {e : Fin D} {sigma z : ℝ} {k a : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hz : z ≠ 0) (ha : a ≠ 0) :
    normalizedOperator e sigma z k a = 0 ↔
      (z = normSq k ∧ a ∈ ordinarySpace e k) ∨
      (sigma * z = extraNumerator e sigma k ∧
        a ∈ longitudinalSpace (extraVector e sigma k)) := by
  constructor
  · intro h
    have hd := weighted_transverse hz h
    have hf := axis_factor hz h
    by_cases hc : z = normSq k
    · left
      refine ⟨hc, (mem_ordinarySpace e k a).mpr ?_⟩
      have hp : (sigma - 1) * perpSq e k * a e = 0 := by
        rw [hc] at hf
        unfold extraNumerator perpSq at hf
        unfold perpSq
        nlinarith
      have he := (mul_eq_zero.mp hp).resolve_left (mul_ne_zero (sub_ne_zero.mpr hs1) hq)
      exact ⟨he, by simpa [he] using hd⟩
    · right
      have he : a e ≠ 0 := by
        intro he
        have hda : dot k a = 0 := by simpa [he] using hd
        rw [ordinary_operator ((mem_ordinarySpace e k a).mpr ⟨he, hda⟩)] at h
        exact (smul_ne_zero (sub_ne_zero.mpr hc) ha) h
      have hfreq : sigma * z = extraNumerator e sigma k :=
        sub_eq_zero.mp ((mul_eq_zero.mp hf).resolve_right he)
      exact ⟨hfreq, Submodule.mem_span_singleton.mpr
        ⟨_, (off_ordinary_reconstruction hs1 hq hz hfreq h).symm⟩⟩
  · rintro (⟨hf, ht⟩ | ⟨hf, hl⟩)
    · rw [ordinary_operator ht, hf, sub_self, zero_smul]
    · obtain ⟨c, rfl⟩ := Submodule.mem_span_singleton.mp hl
      rw [normalizedOperator_smul, extra_operator, hf, sub_self, zero_smul, smul_zero]

theorem ordinary_kernel_iff {e : Fin D} {sigma : ℝ} {k a : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hk : k ≠ 0) :
    normalizedOperator e sigma (normSq k) k a = 0 ↔ a ∈ ordinarySpace e k := by
  constructor
  · intro h
    have hd := weighted_transverse (ne_of_gt (dot_self_pos hk)) h
    have hf := axis_factor (ne_of_gt (dot_self_pos hk)) h
    have hp : (sigma - 1) * perpSq e k * a e = 0 := by
      unfold extraNumerator perpSq at hf
      unfold perpSq
      nlinarith
    have he := (mul_eq_zero.mp hp).resolve_left (mul_ne_zero (sub_ne_zero.mpr hs1) hq)
    exact (mem_ordinarySpace e k a).mpr ⟨he, by simpa [he] using hd⟩
  · intro ha
    rw [ordinary_operator ha, sub_self, zero_smul]

theorem extra_kernel_iff {e : Fin D} {sigma : ℝ} {k a : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hk : k ≠ 0) :
    normalizedOperator e sigma (extraValue e sigma k) k a = 0 ↔
      a ∈ longitudinalSpace (extraVector e sigma k) := by
  have hf : sigma * extraValue e sigma k = extraNumerator e sigma k := by
    unfold extraValue
    field_simp
  constructor
  · intro h
    exact Submodule.mem_span_singleton.mpr
      ⟨_, (off_ordinary_reconstruction hs1 hq (ne_of_gt (extraValue_pos hs hk)) hf h).symm⟩
  · intro hl
    obtain ⟨c, rfl⟩ := Submodule.mem_span_singleton.mp hl
    rw [normalizedOperator_smul, extra_operator, hf, sub_self, zero_smul, smul_zero]

/-- Parallel wavevectors: one positive branch, with D-1 transverse amplitudes. -/
theorem parallel_kernel_iff {e : Fin D} {sigma z : ℝ} {k a : Vec D}
    (hs : sigma ≠ 0) (hz : z ≠ 0) (hk : k ≠ 0) (hq : perpSq e k = 0) :
    normalizedOperator e sigma z k a = 0 ↔
      a = 0 ∨ (z = normSq k ∧ a ∈ ordinarySpace e k) := by
  have hpar := (perpSq_zero_iff e k).mp hq
  have hp : k e ≠ 0 := by
    intro hp
    exact hk (by simpa [hp] using hpar)
  constructor
  · intro h
    have hd := weighted_transverse hz h
    have hdpar : dot k a = k e * a e := by
      calc
        dot k a = dot (k e • unit e) a := congrArg (fun v => dot v a) hpar
        _ = k e * a e := by rw [dot_smul_left, dot_unit_left]
    have he : a e = 0 := by
      rw [hdpar] at hd
      have hm : (sigma * k e) * a e = 0 := by nlinarith
      exact (mul_eq_zero.mp hm).resolve_left (mul_ne_zero hs hp)
    have ht := (mem_ordinarySpace e k a).mpr ⟨he, by rw [hdpar, he, mul_zero]⟩
    rw [ordinary_operator ht] at h
    rcases smul_eq_zero.mp h with hf | ha
    · exact Or.inr ⟨sub_eq_zero.mp hf, ht⟩
    · exact Or.inl ha
  · rintro (rfl | ⟨hf, ht⟩)
    · simp [normalizedOperator, dot_zero]
    · rw [ordinary_operator ht, hf, sub_self, zero_smul]

theorem zero_frequency_iff {e : Fin D} {sigma rho mu : ℝ} {k a : Vec D}
    (hmu : mu ≠ 0) (hk : k ≠ 0) :
    ModalStationary e sigma rho mu 0 k a ↔ a ∈ longitudinalSpace k := by
  rw [modal_stationary_iff]
  simp only [modalOperator, zero_pow (by decide : 2 ≠ 0), mul_zero, zero_mul, zero_smul, add_zero]
  rw [← S10Pilot.modal_stationary_iff]
  exact S10Pilot.zero_frequency_iff hmu hk

end
end S10Anisotropic
