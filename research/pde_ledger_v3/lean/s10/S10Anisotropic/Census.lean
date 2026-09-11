import S10Anisotropic.Spectrum

/-! Nullities and dimensions of intersections with the Euclidean transverse
space. The perpendicular stratum is separate from the parallel stratum. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

def normalizedMap (e : Fin D) (sigma z : ℝ) (k : Vec D) : Vec D →ₗ[ℝ] Vec D where
  toFun := normalizedOperator e sigma z k
  map_add' a b := by ext i; simp [normalizedOperator, dot_add_right]; ring
  map_smul' c a := normalizedOperator_smul e sigma z c k a

def modeSpace (e : Fin D) (sigma z : ℝ) (k : Vec D) : Submodule ℝ (Vec D) :=
  (normalizedMap e sigma z k).ker

theorem mem_modeSpace (e : Fin D) (sigma z : ℝ) (k a : Vec D) :
    a ∈ modeSpace e sigma z k ↔ normalizedOperator e sigma z k a = 0 := Iff.rfl

theorem ordinary_modeSpace {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hk : k ≠ 0) :
    modeSpace e sigma (normSq k) k = ordinarySpace e k := by
  ext a
  exact ordinary_kernel_iff hs1 hq hk

theorem extra_modeSpace {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hk : k ≠ 0) :
    modeSpace e sigma (extraValue e sigma k) k = longitudinalSpace (extraVector e sigma k) := by
  ext a
  exact extra_kernel_iff hs hs1 hq hk

theorem parallel_modeSpace {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : sigma ≠ 0) (hk : k ≠ 0) (hq : perpSq e k = 0) :
    modeSpace e sigma (normSq k) k = ordinarySpace e k := by
  ext a
  rw [mem_modeSpace, parallel_kernel_iff hs (ne_of_gt (dot_self_pos hk)) hk hq]
  simp only [true_and]
  exact or_iff_right_of_imp (fun h => h ▸ (ordinarySpace e k).zero_mem)

theorem zero_modeSpace {e : Fin D} {sigma : ℝ} {k : Vec D} (hk : k ≠ 0) :
    modeSpace e sigma 0 k = longitudinalSpace k := by
  ext a
  have h := zero_frequency_iff (e := e) (sigma := sigma) (rho := 1) (mu := 1) (a := a)
    (by norm_num) hk
  rw [modal_stationary_normalized (by norm_num)] at h
  simpa [mem_modeSpace] using h

theorem ordinary_inf_transverse (e : Fin D) (k : Vec D) :
    ordinarySpace e k ⊓ transverseSpace k = ordinarySpace e k :=
  inf_eq_left.mpr inf_le_right

theorem span_inf_transverse_eq_bot {w k : Vec D} (hd : dot k w ≠ 0) :
    longitudinalSpace w ⊓ transverseSpace k = ⊥ := by
  ext a
  change (a ∈ longitudinalSpace w ∧ dot k a = 0) ↔ a = 0
  constructor
  · rintro ⟨hl, ht⟩
    obtain ⟨c, rfl⟩ := Submodule.mem_span_singleton.mp hl
    rw [dot_smul_right] at ht
    have hc := (mul_eq_zero.mp ht).resolve_right hd
    simp [hc]
  · rintro rfl
    simp

theorem extra_oblique_inf_transverse {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) (hp : k e ≠ 0) :
    longitudinalSpace (extraVector e sigma k) ⊓ transverseSpace k = ⊥ := by
  apply span_inf_transverse_eq_bot
  rw [dot_extraVector]
  exact mul_ne_zero (mul_ne_zero (sub_ne_zero.mpr hs1.symm) hp) hq

theorem extra_perpendicular_inf_transverse {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hp : k e = 0) :
    longitudinalSpace (extraVector e sigma k) ⊓ transverseSpace k =
      longitudinalSpace (extraVector e sigma k) := by
  apply inf_eq_left.mpr
  intro a ha
  obtain ⟨c, rfl⟩ := Submodule.mem_span_singleton.mp ha
  change dot k (c • extraVector e sigma k) = 0
  rw [dot_smul_right, dot_extraVector, hp]
  ring

/-- Generic oblique direction: (N2,N3) = (D-2,D-2) and (1,0). -/
theorem oblique_counts {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hk : k ≠ 0)
    (hq : perpSq e k ≠ 0) (hp : k e ≠ 0) :
    Module.finrank ℝ (modeSpace e sigma (normSq k) k) = D - 2 ∧
    Module.finrank ℝ (modeSpace e sigma (normSq k) k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = D - 2 ∧
    Module.finrank ℝ (modeSpace e sigma (extraValue e sigma k) k) = 1 ∧
    Module.finrank ℝ (modeSpace e sigma (extraValue e sigma k) k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = 0 := by
  rw [ordinary_modeSpace hs1 hq hk, extra_modeSpace hs hs1 hq hk,
    ordinary_inf_transverse, extra_oblique_inf_transverse hs1 hq hp,
    ordinarySpace_finrank hq, longitudinalSpace_finrank (extraVector_ne_zero hq)]
  simp

/-- Perpendicular direction: the split extra branch becomes exactly transverse.
The frequencies remain distinct when sigma != 1. -/
theorem perpendicular_counts {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hk : k ≠ 0) (hp : k e = 0) :
    Module.finrank ℝ (modeSpace e sigma (normSq k) k) = D - 2 ∧
    Module.finrank ℝ (modeSpace e sigma (normSq k) k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = D - 2 ∧
    Module.finrank ℝ (modeSpace e sigma (extraValue e sigma k) k) = 1 ∧
    Module.finrank ℝ (modeSpace e sigma (extraValue e sigma k) k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = 1 := by
  have hq : perpSq e k ≠ 0 := by simpa [perpSq, hp] using ne_of_gt (dot_self_pos hk)
  rw [ordinary_modeSpace hs1 hq hk, extra_modeSpace hs hs1 hq hk,
    ordinary_inf_transverse, extra_perpendicular_inf_transverse hp,
    ordinarySpace_finrank hq, longitudinalSpace_finrank (extraVector_ne_zero hq)]
  simp

/-- Parallel direction: the two formulas coincide and there is one D-1 dimensional root. -/
theorem parallel_counts {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : sigma ≠ 0) (hk : k ≠ 0) (hq : perpSq e k = 0) :
    Module.finrank ℝ (modeSpace e sigma (normSq k) k) = D - 1 ∧
    Module.finrank ℝ (modeSpace e sigma (normSq k) k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = D - 1 := by
  rw [parallel_modeSpace hs hk hq, ordinary_inf_transverse, ordinarySpace_parallel_finrank hq]
  exact ⟨rfl, rfl⟩

theorem zero_counts {e : Fin D} {sigma : ℝ} {k : Vec D} (hk : k ≠ 0) :
    Module.finrank ℝ (modeSpace e sigma 0 k) = 1 ∧
    Module.finrank ℝ (modeSpace e sigma 0 k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = 0 := by
  rw [zero_modeSpace hk, longitudinalSpace_finrank hk,
    span_inf_transverse_eq_bot (ne_of_gt (dot_self_pos hk))]
  simp

/-- Sum over the two distinct roots; valid even at D=2 where the ordinary
candidate frequency has zero nullity and is therefore not an actual mode. -/
theorem split_total_count {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hk : k ≠ 0) (hq : perpSq e k ≠ 0) :
    Module.finrank ℝ (modeSpace e sigma (normSq k) k) +
      Module.finrank ℝ (modeSpace e sigma (extraValue e sigma k) k) = D - 1 := by
  rw [ordinary_modeSpace hs1 hq hk, extra_modeSpace hs hs1 hq hk,
    ordinarySpace_finrank hq, longitudinalSpace_finrank (extraVector_ne_zero hq)]
  have hD := two_le_dimension hq
  omega

end
end S10Anisotropic
