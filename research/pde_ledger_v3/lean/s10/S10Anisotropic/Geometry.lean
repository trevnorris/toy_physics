import S10Anisotropic.Action
import S10Pilot.Spectrum

/-! Geometry of the distinguished inertia axis and the wavevector. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

@[simp] theorem unit_self (e : Fin D) : unit e e = 1 := by simp [unit]

theorem unit_ne_zero (e : Fin D) : unit e ≠ 0 := by
  intro h
  have := congrFun h e
  simp at this

def perpPart (e : Fin D) (k : Vec D) : Vec D := k - k e • unit e
def perpSq (e : Fin D) (k : Vec D) : ℝ := normSq k - k e ^ 2

theorem dot_sub_left (a b c : Vec D) : dot (a - b) c = dot a c - dot b c := by
  simp [dot, sub_mul, Finset.sum_sub_distrib]

theorem dot_sub_right (a b c : Vec D) : dot a (b - c) = dot a b - dot a c := by
  simp [dot, mul_sub, Finset.sum_sub_distrib]

@[simp] theorem perpPart_axis (e : Fin D) (k : Vec D) : perpPart e k e = 0 := by
  simp [perpPart]

theorem dot_perpPart (e : Fin D) (k : Vec D) : dot k (perpPart e k) = perpSq e k := by
  simp [perpPart, perpSq, normSq, sub_eq_add_neg, ← neg_smul, dot_add_right, dot_smul_right, dot_unit_right]
  ring

theorem normSq_perpPart (e : Fin D) (k : Vec D) : normSq (perpPart e k) = perpSq e k := by
  simp only [normSq, perpPart, dot_sub_left, dot_sub_right,
    dot_smul_left, dot_smul_right, dot_unit_left, dot_unit_right,
    Pi.sub_apply, Pi.smul_apply, smul_eq_mul, unit_self, perpSq]
  ring

theorem perpSq_nonneg (e : Fin D) (k : Vec D) : 0 ≤ perpSq e k := by
  rw [← normSq_perpPart]
  by_cases h : perpPart e k = 0
  · simp [h, normSq]
  · exact le_of_lt (dot_self_pos h)

theorem perpSq_zero_iff (e : Fin D) (k : Vec D) :
    perpSq e k = 0 ↔ k = k e • unit e := by
  rw [← sub_eq_zero (a := k), ← perpPart]
  constructor
  · intro h
    by_contra hn
    have hp := dot_self_pos hn
    rw [normSq_perpPart, h] at hp
    exact (lt_irrefl 0) hp
  · intro h
    rw [← normSq_perpPart, h]
    simp [normSq]

def ordinarySpace (e : Fin D) (k : Vec D) : Submodule ℝ (Vec D) :=
  transverseSpace (unit e) ⊓ transverseSpace k

theorem mem_ordinarySpace (e : Fin D) (k a : Vec D) :
    a ∈ ordinarySpace e k ↔ a e = 0 ∧ dot k a = 0 := by
  change (dot (unit e) a = 0 ∧ dot k a = 0) ↔ _
  rw [dot_unit_left]

def constraintMap (e : Fin D) (k : Vec D) : Vec D →ₗ[ℝ] Vec 2 where
  toFun a := ![a e, dot k a]
  map_add' a b := by ext i; fin_cases i <;> simp [dot_add_right]
  map_smul' s a := by ext i; fin_cases i <;> simp [dot_smul_right]

theorem constraintMap_surjective {e : Fin D} {k : Vec D} (hq : perpSq e k ≠ 0) :
    Function.Surjective (constraintMap e k) := by
  intro b
  refine ⟨b 0 • unit e + ((b 1 - k e * b 0) / perpSq e k) • perpPart e k, ?_⟩
  ext i
  fin_cases i
  · simp [constraintMap]
  · change dot k (b 0 • unit e + ((b 1 - k e * b 0) / perpSq e k) • perpPart e k) = b 1
    rw [dot_add_right, dot_smul_right, dot_smul_right, dot_unit_right, dot_perpPart]
    rw [div_mul_cancel₀ _ hq]
    ring

theorem two_le_dimension {e : Fin D} {k : Vec D} (hq : perpSq e k ≠ 0) : 2 ≤ D := by
  have h := LinearMap.finrank_le_finrank_of_surjective (constraintMap_surjective hq)
  simpa using h

theorem ordinarySpace_finrank {e : Fin D} {k : Vec D} (hq : perpSq e k ≠ 0) :
    Module.finrank ℝ (ordinarySpace e k) = D - 2 := by
  have hs := constraintMap_surjective hq
  have heq : (constraintMap e k).ker = ordinarySpace e k := by
    ext a
    rw [mem_ordinarySpace]
    change (![a e, dot k a] : Vec 2) = 0 ↔ _
    simp [funext_iff, Fin.forall_fin_two]
  have h := (constraintMap e k).finrank_range_add_finrank_ker
  rw [LinearMap.range_eq_top.mpr hs, heq] at h
  simp at h
  omega

theorem ordinarySpace_parallel {e : Fin D} {k : Vec D} (hq : perpSq e k = 0) :
    ordinarySpace e k = transverseSpace (unit e) := by
  ext a
  rw [mem_ordinarySpace]
  change (a e = 0 ∧ dot k a = 0) ↔ dot (unit e) a = 0
  rw [dot_unit_left, (perpSq_zero_iff e k).mp hq, dot_smul_left, dot_unit_left]
  simp +contextual

theorem ordinarySpace_parallel_finrank {e : Fin D} {k : Vec D} (hq : perpSq e k = 0) :
    Module.finrank ℝ (ordinarySpace e k) = D - 1 := by
  rw [ordinarySpace_parallel hq]
  exact transverseSpace_finrank (unit_ne_zero e)

def extraNumerator (e : Fin D) (sigma : ℝ) (k : Vec D) : ℝ :=
  perpSq e k + sigma * k e ^ 2

def extraValue (e : Fin D) (sigma : ℝ) (k : Vec D) : ℝ := extraNumerator e sigma k / sigma

def extraVector (e : Fin D) (sigma : ℝ) (k : Vec D) : Vec D :=
  extraNumerator e sigma k • unit e - (sigma * k e) • k

theorem extraVector_axis (e : Fin D) (sigma : ℝ) (k : Vec D) :
    extraVector e sigma k e = perpSq e k := by
  simp [extraVector, extraNumerator]
  ring

theorem dot_extraVector (e : Fin D) (sigma : ℝ) (k : Vec D) :
    dot k (extraVector e sigma k) = (1 - sigma) * k e * perpSq e k := by
  simp only [extraVector, extraNumerator, perpSq, normSq,
    dot_sub_right, dot_smul_right, dot_unit_right]
  ring

theorem extraVector_ne_zero {e : Fin D} {sigma : ℝ} {k : Vec D} (hq : perpSq e k ≠ 0) :
    extraVector e sigma k ≠ 0 := by
  intro h
  have hi := congrFun h e
  rw [extraVector_axis] at hi
  exact hq hi

theorem extraValue_pos {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hk : k ≠ 0) : 0 < extraValue e sigma k := by
  have hq := perpSq_nonneg e k
  have hn := dot_self_pos hk
  have ha : 0 < extraNumerator e sigma k := by
    unfold extraNumerator
    by_cases hz : k e = 0
    · simpa [perpSq, hz] using hn
    · exact add_pos_of_nonneg_of_pos hq (mul_pos hs (sq_pos_of_ne_zero hz))
  exact div_pos ha hs

theorem extraValue_ne_ordinary {e : Fin D} {sigma : ℝ} {k : Vec D}
    (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) (hq : perpSq e k ≠ 0) :
    extraValue e sigma k ≠ normSq k := by
  intro h
  have hm := (div_eq_iff hs).mp h
  have hz : (sigma - 1) * perpSq e k = 0 := by
    unfold extraNumerator perpSq at hm
    unfold perpSq
    nlinarith
  exact (mul_ne_zero (sub_ne_zero.mpr hs1) hq) hz

end
end S10Anisotropic
