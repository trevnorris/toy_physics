import S10Audit.CAS.MinorSupport

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

def staticRankMap : Fin 3 → Fin 3 → Fin 6 := ![![0, 1, 2], ![1, 3, 4], ![2, 4, 5]]
theorem staticRankMap_onto : ∀ c, ∃ i j, staticRankMap i j = c := by decide

theorem staticRank_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, staticRankMinors mu k (staticRankMap i j) =
      Matrix.det ((rootMatrix rho mu sigma k 0).submatrix (choose23 i) (choose23 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 0 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 1) - (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 1 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 2) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 2 = (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 2) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 1 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 2 1) - (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 3 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 2 2) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 4 = (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 2 2) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 2 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 2 = (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 1) - (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 4 = (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 2) - (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change staticRankMinors mu k 5 = (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 2) - (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 1)
      minor_equal

theorem staticRank_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, staticRankMinors mu k c = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 0) 2 := by
  rw [← selection_table_complete (rootMatrix rho mu sigma k 0) choose23 choose23
    choose23_complete choose23_complete choose23_strict choose23_strict]
  constructor
  · intro h i j
    rw [← staticRank_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := staticRankMap_onto c
    rw [staticRank_det rho mu sigma k hr hs i j]
    exact h i j

def staticTransverseMap : Fin 4 → Fin 1 → Fin 4 := ![![0], ![1], ![2], ![3]]
theorem staticTransverseMap_onto : ∀ c, ∃ i j, staticTransverseMap i j = c := by decide

theorem staticTransverse_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, staticTransverseMinors mu k (staticTransverseMap i j) =
      Matrix.det ((rootStack rho mu sigma k 0).submatrix (choose34 i) (choose33 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change staticTransverseMinors mu k 0 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 2) - (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 1) - (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 2) + (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 0) + (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 1) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change staticTransverseMinors mu k 1 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 1) * (k 2) - (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 1 2) * (k 1) - (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 0) * (k 2) + (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 1 2) * (k 0) + (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 0) * (k 1) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 1 1) * (k 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change staticTransverseMinors mu k 2 = (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 2 1) * (k 2) - (rootMatrix rho mu sigma k 0 0 0) * (rootMatrix rho mu sigma k 0 2 2) * (k 1) - (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 2 0) * (k 2) + (rootMatrix rho mu sigma k 0 0 1) * (rootMatrix rho mu sigma k 0 2 2) * (k 0) + (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 2 0) * (k 1) - (rootMatrix rho mu sigma k 0 0 2) * (rootMatrix rho mu sigma k 0 2 1) * (k 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change staticTransverseMinors mu k 3 = (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 1) * (k 2) - (rootMatrix rho mu sigma k 0 1 0) * (rootMatrix rho mu sigma k 0 2 2) * (k 1) - (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 0) * (k 2) + (rootMatrix rho mu sigma k 0 1 1) * (rootMatrix rho mu sigma k 0 2 2) * (k 0) + (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 0) * (k 1) - (rootMatrix rho mu sigma k 0 1 2) * (rootMatrix rho mu sigma k 0 2 1) * (k 0)
      minor_equal

theorem staticTransverse_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, staticTransverseMinors mu k c = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 0) 3 := by
  rw [← selection_table_complete (rootStack rho mu sigma k 0) choose34 choose33
    choose34_complete choose33_complete choose34_strict choose33_strict]
  constructor
  · intro h i j
    rw [← staticTransverse_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := staticTransverseMap_onto c
    rw [staticTransverse_det rho mu sigma k hr hs i j]
    exact h i j

def ordinaryRankMap : Fin 3 → Fin 3 → Fin 4 := ![![0, 1, 2], ![1, 3, 2], ![2, 2, 2]]
theorem ordinaryRankMap_onto : ∀ c, ∃ i j, ordinaryRankMap i j = c := by decide

theorem ordinaryRank_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, ordinaryRankMinors mu sigma k (ordinaryRankMap i j) =
      Matrix.det ((rootMatrix rho mu sigma k 1).submatrix (choose23 i) (choose23 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 0 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 1 1) - (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 1 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 1 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 1 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 1 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 1 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 2 1) - (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 3 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 2 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (rootMatrix rho mu sigma k 1 2 1) - (rootMatrix rho mu sigma k 1 1 1) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 1 2) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 1) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 1 2) * (rootMatrix rho mu sigma k 1 2 1)
      minor_equal

theorem ordinaryRank_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, ordinaryRankMinors mu sigma k c = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 1) 2 := by
  rw [← selection_table_complete (rootMatrix rho mu sigma k 1) choose23 choose23
    choose23_complete choose23_complete choose23_strict choose23_strict]
  constructor
  · intro h i j
    rw [← ordinaryRank_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := ordinaryRankMap_onto c
    rw [ordinaryRank_det rho mu sigma k hr hs i j]
    exact h i j

def ordinaryTransverseMap : Fin 6 → Fin 3 → Fin 6 := ![![0, 1, 2], ![1, 3, 2], ![4, 5, 2], ![2, 2, 2], ![2, 2, 2], ![2, 2, 2]]
theorem ordinaryTransverseMap_onto : ∀ c, ∃ i j, ordinaryTransverseMap i j = c := by decide

theorem ordinaryTransverse_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, ordinaryTransverseMinors mu sigma k (ordinaryTransverseMap i j) =
      Matrix.det ((rootStack rho mu sigma k 1).submatrix (choose24 i) (choose23 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 0 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 1 1) - (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 1 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 1 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 1 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 1 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 1 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 2 1) - (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 3 = (rootMatrix rho mu sigma k 1 0 0) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 0 1) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 0 2) * (rootMatrix rho mu sigma k 1 2 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 4 = (rootMatrix rho mu sigma k 1 0 0) * (k 1) - (rootMatrix rho mu sigma k 1 0 1) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 5 = (rootMatrix rho mu sigma k 1 0 0) * (k 2) - (rootMatrix rho mu sigma k 1 0 2) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 0 1) * (k 2) - (rootMatrix rho mu sigma k 1 0 2) * (k 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (rootMatrix rho mu sigma k 1 2 1) - (rootMatrix rho mu sigma k 1 1 1) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 1 2) * (rootMatrix rho mu sigma k 1 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 1) * (rootMatrix rho mu sigma k 1 2 2) - (rootMatrix rho mu sigma k 1 1 2) * (rootMatrix rho mu sigma k 1 2 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (k 1) - (rootMatrix rho mu sigma k 1 1 1) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 0) * (k 2) - (rootMatrix rho mu sigma k 1 1 2) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 1 1) * (k 2) - (rootMatrix rho mu sigma k 1 1 2) * (k 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 2 0) * (k 1) - (rootMatrix rho mu sigma k 1 2 1) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 2 0) * (k 2) - (rootMatrix rho mu sigma k 1 2 2) * (k 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change ordinaryTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 1 2 1) * (k 2) - (rootMatrix rho mu sigma k 1 2 2) * (k 1)
      minor_equal

theorem ordinaryTransverse_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, ordinaryTransverseMinors mu sigma k c = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 1) 2 := by
  rw [← selection_table_complete (rootStack rho mu sigma k 1) choose24 choose23
    choose24_complete choose23_complete choose24_strict choose23_strict]
  constructor
  · intro h i j
    rw [← ordinaryTransverse_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := ordinaryTransverseMap_onto c
    rw [ordinaryTransverse_det rho mu sigma k hr hs i j]
    exact h i j

def extraRankMap : Fin 3 → Fin 3 → Fin 6 := ![![0, 1, 2], ![1, 3, 4], ![2, 4, 5]]
theorem extraRankMap_onto : ∀ c, ∃ i j, extraRankMap i j = c := by decide

theorem extraRank_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, extraRankMinors mu sigma k (extraRankMap i j) =
      Matrix.det ((rootMatrix rho mu sigma k 2).submatrix (choose23 i) (choose23 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 0 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 1) - (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 1 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 2) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 2) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 1 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 2 1) - (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 3 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 2 2) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 4 = (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 2 2) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 2 1)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 2 = (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 1) - (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 4 = (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 2) - (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 0)
      minor_equal
    · rw [Matrix.det_fin_two]
      change extraRankMinors mu sigma k 5 = (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 2) - (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 1)
      minor_equal

theorem extraRank_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, extraRankMinors mu sigma k c = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 2) 2 := by
  rw [← selection_table_complete (rootMatrix rho mu sigma k 2) choose23 choose23
    choose23_complete choose23_complete choose23_strict choose23_strict]
  constructor
  · intro h i j
    rw [← extraRank_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := extraRankMap_onto c
    rw [extraRank_det rho mu sigma k hr hs i j]
    exact h i j

def extraTransverseMap : Fin 4 → Fin 1 → Fin 4 := ![![0], ![1], ![2], ![3]]
theorem extraTransverseMap_onto : ∀ c, ∃ i j, extraTransverseMap i j = c := by decide

theorem extraTransverse_det (rho mu sigma : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    ∀ i j, extraTransverseMinors mu sigma k (extraTransverseMap i j) =
      Matrix.det ((rootStack rho mu sigma k 2).submatrix (choose34 i) (choose33 j)) := by
  intro i j
  fin_cases i
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change extraTransverseMinors mu sigma k 0 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 2) - (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 1) - (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 2) + (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 0) + (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 1) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change extraTransverseMinors mu sigma k 1 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 1) * (k 2) - (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 1 2) * (k 1) - (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 0) * (k 2) + (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 1 2) * (k 0) + (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 0) * (k 1) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 1 1) * (k 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change extraTransverseMinors mu sigma k 2 = (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 2 1) * (k 2) - (rootMatrix rho mu sigma k 2 0 0) * (rootMatrix rho mu sigma k 2 2 2) * (k 1) - (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 2 0) * (k 2) + (rootMatrix rho mu sigma k 2 0 1) * (rootMatrix rho mu sigma k 2 2 2) * (k 0) + (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 2 0) * (k 1) - (rootMatrix rho mu sigma k 2 0 2) * (rootMatrix rho mu sigma k 2 2 1) * (k 0)
      minor_equal
  · fin_cases j
    · rw [Matrix.det_fin_three]
      change extraTransverseMinors mu sigma k 3 = (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 1) * (k 2) - (rootMatrix rho mu sigma k 2 1 0) * (rootMatrix rho mu sigma k 2 2 2) * (k 1) - (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 0) * (k 2) + (rootMatrix rho mu sigma k 2 1 1) * (rootMatrix rho mu sigma k 2 2 2) * (k 0) + (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 0) * (k 1) - (rootMatrix rho mu sigma k 2 1 2) * (rootMatrix rho mu sigma k 2 2 1) * (k 0)
      minor_equal

theorem extraTransverse_complete (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ c, extraTransverseMinors mu sigma k c = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 2) 3 := by
  rw [← selection_table_complete (rootStack rho mu sigma k 2) choose34 choose33
    choose34_complete choose33_complete choose34_strict choose33_strict]
  constructor
  · intro h i j
    rw [← extraTransverse_det rho mu sigma k hr hs i j]
    exact h _
  · intro h c
    obtain ⟨i, j, rfl⟩ := extraTransverseMap_onto c
    rw [extraTransverse_det rho mu sigma k hr hs i j]
    exact h i j

end
end S10Audit.CAS
