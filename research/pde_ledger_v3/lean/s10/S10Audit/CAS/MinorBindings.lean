import S10Audit.CAS.PYMinors
import S10Audit.CAS.WLMinors
import S10Audit.CAS.MinorReference

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYMinors

def staticRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n13, n19, n23, n26, n31, n34] : Fin 6 → Expr Symbol) i)
def staticRankIndex : Fin 6 → Fin 6 := ![0, 1, 2, 3, 4, 5]
theorem staticRankIndex_onto : Function.Surjective staticRankIndex := by decide

theorem staticRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, staticRank rho mu sigma z k i = staticRankMinors mu k (staticRankIndex i) := by
  intro i
  fin_cases i
  · exact root1_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell2 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell3 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell4 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell5 rho mu sigma z k

theorem staticRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, staticRank rho mu sigma z k i = 0) ↔ (∀ c, staticRankMinors mu k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := staticRankIndex_onto c
    rw [← staticRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [staticRank_reference rho mu sigma z k hs i]; exact h _

theorem staticRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, staticRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 0) 2 := by
  rw [staticRank_zero_iff rho mu sigma z k hs, CAS.staticRank_complete rho mu sigma k hr hs]

def staticTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n43, n47, n50, n53] : Fin 4 → Expr Symbol) i)
def staticTransverseIndex : Fin 4 → Fin 4 := ![0, 1, 2, 3]
theorem staticTransverseIndex_onto : Function.Surjective staticTransverseIndex := by decide

theorem staticTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, staticTransverse rho mu sigma z k i = staticTransverseMinors mu k (staticTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root1_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k

theorem staticTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, staticTransverse rho mu sigma z k i = 0) ↔ (∀ c, staticTransverseMinors mu k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := staticTransverseIndex_onto c
    rw [← staticTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [staticTransverse_reference rho mu sigma z k hs i]; exact h _

theorem staticTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, staticTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 0) 3 := by
  rw [staticTransverse_zero_iff rho mu sigma z k hs, CAS.staticTransverse_complete rho mu sigma k hr hs]

def ordinaryRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n58, n63, n68, n71] : Fin 4 → Expr Symbol) i)
def ordinaryRankIndex : Fin 4 → Fin 4 := ![0, 1, 2, 3]
theorem ordinaryRankIndex_onto : Function.Surjective ordinaryRankIndex := by decide

theorem ordinaryRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, ordinaryRank rho mu sigma z k i = ordinaryRankMinors mu sigma k (ordinaryRankIndex i) := by
  intro i
  fin_cases i
  · exact root2_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell2 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell3 rho mu sigma z k

theorem ordinaryRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, ordinaryRank rho mu sigma z k i = 0) ↔ (∀ c, ordinaryRankMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := ordinaryRankIndex_onto c
    rw [← ordinaryRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [ordinaryRank_reference rho mu sigma z k hs i]; exact h _

theorem ordinaryRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, ordinaryRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 1) 2 := by
  rw [ordinaryRank_zero_iff rho mu sigma z k hs, CAS.ordinaryRank_complete rho mu sigma k hr hs]

def ordinaryTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n58, n63, n68, n71, n76, n80] : Fin 6 → Expr Symbol) i)
def ordinaryTransverseIndex : Fin 6 → Fin 6 := ![0, 1, 2, 3, 4, 5]
theorem ordinaryTransverseIndex_onto : Function.Surjective ordinaryTransverseIndex := by decide

theorem ordinaryTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, ordinaryTransverse rho mu sigma z k i = ordinaryTransverseMinors mu sigma k (ordinaryTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root2_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell4 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell5 rho mu sigma z k

theorem ordinaryTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, ordinaryTransverse rho mu sigma z k i = 0) ↔ (∀ c, ordinaryTransverseMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := ordinaryTransverseIndex_onto c
    rw [← ordinaryTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [ordinaryTransverse_reference rho mu sigma z k hs i]; exact h _

theorem ordinaryTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, ordinaryTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 1) 2 := by
  rw [ordinaryTransverse_zero_iff rho mu sigma z k hs, CAS.ordinaryTransverse_complete rho mu sigma k hr hs]

def extraRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n85, n90, n95, n99, n102, n109] : Fin 6 → Expr Symbol) i)
def extraRankIndex : Fin 6 → Fin 6 := ![0, 1, 2, 3, 4, 5]
theorem extraRankIndex_onto : Function.Surjective extraRankIndex := by decide

theorem extraRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, extraRank rho mu sigma z k i = extraRankMinors mu sigma k (extraRankIndex i) := by
  intro i
  fin_cases i
  · exact root3_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell2 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell3 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell4 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell5 rho mu sigma z k _hs

theorem extraRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, extraRank rho mu sigma z k i = 0) ↔ (∀ c, extraRankMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := extraRankIndex_onto c
    rw [← extraRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [extraRank_reference rho mu sigma z k hs i]; exact h _

theorem extraRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, extraRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 2) 2 := by
  rw [extraRank_zero_iff rho mu sigma z k hs, CAS.extraRank_complete rho mu sigma k hr hs]

def extraTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n43, n115, n119, n122] : Fin 4 → Expr Symbol) i)
def extraTransverseIndex : Fin 4 → Fin 4 := ![0, 1, 2, 3]
theorem extraTransverseIndex_onto : Function.Surjective extraTransverseIndex := by decide

theorem extraTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, extraTransverse rho mu sigma z k i = extraTransverseMinors mu sigma k (extraTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root3_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root3_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k _hs
  · exact root3_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k _hs
  · exact root3_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k _hs

theorem extraTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, extraTransverse rho mu sigma z k i = 0) ↔ (∀ c, extraTransverseMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := extraTransverseIndex_onto c
    rw [← extraTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [extraTransverse_reference rho mu sigma z k hs i]; exact h _

theorem extraTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, extraTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 2) 3 := by
  rw [extraTransverse_zero_iff rho mu sigma z k hs, CAS.extraTransverse_complete rho mu sigma k hr hs]

end PYMinors

namespace WLMinors

def staticRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 9 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n13, n21, n25, n21, n28, n32, n25, n32, n35] : Fin 9 → Expr Symbol) i)
def staticRankIndex : Fin 9 → Fin 6 := ![0, 1, 2, 1, 3, 4, 2, 4, 5]
theorem staticRankIndex_onto : Function.Surjective staticRankIndex := by decide

theorem staticRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, staticRank rho mu sigma z k i = staticRankMinors mu k (staticRankIndex i) := by
  intro i
  fin_cases i
  · exact root1_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell2 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell3 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell4 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell5 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell6 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell7 rho mu sigma z k
  · exact root1_q8_rank_drop_minors_cell8 rho mu sigma z k

theorem staticRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, staticRank rho mu sigma z k i = 0) ↔ (∀ c, staticRankMinors mu k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := staticRankIndex_onto c
    rw [← staticRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [staticRank_reference rho mu sigma z k hs i]; exact h _

theorem staticRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, staticRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 0) 2 := by
  rw [staticRank_zero_iff rho mu sigma z k hs, CAS.staticRank_complete rho mu sigma k hr hs]

def staticTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n43, n47, n50, n53] : Fin 4 → Expr Symbol) i)
def staticTransverseIndex : Fin 4 → Fin 4 := ![0, 1, 2, 3]
theorem staticTransverseIndex_onto : Function.Surjective staticTransverseIndex := by decide

theorem staticTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, staticTransverse rho mu sigma z k i = staticTransverseMinors mu k (staticTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root1_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k
  · exact root1_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k

theorem staticTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, staticTransverse rho mu sigma z k i = 0) ↔ (∀ c, staticTransverseMinors mu k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := staticTransverseIndex_onto c
    rw [← staticTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [staticTransverse_reference rho mu sigma z k hs i]; exact h _

theorem staticTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, staticTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 0) 3 := by
  rw [staticTransverse_zero_iff rho mu sigma z k hs, CAS.staticTransverse_complete rho mu sigma k hr hs]

def ordinaryRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 9 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n57, n59, n64, n59, n66, n64, n64, n64, n64] : Fin 9 → Expr Symbol) i)
def ordinaryRankIndex : Fin 9 → Fin 4 := ![0, 1, 2, 1, 3, 2, 2, 2, 2]
theorem ordinaryRankIndex_onto : Function.Surjective ordinaryRankIndex := by decide

theorem ordinaryRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, ordinaryRank rho mu sigma z k i = ordinaryRankMinors mu sigma k (ordinaryRankIndex i) := by
  intro i
  fin_cases i
  · exact root2_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell2 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell3 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell4 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell5 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell6 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell7 rho mu sigma z k
  · exact root2_q8_rank_drop_minors_cell8 rho mu sigma z k

theorem ordinaryRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, ordinaryRank rho mu sigma z k i = 0) ↔ (∀ c, ordinaryRankMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := ordinaryRankIndex_onto c
    rw [← ordinaryRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [ordinaryRank_reference rho mu sigma z k hs i]; exact h _

theorem ordinaryRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, ordinaryRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 1) 2 := by
  rw [ordinaryRank_zero_iff rho mu sigma z k hs, CAS.ordinaryRank_complete rho mu sigma k hr hs]

def ordinaryTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 18 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n57, n59, n64, n59, n66, n64, n71, n75, n82, n64, n64, n64, n82, n82, n82, n82, n82, n82] : Fin 18 → Expr Symbol) i)
def ordinaryTransverseIndex : Fin 18 → Fin 6 := ![0, 1, 2, 1, 3, 2, 4, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
theorem ordinaryTransverseIndex_onto : Function.Surjective ordinaryTransverseIndex := by decide

theorem ordinaryTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, ordinaryTransverse rho mu sigma z k i = ordinaryTransverseMinors mu sigma k (ordinaryTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root2_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell4 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell5 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell6 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell7 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell8 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell9 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell10 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell11 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell12 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell13 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell14 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell15 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell16 rho mu sigma z k
  · exact root2_q8_transverse_rank_drop_minors_cell17 rho mu sigma z k

theorem ordinaryTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, ordinaryTransverse rho mu sigma z k i = 0) ↔ (∀ c, ordinaryTransverseMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := ordinaryTransverseIndex_onto c
    rw [← ordinaryTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [ordinaryTransverse_reference rho mu sigma z k hs i]; exact h _

theorem ordinaryTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, ordinaryTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 1) 2 := by
  rw [ordinaryTransverse_zero_iff rho mu sigma z k hs, CAS.ordinaryTransverse_complete rho mu sigma k hr hs]

def extraRank (rho mu sigma z : ℝ) (k : Vec 3) : Fin 9 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n86, n91, n97, n91, n101, n106, n97, n106, n112] : Fin 9 → Expr Symbol) i)
def extraRankIndex : Fin 9 → Fin 6 := ![0, 1, 2, 1, 3, 4, 2, 4, 5]
theorem extraRankIndex_onto : Function.Surjective extraRankIndex := by decide

theorem extraRank_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, extraRank rho mu sigma z k i = extraRankMinors mu sigma k (extraRankIndex i) := by
  intro i
  fin_cases i
  · exact root3_q8_rank_drop_minors_cell0 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell1 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell2 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell3 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell4 rho mu sigma z k
  · exact root3_q8_rank_drop_minors_cell5 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell6 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell7 rho mu sigma z k _hs
  · exact root3_q8_rank_drop_minors_cell8 rho mu sigma z k _hs

theorem extraRank_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, extraRank rho mu sigma z k i = 0) ↔ (∀ c, extraRankMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := extraRankIndex_onto c
    rw [← extraRank_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [extraRank_reference rho mu sigma z k hs i]; exact h _

theorem extraRank_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, extraRank rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootMatrix rho mu sigma k 2) 2 := by
  rw [extraRank_zero_iff rho mu sigma z k hs, CAS.extraRank_complete rho mu sigma k hr hs]

def extraTransverse (rho mu sigma z : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  fun i => Expr.eval (values rho mu sigma z k) ((![n43, n119, n123, n128] : Fin 4 → Expr Symbol) i)
def extraTransverseIndex : Fin 4 → Fin 4 := ![0, 1, 2, 3]
theorem extraTransverseIndex_onto : Function.Surjective extraTransverseIndex := by decide

theorem extraTransverse_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :
    ∀ i, extraTransverse rho mu sigma z k i = extraTransverseMinors mu sigma k (extraTransverseIndex i) := by
  intro i
  fin_cases i
  · exact root3_q8_transverse_rank_drop_minors_cell0 rho mu sigma z k
  · exact root3_q8_transverse_rank_drop_minors_cell1 rho mu sigma z k _hs
  · exact root3_q8_transverse_rank_drop_minors_cell2 rho mu sigma z k _hs
  · exact root3_q8_transverse_rank_drop_minors_cell3 rho mu sigma z k _hs

theorem extraTransverse_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, extraTransverse rho mu sigma z k i = 0) ↔ (∀ c, extraTransverseMinors mu sigma k c = 0) := by
  constructor
  · intro h c
    obtain ⟨i, rfl⟩ := extraTransverseIndex_onto c
    rw [← extraTransverse_reference rho mu sigma z k hs i]
    exact h i
  · intro h i; rw [extraTransverse_reference rho mu sigma z k hs i]; exact h _

theorem extraTransverse_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (∀ i, extraTransverse rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero (rootStack rho mu sigma k 2) 3 := by
  rw [extraTransverse_zero_iff rho mu sigma z k hs, CAS.extraTransverse_complete rho mu sigma k hr hs]

end WLMinors

theorem staticRank_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.staticRank rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.staticRank rho mu sigma z k i = 0) := by
  rw [PYMinors.staticRank_zero_iff rho mu sigma z k hs, WLMinors.staticRank_zero_iff rho mu sigma z k hs]

theorem staticTransverse_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.staticTransverse rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.staticTransverse rho mu sigma z k i = 0) := by
  rw [PYMinors.staticTransverse_zero_iff rho mu sigma z k hs, WLMinors.staticTransverse_zero_iff rho mu sigma z k hs]

theorem ordinaryRank_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.ordinaryRank rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.ordinaryRank rho mu sigma z k i = 0) := by
  rw [PYMinors.ordinaryRank_zero_iff rho mu sigma z k hs, WLMinors.ordinaryRank_zero_iff rho mu sigma z k hs]

theorem ordinaryTransverse_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.ordinaryTransverse rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.ordinaryTransverse rho mu sigma z k i = 0) := by
  rw [PYMinors.ordinaryTransverse_zero_iff rho mu sigma z k hs, WLMinors.ordinaryTransverse_zero_iff rho mu sigma z k hs]

theorem extraRank_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.extraRank rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.extraRank rho mu sigma z k i = 0) := by
  rw [PYMinors.extraRank_zero_iff rho mu sigma z k hs, WLMinors.extraRank_zero_iff rho mu sigma z k hs]

theorem extraTransverse_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :
    (∀ i, PYMinors.extraTransverse rho mu sigma z k i = 0) ↔
      (∀ i, WLMinors.extraTransverse rho mu sigma z k i = 0) := by
  rw [PYMinors.extraTransverse_zero_iff rho mu sigma z k hs, WLMinors.extraTransverse_zero_iff rho mu sigma z k hs]

end
end S10Audit.CAS
