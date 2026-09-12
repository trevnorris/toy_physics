import S10Audit.CAS.Support

/-! Explicit, complete selection tables for the D3 rank-drop minors.
The small canonical families below are checked against every ordered minor. -/

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section

def choose23 : Fin 3 → Fin 2 → Fin 3 := ![![0, 1], ![0, 2], ![1, 2]]
def choose24 : Fin 6 → Fin 2 → Fin 4 :=
  ![![0, 1], ![0, 2], ![0, 3], ![1, 2], ![1, 3], ![2, 3]]
def choose33 : Fin 1 → Fin 3 → Fin 3 := ![![0, 1, 2]]
def choose34 : Fin 4 → Fin 3 → Fin 4 :=
  ![![0, 1, 2], ![0, 1, 3], ![0, 2, 3], ![1, 2, 3]]

theorem choose23_complete : ∀ f : Fin 2 → Fin 3, StrictMono f → ∃ i, choose23 i = f := by decide
theorem choose24_complete : ∀ f : Fin 2 → Fin 4, StrictMono f → ∃ i, choose24 i = f := by decide
theorem choose33_complete : ∀ f : Fin 3 → Fin 3, StrictMono f → ∃ i, choose33 i = f := by decide
theorem choose34_complete : ∀ f : Fin 3 → Fin 4, StrictMono f → ∃ i, choose34 i = f := by decide
theorem choose23_strict : ∀ i, StrictMono (choose23 i) := by decide
theorem choose24_strict : ∀ i, StrictMono (choose24 i) := by decide
theorem choose33_strict : ∀ i, StrictMono (choose33 i) := by decide
theorem choose34_strict : ∀ i, StrictMono (choose34 i) := by decide

def AllOrderedMinorsZero {m n : ℕ} (M : Matrix (Fin m) (Fin n) ℝ) (q : ℕ) : Prop :=
  ∀ rows : Fin q → Fin m, ∀ cols : Fin q → Fin n,
    StrictMono rows → StrictMono cols → Matrix.det (M.submatrix rows cols) = 0

theorem selection_table_complete {m n q a b : ℕ} (M : Matrix (Fin m) (Fin n) ℝ)
    (rows : Fin a → Fin q → Fin m) (cols : Fin b → Fin q → Fin n)
    (hr : ∀ f, StrictMono f → ∃ i, rows i = f)
    (hc : ∀ f, StrictMono f → ∃ i, cols i = f)
    (hrs : ∀ i, StrictMono (rows i)) (hcs : ∀ j, StrictMono (cols j)) :
    (∀ i j, Matrix.det (M.submatrix (rows i) (cols j)) = 0) ↔ AllOrderedMinorsZero M q := by
  constructor
  · intro h r c hsr hsc
    obtain ⟨i, rfl⟩ := hr r hsr
    obtain ⟨j, rfl⟩ := hc c hsc
    exact h i j
  · intro h i j
    exact h (rows i) (cols j) (hrs i) (hcs j)

def rootMatrix (rho mu sigma : ℝ) (k : Vec 3) (r : Fin 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  (1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k r) k

def rootStack (rho mu sigma : ℝ) (k : Vec 3) (r : Fin 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  ![rootMatrix rho mu sigma k r 0, rootMatrix rho mu sigma k r 1,
    rootMatrix rho mu sigma k r 2, k]

def staticRankMinors (mu : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  ![k 2 ^ 2 * mu ^ 2 * normSq k / 4,
    -(k 1 * k 2 * mu ^ 2 * normSq k) / 4,
    k 0 * k 2 * mu ^ 2 * normSq k / 4,
    k 1 ^ 2 * mu ^ 2 * normSq k / 4,
    -(k 0 * k 1 * mu ^ 2 * normSq k) / 4,
    k 0 ^ 2 * mu ^ 2 * normSq k / 4]

def staticTransverseMinors (mu : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  ![0, k 2 * mu ^ 2 * normSq k ^ 2 / 4,
    -(k 1 * mu ^ 2 * normSq k ^ 2) / 4,
    k 0 * mu ^ 2 * normSq k ^ 2 / 4]

def ordinaryRankMinors (mu sigma : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  ![k 1 ^ 2 * mu ^ 2 * (sigma - 1) * normSq k / 4,
    k 1 * k 2 * mu ^ 2 * (sigma - 1) * normSq k / 4, 0,
    k 2 ^ 2 * mu ^ 2 * (sigma - 1) * normSq k / 4]

def ordinaryTransverseMinors (mu sigma : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  ![ordinaryRankMinors mu sigma k 0, ordinaryRankMinors mu sigma k 1, 0,
    ordinaryRankMinors mu sigma k 3,
    k 1 * mu * (sigma - 1) * normSq k / 2,
    k 2 * mu * (sigma - 1) * normSq k / 2]

def extraRankMinors (mu sigma : ℝ) (k : Vec 3) : Fin 6 → ℝ :=
  ![-(k 0 ^ 2 * k 2 ^ 2 * mu ^ 2 * (sigma - 1)) / 4,
    k 0 ^ 2 * k 1 * k 2 * mu ^ 2 * (sigma - 1) / 4,
    k 0 * k 2 * mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) * (sigma - 1) / (4 * sigma),
    -(k 0 ^ 2 * k 1 ^ 2 * mu ^ 2 * (sigma - 1)) / 4,
    -(k 0 * k 1 * mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) * (sigma - 1)) / (4 * sigma),
    -(mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) ^ 2 * (sigma - 1)) / (4 * sigma ^ 2)]

def extraTransverseMinors (mu sigma : ℝ) (k : Vec 3) : Fin 4 → ℝ :=
  ![0, -(k 0 ^ 2 * k 2 * mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) * (sigma - 1) ^ 2) / (4 * sigma),
    k 0 ^ 2 * k 1 * mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) * (sigma - 1) ^ 2 / (4 * sigma),
    k 0 * mu ^ 2 * (k 1 ^ 2 + k 2 ^ 2) ^ 2 * (sigma - 1) ^ 2 / (4 * sigma ^ 2)]

open Lean Elab Tactic in
elab "minor_equal" : tactic => do
  evalTactic (← `(tactic| simp only [staticRankMinors, staticTransverseMinors,
    ordinaryRankMinors, ordinaryTransverseMinors, extraRankMinors, extraTransverseMinors,
    rootMatrix, rootStack, choose23, choose24, choose33, choose34, Matrix.submatrix,
    Matrix.det_fin_two, Matrix.det_fin_three, Matrix.cons_val,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.cons_val_four, Matrix.cons_val_succ',
    Matrix.head_cons, Matrix.tail_cons]))
  evalTactic (← `(tactic| cas_equal))

open Lean Elab Tactic in
elab "minor_eval" h:Lean.Parser.Tactic.rwRule : tactic => do
  evalTactic (← `(tactic| rw [$h]))
  if !(← getGoals).isEmpty then evalTactic (← `(tactic| minor_equal))

end
end S10Audit.CAS
