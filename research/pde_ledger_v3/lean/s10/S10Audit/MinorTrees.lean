import S10Audit.MatrixTrees

/-! Determinants are expression trees, with a different permitted unit for each
row. This covers both the modal matrix and the mixed-unit N3 stack. -/

set_option backward.isDefEq.respectTransparency false

namespace S10Audit
noncomputable section

def detTree {α : Type} : (q : ℕ) → (Fin q → Fin q → Expr α) → Expr α
  | 0, _ => .scalar 1
  | q + 1, M => .sum q fun j =>
      .mul (.mul (.scalar ((-1 : ℝ) ^ j.val)) (M 0 j))
        (detTree q fun i l => M i.succ (j.succAbove l))

theorem detTree_eval {α : Type} (v : α → ℝ) (q : ℕ) (M : Fin q → Fin q → Expr α) :
    Expr.eval v (detTree q M) = Matrix.det (fun i j => Expr.eval v (M i j)) := by
  induction q with
  | zero => exact Matrix.det_isEmpty.symm
  | succ q ih =>
    rw [Matrix.det_succ_row_zero]
    simp only [detTree, Expr.eval, ih, Matrix.submatrix]
    rfl

theorem detTree_hasDim {α : Type} (u : α → Dim) (q : ℕ)
    (M : Fin q → Fin q → Expr α) (row : Fin q → Dim)
    (h : ∀ i j, Expr.HasDim u (M i j) (row i)) :
    Expr.HasDim u (detTree q M) (∏ i, row i) := by
  induction q with
  | zero => simpa only [detTree, Fin.prod_univ_zero] using (Expr.HasDim.scalar (u := u) 1)
  | succ q ih =>
    rw [Fin.prod_univ_succ]
    apply Expr.HasDim.sum
    intro j
    simpa only [one_mul] using Expr.HasDim.mul
      (Expr.HasDim.mul (.scalar ((-1 : ℝ) ^ j.val)) (h 0 j))
      (ih (fun i l => M i.succ (j.succAbove l)) (fun i => row i.succ)
        (fun i l => h i.succ (j.succAbove l)))

/-- Row and column selections are explicit; repeated selections are allowed and
produce the usual zero determinant, still carrying the declared slot units. -/
def minorTree {α ι κ : Type} (M : ι → κ → Expr α) {q : ℕ}
    (rows : Fin q → ι) (cols : Fin q → κ) : Expr α :=
  detTree q fun i j => M (rows i) (cols j)

theorem minorTree_eval {α ι κ : Type} (v : α → ℝ) (M : ι → κ → Expr α) {q : ℕ}
    (rows : Fin q → ι) (cols : Fin q → κ) :
    Expr.eval v (minorTree M rows cols) =
      Matrix.det (Matrix.submatrix (fun i j => Expr.eval v (M i j)) rows cols) :=
  detTree_eval v q _

theorem minorTree_hasDim {α ι κ : Type} (u : α → Dim) (M : ι → κ → Expr α)
    (row : ι → Dim) (h : ∀ i j, Expr.HasDim u (M i j) (row i)) {q : ℕ}
    (rows : Fin q → ι) (cols : Fin q → κ) :
    Expr.HasDim u (minorTree M rows cols) (∏ i, row (rows i)) :=
  detTree_hasDim u q _ _ (fun i j => h (rows i) (cols j))

theorem matrix_minor_hasDim {n q : ℕ} (p : Package) (e : Fin (n + 1)) (field : Dim)
    {z : Tree (n + 1)}
    (hz : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) z (timeDim⁻¹ ^ (2 : ℕ)))
    (rows cols : Fin q → Fin (n + 1)) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field))
      (minorTree (matrixTree p e z) rows cols) (matrixDim (n + 1) field ^ q) := by
  simpa using minorTree_hasDim _ (matrixTree p e z) (fun _ => matrixDim (n + 1) field)
    (fun i j => matrixTree_hasDim p e i j field hz) rows cols

theorem stacked_minor_hasDim {n q : ℕ} (p : Package) (e : Fin (n + 1)) (field : Dim)
    {z : Tree (n + 1)}
    (hz : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) z (timeDim⁻¹ ^ (2 : ℕ)))
    (rows : Fin q → Fin (n + 2)) (cols : Fin q → Fin (n + 1)) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field))
      (minorTree (stackedTree (matrixTree p e z)) rows cols)
      (∏ i, stackedRowDim (n + 1) field (rows i)) :=
  minorTree_hasDim _ _ _ (stackedTree_hasDim p e field hz) rows cols

end
end S10Audit
