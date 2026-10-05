import S10Audit.MinorTrees

/-! Unit changes preserve the entire kernel and rank, including for a stack
whose final wavevector row scales differently from the dynamical rows. -/

set_option backward.isDefEq.respectTransparency false

namespace S10Audit
noncomputable section

def rowScale {ι κ : Type} (c : ι → ℝ) (M : Matrix ι κ ℝ) : Matrix ι κ ℝ :=
  fun i j => c i * M i j

theorem rowScale_mulVec {ι κ : Type} [Fintype κ] (c : ι → ℝ) (M : Matrix ι κ ℝ)
    (a : κ → ℝ) (i : ι) : (rowScale c M).mulVec a i = c i * M.mulVec a i := by
  simp only [rowScale, Matrix.mulVec, dotProduct, mul_assoc, Finset.mul_sum]

theorem rowScale_kernel_iff {ι κ : Type} [Fintype κ] (c : ι → ℝ)
    (hc : ∀ i, c i ≠ 0) (M : Matrix ι κ ℝ) (a : κ → ℝ) :
    (rowScale c M).mulVec a = 0 ↔ M.mulVec a = 0 := by
  simp only [funext_iff, rowScale_mulVec, Pi.zero_apply, mul_eq_zero, hc, false_or]

theorem rowScale_rank {ι κ : Type} [Fintype ι] [Fintype κ] [DecidableEq ι]
    (c : ι → ℝ) (hc : ∀ i, c i ≠ 0) (M : Matrix ι κ ℝ) :
    (rowScale c M).rank = M.rank := by
  have he : rowScale c M = Matrix.diagonal c * M := by
    ext i j
    exact (Matrix.diagonal_mul c M i j).symm
  rw [he]
  apply Matrix.rank_mul_eq_right_of_det_ne_zero
  simpa only [Matrix.det_diagonal, Finset.prod_ne_zero_iff] using fun i _ => hc i

theorem matrixTree_unit_covariance {α ι κ : Type} (u : α → Dim) (v : α → ℝ)
    (M : ι → κ → Expr α) (row : ι → Dim)
    (h : ∀ i j, Expr.HasDim u (M i j) (row i)) (χ : Dim →* ℝˣ) :
    (fun i j => Expr.eval (fun a => (χ (u a) : ℝ) * v a) (M i j)) =
      rowScale (fun i => (χ (row i) : ℝ)) (fun i j => Expr.eval v (M i j)) := by
  funext i j
  exact (h i j).rescale χ v

theorem matrixTree_unit_kernel_iff {α ι κ : Type} [Fintype κ]
    (u : α → Dim) (v : α → ℝ) (M : ι → κ → Expr α) (row : ι → Dim)
    (h : ∀ i j, Expr.HasDim u (M i j) (row i)) (χ : Dim →* ℝˣ) (a : κ → ℝ) :
    Matrix.mulVec (fun i j => Expr.eval (fun b => (χ (u b) : ℝ) * v b) (M i j)) a = 0 ↔
      Matrix.mulVec (fun i j => Expr.eval v (M i j)) a = 0 := by
  rw [matrixTree_unit_covariance u v M row h χ]
  exact rowScale_kernel_iff _ (fun i => Units.ne_zero _) _ a

theorem matrixTree_unit_rank {α ι κ : Type} [Fintype ι] [Fintype κ] [DecidableEq ι]
    (u : α → Dim) (v : α → ℝ) (M : ι → κ → Expr α) (row : ι → Dim)
    (h : ∀ i j, Expr.HasDim u (M i j) (row i)) (χ : Dim →* ℝˣ) :
    Matrix.rank (fun i j => Expr.eval (fun b => (χ (u b) : ℝ) * v b) (M i j)) =
      Matrix.rank (fun i j => Expr.eval v (M i j)) := by
  rw [matrixTree_unit_covariance u v M row h χ]
  exact rowScale_rank _ (fun i => Units.ne_zero _) _

theorem minorTree_unit_zero_iff {α ι κ : Type} (u : α → Dim) (v : α → ℝ)
    (M : ι → κ → Expr α) (row : ι → Dim)
    (h : ∀ i j, Expr.HasDim u (M i j) (row i)) (χ : Dim →* ℝˣ) {q : ℕ}
    (rows : Fin q → ι) (cols : Fin q → κ) :
    Expr.eval (fun a => (χ (u a) : ℝ) * v a) (minorTree M rows cols) = 0 ↔
      Expr.eval v (minorTree M rows cols) = 0 := by
  rw [(minorTree_hasDim u M row h rows cols).rescale χ v]
  exact mul_eq_zero.trans (or_iff_right (Units.ne_zero _))

end
end S10Audit
