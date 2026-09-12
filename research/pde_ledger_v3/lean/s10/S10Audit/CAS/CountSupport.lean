import S10Audit.CAS.RerunBindings
import Mathlib.LinearAlgebra.Matrix.Rank

/-! Counts refer to the actual imported matrix maps and basis index sets.
Residual subtraction takes place in ℤ, so a negative mismatch cannot be
silently truncated by natural-number subtraction. -/

namespace S10Audit.CAS
open S10Pilot
noncomputable section
set_option backward.isDefEq.respectTransparency false

def matrixNullity {m : ℕ} (M : Matrix (Fin m) (Fin 3) ℝ) : ℕ :=
  Module.finrank ℝ (LinearMap.ker M.mulVecLin)

def basisCount {n : ℕ} (_b : Fin n → Vec 3) : ℕ := Fintype.card (Fin n)

def nullityDifference (M : Matrix (Fin 3) (Fin 3) ℝ)
    (S : Matrix (Fin 4) (Fin 3) ℝ) : ℤ :=
  (matrixNullity M : ℤ) - (matrixNullity S : ℤ)

def basisCountResidual {n : ℕ} (b : Fin n → Vec 3) (M : Matrix (Fin 3) (Fin 3) ℝ) : ℤ :=
  (basisCount b : ℤ) - (matrixNullity M : ℤ)

def appendConstraint (M : Matrix (Fin 3) (Fin 3) ℝ) (p : Vec 3) :
    Matrix (Fin 4) (Fin 3) ℝ := ![M 0, M 1, M 2, p]

theorem appendConstraint_kernel (M : Matrix (Fin 3) (Fin 3) ℝ) (p : Vec 3) :
    LinearMap.ker (appendConstraint M p).mulVecLin =
      LinearMap.ker M.mulVecLin ⊓ transverseSpace p := by
  ext a
  change (appendConstraint M p).mulVec a = 0 ↔ M.mulVec a = 0 ∧ dot p a = 0
  constructor
  · intro h
    refine ⟨?_, ?_⟩
    · ext i
      fin_cases i
      · exact congrFun h (0 : Fin 4)
      · exact congrFun h (1 : Fin 4)
      · exact congrFun h (2 : Fin 4)
    · exact congrFun h (3 : Fin 4)
  · rintro ⟨hm, hp⟩
    ext i
    fin_cases i
    · exact congrFun hm 0
    · exact congrFun hm 1
    · exact congrFun hm 2
    · exact hp

theorem matrix_rank_nullity {m : ℕ} (M : Matrix (Fin m) (Fin 3) ℝ) :
    M.rank + matrixNullity M = 3 := by
  simpa [Matrix.rank, matrixNullity] using M.mulVecLin.finrank_range_add_finrank_ker

theorem matrix_rank_of_nullity {m n : ℕ} (M : Matrix (Fin m) (Fin 3) ℝ)
    (h : matrixNullity M = n) : M.rank = 3 - n := by
  have hn := matrix_rank_nullity M
  omega

theorem matrix_counts_of_kernel {m n : ℕ} (M : Matrix (Fin m) (Fin 3) ℝ)
    (N : Matrix (Fin n) (Fin 3) ℝ) (h : LinearMap.ker M.mulVecLin = LinearMap.ker N.mulVecLin) :
    matrixNullity M = matrixNullity N ∧ M.rank = N.rank := by
  have hn : matrixNullity M = matrixNullity N := by rw [matrixNullity, matrixNullity, h]
  refine ⟨hn, ?_⟩
  have hM := matrix_rank_nullity M
  have hN := matrix_rank_nullity N
  omega

theorem physical_modal_kernel (rho mu sigma z c : ℝ) (p : Vec 3) (hc : c ≠ 0) :
    LinearMap.ker (((1/2 : ℝ) • referenceMatrix rho mu sigma (c ^ 2 * z) (c • p)).mulVecLin) =
    LinearMap.ker (((1/2 : ℝ) • referenceMatrix rho mu sigma z p).mulVecLin) := by
  ext a
  exact coordinate_kernel_scale rho mu sigma z c p a hc

theorem physical_constraint_kernel (rho mu sigma z c : ℝ) (p : Vec 3) (hc : c ≠ 0) :
    LinearMap.ker (appendConstraint
      ((1/2 : ℝ) • referenceMatrix rho mu sigma (c ^ 2 * z) (c • p)) (c • p)).mulVecLin =
    LinearMap.ker (appendConstraint
      ((1/2 : ℝ) • referenceMatrix rho mu sigma z p) p).mulVecLin := by
  rw [appendConstraint_kernel, appendConstraint_kernel]
  ext a
  change (_ = 0 ∧ dot (c • p) a = 0) ↔ (_ = 0 ∧ dot p a = 0)
  simp only [Matrix.mulVecLin_apply]
  rw [coordinate_kernel_scale rho mu sigma z c p a hc, coordinate_dot_scale,
    mul_eq_zero, or_iff_right hc]

open Lean Elab Tactic in
/-- Integral count expressions need only exact arithmetic after evaluation. -/
elab "count_eval" h:Lean.Parser.Tactic.rwRule : tactic => do
  evalTactic (← `(tactic| rw [$h]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| norm_num))

end
end S10Audit.CAS
