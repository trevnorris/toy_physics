import S11NonlinearPole.Laurent
import S11NonlinearPole.Scalar

/-! NP3: a defective affine state pencil still has a genuine projection.
Its specified physical transfer has a double pole and zero residue. -/
namespace S11NonlinearPole
noncomputable section
open Matrix
open scoped Matrix.Norms.Elementwise

abbrev Mat2 := Matrix (Fin 2) (Fin 2) ℂ

def jordanN : Mat2 := Matrix.of ![![0,1],![0,0]]

def jordanPencil (z : ℂ) : Mat2 := z • 1 - jordanN

def jordanResolvent (z : ℂ) : Mat2 := doublePrincipal 1 jordanN z

theorem jordan_nilpotent : jordanN * jordanN = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;> norm_num [jordanN, Matrix.mul_apply, Fin.sum_univ_two]

theorem jordan_nonzero : jordanN ≠ 0 := by
  intro h
  have := congrArg (fun M : Mat2 => M 0 1) h
  norm_num [jordanN] at this

theorem jordan_derivative (z : ℂ) : deriv jordanPencil z = 1 := by
  convert! (((hasDerivAt_id z).smul_const (1 : Mat2)).sub_const jordanN).deriv using 1
  simp

theorem jordan_inverse_left {z : ℂ} (hz : z ≠ 0) :
    jordanPencil z * jordanResolvent z = 1 := by
  simp only [jordanPencil, jordanResolvent, doublePrincipal, Matrix.sub_mul, Matrix.mul_add,
    Matrix.smul_mul, Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, jordan_nilpotent]
  ext i j
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply, Matrix.zero_apply, smul_eq_mul]
  field_simp
  ring

theorem jordan_inverse_right {z : ℂ} (hz : z ≠ 0) :
    jordanResolvent z * jordanPencil z = 1 := by
  simp only [jordanPencil, jordanResolvent, doublePrincipal, Matrix.mul_sub, Matrix.add_mul,
    Matrix.smul_mul, Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, jordan_nilpotent]
  ext i j
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply, Matrix.zero_apply, smul_eq_mul]
  field_simp
  ring

theorem jordan_determinant (z : ℂ) : (jordanPencil z).det = z ^ 2 := by
  simp [jordanPencil, jordanN, Matrix.det_fin_two, pow_two]

/-- The kernel consists exactly of the first coordinate axis, despite the
second-order determinant zero and nonzero nilpotent part. -/
theorem jordan_kernel (v : Fin 2 → ℂ) : jordanPencil 0 *ᵥ v = 0 ↔ v 1 = 0 := by
  simp [jordanPencil, jordanN, dotProduct, Fin.sum_univ_two, funext_iff,
    Fin.forall_fin_two]

theorem jordan_chain :
    jordanN *ᵥ ![0,1] = ![1,0] ∧ jordanN *ᵥ ![1,0] = 0 := by
  constructor <;> ext i <;> fin_cases i <;>
    norm_num [jordanN, Matrix.mulVec, dotProduct, Fin.sum_univ_two]

theorem jordan_residue (R : ℝ) (hR : 0 < R) : moment R jordanResolvent = 1 := by
  have heq : jordanResolvent = (fun z => ∑ k : Fin 2,
      z ^ (![-1,-2] k : ℤ) • (![1,jordanN] k : Mat2)) := by
    funext z
    simp [jordanResolvent, doublePrincipal, Fin.sum_univ_two, zpow_neg]
  rw [heq, moment_finite_laurent Finset.univ R hR]
  simp [Fin.sum_univ_two]

theorem jordan_higher_coefficient (R : ℝ) (hR : 0 < R) :
    moment R (fun z => z • jordanResolvent z) = jordanN := by
  have heq : ∀ z : ℂ, z ≠ 0 → z • jordanResolvent z =
      ∑ k : Fin 2, z ^ (![0,-1] k : ℤ) • (![1,jordanN] k : Mat2) := by
    intro z hz
    have h : z * (z ^ 2)⁻¹ = z⁻¹ := by field_simp
    simp [jordanResolvent, doublePrincipal, smul_add, smul_smul,
      Fin.sum_univ_two, mul_inv_cancel₀ hz, h]
  rw [moment_congr hR.le (fun z hz => heq z (nonzero_on_circle hR hz)),
    moment_finite_laurent Finset.univ R hR]
  simp [Fin.sum_univ_two]

def jordanLog (z : ℂ) : Mat2 := jordanResolvent z * deriv jordanPencil z

theorem jordan_log_moment (R : ℝ) (hR : 0 < R) : moment R jordanLog = 1 := by
  have heq : jordanLog = jordanResolvent := by
    funext z
    simp [jordanLog, jordan_derivative]
  rw [heq, jordan_residue R hR]

theorem jordan_projection_idempotent (R : ℝ) (hR : 0 < R) :
    moment R jordanLog * moment R jordanLog = moment R jordanLog := by
  rw [jordan_log_moment R hR, mul_one]

theorem jordan_count (R : ℝ) (hR : 0 < R) : (moment R jordanLog).trace = 2 := by
  rw [jordan_log_moment R hR]
  norm_num [Matrix.trace, Fin.sum_univ_two]

/-- Observation is the first state coordinate, forcing is the second coordinate.
This is the actual entry of the inverse just proved, not an arbitrary compression. -/
def jordanTransfer (z : ℂ) : ℂ := jordanResolvent z 0 1

theorem jordan_transfer_exact (z : ℂ) : jordanTransfer z = squareInverse z := by
  simp [jordanTransfer, jordanResolvent, doublePrincipal, jordanN, squareInverse, squarePencil]

theorem jordan_transfer_residue (R : ℝ) (hR : 0 < R) : moment R jordanTransfer = 0 := by
  have heq : jordanTransfer = squareInverse := funext jordan_transfer_exact
  rw [heq, square_residue R hR]

theorem jordan_transfer_nonzero : jordanTransfer 1 = 1 := by
  rw [jordan_transfer_exact, square_response_nonzero]

end
end S11NonlinearPole
