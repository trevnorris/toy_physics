import S11ScatteringSensitivity.Pipeline
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Algebra.BigOperators.Fin

/-! Exact witnesses and admissible domain cases; no measured physical inverse bound. -/
namespace S11ScatteringSensitivity
noncomputable section
open S11ScatteringFlux Matrix

def scalarA (z : ℂ) : Fin 1 → ℂ := fun _ => z
def scalarJ (z : ℂ) : Matrix (Fin 1) (Fin 1) ℂ := fun _ _ => z
def scalarC : Fin 1 → ℂ →L[ℂ] ℂ := fun _ => (3 : ℂ) • ContinuousLinearMap.id ℂ ℂ
def quarter : ℂ →L[ℂ] ℂ := (1 / 4 : ℂ) • ContinuousLinearMap.id ℂ ℂ

theorem residual_sign_witness : residual quarter 1 0 = (-1 : ℂ) := by
  norm_num [residual, quarter]

theorem inverse_factor_witness : ‖(4 : ℂ) - 0‖ = 4 * ‖residual quarter 1 0‖ := by
  norm_num [residual_sign_witness]

theorem small_residual_large_error :
    |(1 / 1024 : ℝ) * 1 - 0| = 1 / 1024 ∧ |(1 : ℝ) - 0| = 1 := by norm_num

theorem scaling_witness : (1 / 2 : ℝ) * 8 * (3 / 4) - (1 / 2) * 16 = -5 := by norm_num

theorem unscale_witness : (1 / 4 : ℝ) * (3 - 1) = 1 / 2 := by norm_num

theorem observed_error_witness :
    mass (amplitude scalarC (scalarA 1) 1 - amplitude scalarC (scalarA 1) 0) = 3 := by
  norm_num [mass, amplitude, scalarC, scalarA]

theorem affine_current_witness :
    flux (scalarJ 1) (amplitude scalarC (scalarA 1) 1) -
      flux (scalarJ 1) (amplitude scalarC (scalarA 1) 0) = 15 := by
  norm_num [flux, pair, amplitude, scalarC, scalarA, scalarJ, mulVec, dotProduct]

theorem full_error_witness :
    flux (scalarJ 1) (scalarA 3) - flux (scalarJ 1) (scalarA 2) = 5 := by
  norm_num [flux, pair, scalarA, scalarJ, mulVec, dotProduct]

theorem current_error_witness :
    flux (scalarJ 2) (scalarA 3) - flux (scalarJ 1) (scalarA 3) = 9 := by
  norm_num [flux, pair, scalarA, scalarJ, mulVec, dotProduct]

def fullJ : Matrix (Fin 2) (Fin 2) ℂ := !![1, Complex.I; -Complex.I, 1]
def complexA : Fin 2 → ℂ := ![1, -Complex.I]

theorem complex_current_witness : flux fullJ complexA = 4 := by
  norm_num [flux, pair, fullJ, complexA, mulVec, dotProduct, Fin.sum_univ_two]

theorem normalized_error_witness : |(3 : ℝ) / 2 - 2 / 4| = 1 := by norm_num
theorem denominator_only_witness : |(2 : ℝ) / 2 - 2 / 4| = 1 / 2 := by norm_num
theorem denominator_scale_witness : |(1 : ℝ) / (1 / 4) - 0 / (1 / 4)| = 4 := by norm_num
theorem signed_fraction_witness : |(3 : ℝ) / (-2) - 2 / (-4)| = 1 := by norm_num
theorem positive_margin_witness : |(2 : ℝ)| - 1 = 1 := by norm_num
theorem zero_margin_witness : |(0 : ℝ) - 1| = |(1 : ℝ)| := by norm_num
theorem critical_inverse_margin : ¬ ∃ x : ℂ, (1 + (-1 : ℂ)) * x = 1 := by norm_num

theorem identity_residual_positive (b x : ℂ) :
    ‖x - b‖ ≤ 1 * ‖residual (ContinuousLinearMap.id ℂ ℂ) b x‖ := by
  have h : ‖((ContinuousLinearEquiv.refl ℂ ℂ).symm : ℂ →L[ℂ] ℂ)‖ ≤ 1 :=
    ContinuousLinearMap.norm_id_le
  simpa only [ContinuousLinearEquiv.refl_symm, ContinuousLinearEquiv.refl_apply,
    ContinuousLinearEquiv.coe_refl] using residual_bound (ContinuousLinearEquiv.refl ℂ ℂ) h b x

theorem nonzero_margin_positive : (1 : ℝ) ≠ 0 :=
  denominator_nonzero 2 1 1 (by norm_num) (by norm_num)

theorem signed_fraction_positive : |(3 : ℝ) / (-2) - 2 / (-4)| ≤
    |(3 : ℝ) - 2| / 2 + |(2 : ℝ)| * |(-2 : ℝ) - (-4)| / 2 ^ 2 :=
  fraction_error_bound 2 3 (-4) (-2) 2 (by norm_num) (by norm_num) (by norm_num)

end
end S11ScatteringSensitivity
