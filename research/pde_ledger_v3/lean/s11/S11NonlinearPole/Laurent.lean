import S11NonlinearPole.Moments
import Mathlib.Analysis.Matrix.Normed

/-! NP2: ordered products are retained in actual finite Laurent contour operands. -/
namespace S11NonlinearPole
noncomputable section
open Complex Matrix
open scoped Matrix.Norms.Elementwise

variable {n : Type*} [Fintype n]

/-- The inverse principal part at a pole of order at most two. No inverse
existence is inferred from this definition. -/
def doublePrincipal (C₁ C₂ : Matrix n n ℂ) (z : ℂ) : Matrix n n ℂ :=
  z⁻¹ • C₁ + (z ^ 2)⁻¹ • C₂

/-- For an affine derivative A+zB, B is the second derivative of the pencil. -/
theorem double_log_expansion (C₁ C₂ A B : Matrix n n ℂ) {z : ℂ} (hz : z ≠ 0) :
    doublePrincipal C₁ C₂ z * (A + z • B) =
      ∑ k : Fin 4, z ^ (![-1,0,-2,-1] k : ℤ) •
        (![C₁*A,C₁*B,C₂*A,C₂*B] k) := by
  simp only [doublePrincipal, Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.mul_smul]
  simp [Fin.sum_univ_succ, zpow_neg]
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  field_simp
  ring

theorem double_log_moment (C₁ C₂ A B : Matrix n n ℂ) (R : ℝ) (hR : 0 < R) :
    moment R (fun z => doublePrincipal C₁ C₂ z * (A + z • B)) = C₁*A + C₂*B := by
  rw [moment_congr hR.le (fun z hz => double_log_expansion C₁ C₂ A B (nonzero_on_circle hR hz))]
  rw [moment_finite_laurent Finset.univ R hR]
  simp [Fin.sum_univ_succ]

variable {o u : Type*} [Fintype o] [Fintype u]

omit [Fintype o] [Fintype u] in
/-- Derivatives of both observation and forcing can contribute to the residue.
The six factors retain their matrix order. -/
theorem double_response_expansion (O₀ O₁ : Matrix o n ℂ) (C₁ C₂ : Matrix n n ℂ) (B₀ B₁ : Matrix n u ℂ)
    {z : ℂ} (hz : z ≠ 0) :
    (O₀ + z • O₁) * doublePrincipal C₁ C₂ z * (B₀ + z • B₁) =
      ∑ k : Fin 8, z ^ (![-1,0,-2,-1,0,1,-1,0] k : ℤ) •
        (![O₀*C₁*B₀,O₀*C₁*B₁,O₀*C₂*B₀,O₀*C₂*B₁,
           O₁*C₁*B₀,O₁*C₁*B₁,O₁*C₂*B₀,O₁*C₂*B₁] k) := by
  simp only [doublePrincipal, Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.mul_smul]
  simp [Fin.sum_univ_succ, zpow_neg]
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  field_simp
  ring

theorem double_response_moment (O₀ O₁ : Matrix o n ℂ) (C₁ C₂ : Matrix n n ℂ) (B₀ B₁ : Matrix n u ℂ)
    (R : ℝ) (hR : 0 < R) :
    moment R (fun z => (O₀ + z • O₁) * doublePrincipal C₁ C₂ z * (B₀ + z • B₁)) =
      O₀*C₁*B₀ + O₀*C₂*B₁ + O₁*C₂*B₀ := by
  rw [moment_congr hR.le (fun z hz =>
    double_response_expansion O₀ O₁ C₁ C₂ B₀ B₁ (nonzero_on_circle hR hz))]
  rw [moment_finite_laurent Finset.univ R hR]
  simp [Fin.sum_univ_succ, add_assoc]

end
end S11NonlinearPole
