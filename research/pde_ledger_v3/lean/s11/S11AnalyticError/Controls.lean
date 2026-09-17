import S11AnalyticError.Abel
import S11AnalyticError.Stability
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

/-! T4: nonzero exact integral and inverse witnesses for the load-bearing bounds. -/
namespace S11AnalyticError
noncomputable section
open MeasureTheory Set

def atomTailError : ℝ :=
  ‖(∫ _ : ℝ, (1 : ℂ) ∂Measure.dirac 2) -
    ∫ _ : ℝ in Icc (-1) 1, (1 : ℂ) ∂Measure.dirac 2‖

theorem atomTailError_eq : atomTailError = 1 := by
  unfold atomTailError
  rw [setIntegral_dirac]
  norm_num

def complexTailError : ℝ :=
  ‖(∫ _ : ℝ, Complex.I ∂Measure.dirac 2) -
    ∫ _ : ℝ in Icc (-1) 1, Complex.I ∂Measure.dirac 2‖

theorem complexTailError_eq : complexTailError = 1 := by
  unfold complexTailError
  rw [setIntegral_dirac]
  norm_num

def atomAbelError : ℝ :=
  ‖(∫ _ : ℝ, (1 : ℂ) ∂Measure.dirac 1) -
    ∫ x : ℝ, (abelFactor 1 x : ℂ) ∂Measure.dirac 1‖

theorem atomAbelError_eq : atomAbelError = 1 - Real.exp (-1) := by
  simp only [atomAbelError, integral_dirac, abelFactor, neg_mul, one_mul, abs_one]
  rw [← Complex.ofReal_one, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
  exact abs_of_nonneg (sub_nonneg.mpr (Real.exp_le_one_iff.mpr (by norm_num)))

theorem atomAbelError_pos : 0 < atomAbelError := by
  rw [atomAbelError_eq]
  exact sub_pos.mpr (Real.exp_lt_one_iff.mpr (by norm_num))

theorem atomic_moment : firstMoment (Measure.dirac 1) (fun _ => (1 : ℂ)) = 1 := by
  simp [firstMoment]

theorem atomic_integrability : Integrable (fun _ : ℝ => (1 : ℂ)) (Measure.dirac 1) ∧
    Integrable (fun x : ℝ => |x| * ‖(1 : ℂ)‖) (Measure.dirac 1) := by
  constructor <;> exact integrable_dirac (by simp)

theorem negative_regulator_grows : 1 < abelFactor (-1) 1 := by
  simp [abelFactor, Real.one_lt_exp_iff]

theorem width_control : physicalWidth 1 10 = (1 / 10 : ℝ) := by
  norm_num [physicalWidth]

theorem conditioning_witness : |(1 + (-1 / 2 : ℝ))⁻¹ - 1| = 1 := by norm_num

theorem source_error_witness : |(1 + (1 / 4 : ℝ)) / (1 - 1 / 2) - 1| = 3 / 2 := by
  norm_num

theorem observation_witness : |3 * ((1 + (1 / 4 : ℝ)) / (1 - 1 / 2) - 1)| = 9 / 2 := by
  norm_num

theorem critical_margin_singular : ¬∃ r : ℝ, r * (1 + (-1)) = 1 := by norm_num

theorem strict_margin_admissible : (1 : ℝ) * (1 / 2) < 1 := by norm_num

theorem zero_abel_on_constant :
    (∫ x : ℝ, (abelFactor 0 x : ℂ) ∂Measure.dirac 1) = 1 := by simp [abel_zero]

theorem zero_amplitude_pairing (a : ℝ) :
    (∫ x : ℝ, (abelFactor a x : ℂ) * (0 : ℂ)) = 0 := by simp

end
end S11AnalyticError
