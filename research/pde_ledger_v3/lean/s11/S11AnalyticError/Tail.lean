import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.Complex.Norm
import Mathlib.Tactic

/-! T1: actual integral tails, with explicit integrability and multiplier bounds. -/
namespace S11AnalyticError
noncomputable section
open MeasureTheory Set

def tailMass (μ : Measure ℝ) (s : Set ℝ) (g : ℝ → ℂ) : ℝ :=
  ∫ x in sᶜ, ‖g x‖ ∂μ

def firstMoment (μ : Measure ℝ) (g : ℝ → ℂ) : ℝ :=
  ∫ x, |x| * ‖g x‖ ∂μ

theorem tailMass_nonneg (μ : Measure ℝ) (s : Set ℝ) (g : ℝ → ℂ) :
    0 ≤ tailMass μ s g := integral_nonneg (fun _ => norm_nonneg _)

theorem firstMoment_nonneg (μ : Measure ℝ) (g : ℝ → ℂ) :
    0 ≤ firstMoment μ g := integral_nonneg (fun x => mul_nonneg (abs_nonneg x) (norm_nonneg _))

theorem bounded_product_integrable {μ : Measure ℝ} {b g : ℝ → ℂ} {B : ℝ}
    (hg : Integrable g μ) (hb : AEStronglyMeasurable b μ)
    (hB : ∀ x, ‖b x‖ ≤ B) : Integrable (fun x => b x * g x) μ :=
  hg.bdd_mul hb (ae_of_all _ hB)

theorem omitted_integral_bound {μ : Measure ℝ} {b g : ℝ → ℂ} {B : ℝ}
    (s : Set ℝ) (hg : Integrable g μ) (hB : ∀ x, ‖b x‖ ≤ B) :
    ‖∫ x in sᶜ, b x * g x ∂μ‖ ≤ B * tailMass μ s g := by
  calc
    _ ≤ ∫ x in sᶜ, B * ‖g x‖ ∂μ :=
      norm_integral_le_of_norm_le (hg.norm.restrict.const_mul B)
        (ae_of_all _ (fun x => by rw [norm_mul]; exact mul_le_mul_of_nonneg_right (hB x) (norm_nonneg _)))
    _ = _ := by rw [integral_const_mul]; rfl

theorem truncation_error_bound {μ : Measure ℝ} {b g : ℝ → ℂ} {B : ℝ}
    {s : Set ℝ} (hs : MeasurableSet s) (hg : Integrable g μ)
    (hb : AEStronglyMeasurable b μ) (hB : ∀ x, ‖b x‖ ≤ B) :
    ‖(∫ x, b x * g x ∂μ) - ∫ x in s, b x * g x ∂μ‖ ≤ B * tailMass μ s g := by
  rw [← setIntegral_compl hs (bounded_product_integrable hg hb hB)]
  exact omitted_integral_bound s hg hB

theorem measurable_cutoff (R : ℝ) : MeasurableSet {x : ℝ | |x| ≤ R} :=
  isClosed_le continuous_abs continuous_const |>.measurableSet

theorem real_phase_norm (s x : ℝ) : ‖Complex.exp (-(s * x : ℝ) * Complex.I)‖ = 1 := by
  simp [Complex.norm_exp]

end
end S11AnalyticError
