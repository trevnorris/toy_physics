import S11AnalyticError.Tail

/-! T1/T3: coordinate Abel pairing, retaining constant and step contributions. -/
namespace S11AnalyticError
noncomputable section
open MeasureTheory Set

def abelFactor (a x : ℝ) : ℝ := Real.exp (-a * |x|)

theorem abelFactor_complex_continuous (a : ℝ) :
    Continuous (fun x => (abelFactor a x : ℂ)) := by
  unfold abelFactor
  fun_prop

theorem abelFactor_bounds {a : ℝ} (ha : 0 ≤ a) (x : ℝ) :
    0 ≤ abelFactor a x ∧ abelFactor a x ≤ 1 := by
  constructor
  · exact Real.exp_nonneg _
  · apply Real.exp_le_one_iff.mpr
    exact mul_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr ha) (abs_nonneg x)

theorem abel_error_bound {a : ℝ} (ha : 0 ≤ a) (x : ℝ) :
    |1 - abelFactor a x| ≤ a * |x| := by
  rw [abs_of_nonneg (sub_nonneg.mpr (abelFactor_bounds ha x).2)]
  have h := Real.add_one_le_exp (-a * |x|)
  dsimp [abelFactor]
  linarith

theorem abel_integrable {μ : Measure ℝ} {a : ℝ} (ha : 0 ≤ a)
    {f : ℝ → ℂ} (hf : Integrable f μ) :
    Integrable (fun x => (abelFactor a x : ℂ) * f x) μ := by
  apply hf.bdd_mul (abelFactor_complex_continuous a).aestronglyMeasurable
  exact ae_of_all _ (fun x => by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (abelFactor_bounds ha x).1]
    exact (abelFactor_bounds ha x).2)

theorem abel_pairing_bound {μ : Measure ℝ} {a B : ℝ} {b g : ℝ → ℂ}
    (ha : 0 ≤ a) (_hB0 : 0 ≤ B) (hg : Integrable g μ)
    (hm : Integrable (fun x => |x| * ‖g x‖) μ)
    (hb : AEStronglyMeasurable b μ) (hB : ∀ x, ‖b x‖ ≤ B) :
    ‖(∫ x, b x * g x ∂μ) - ∫ x, (abelFactor a x : ℂ) * (b x * g x) ∂μ‖
      ≤ B * a * firstMoment μ g := by
  have hf := bounded_product_integrable hg hb hB
  rw [← integral_sub hf (abel_integrable ha hf)]
  calc
    _ ≤ ∫ x, (B * a) * (|x| * ‖g x‖) ∂μ := by
      apply norm_integral_le_of_norm_le (hm.const_mul (B * a))
      apply ae_of_all
      intro x
      have hid : b x * g x - (abelFactor a x : ℂ) * (b x * g x) =
          ((1 - abelFactor a x : ℝ) : ℂ) * (b x * g x) := by push_cast; ring
      rw [hid, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs]
      calc
        _ ≤ (a * |x|) * (B * ‖g x‖) :=
          mul_le_mul (abel_error_bound ha x)
            (mul_le_mul_of_nonneg_right (hB x) (norm_nonneg _))
            (mul_nonneg (norm_nonneg _) (norm_nonneg _))
            (mul_nonneg ha (abs_nonneg x))
        _ = _ := by ring
    _ = _ := by rw [integral_const_mul]; rfl

theorem abel_truncation_bound {μ : Measure ℝ} {a B : ℝ} {b g : ℝ → ℂ}
    {s : Set ℝ} (hs : MeasurableSet s) (ha : 0 ≤ a) (hB0 : 0 ≤ B)
    (hg : Integrable g μ) (hm : Integrable (fun x => |x| * ‖g x‖) μ)
    (hb : AEStronglyMeasurable b μ) (hB : ∀ x, ‖b x‖ ≤ B) :
    ‖(∫ x, b x * g x ∂μ) - ∫ x in s, (abelFactor a x : ℂ) * (b x * g x) ∂μ‖
      ≤ B * (tailMass μ s g + a * firstMoment μ g) := by
  let c : ℝ → ℂ := fun x => (abelFactor a x : ℂ) * b x
  have hc : AEStronglyMeasurable c μ :=
    (abelFactor_complex_continuous a).aestronglyMeasurable.mul hb
  have hC : ∀ x, ‖c x‖ ≤ B := by
    intro x
    dsimp [c]
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg (abelFactor_bounds ha x).1]
    calc
      _ ≤ 1 * B := mul_le_mul (abelFactor_bounds ha x).2 (hB x) (norm_nonneg _) (by norm_num)
      _ = B := one_mul _
  have ht := truncation_error_bound hs hg hc hC
  simp only [c, mul_assoc] at ht
  have hab := abel_pairing_bound ha hB0 hg hm hb hB
  have htriangle := norm_sub_le_norm_sub_add_norm_sub
    (∫ x, b x * g x ∂μ)
    (∫ x, (abelFactor a x : ℂ) * (b x * g x) ∂μ)
    (∫ x in s, (abelFactor a x : ℂ) * (b x * g x) ∂μ)
  nlinarith

theorem abel_zero (x : ℝ) : abelFactor 0 x = 1 := by simp [abelFactor]

def physicalWidth (a ell : ℝ) : ℝ := a / ell

theorem physicalWidth_positive {a ell : ℝ} (ha : 0 < a) (hl : 0 < ell) :
    0 < physicalWidth a ell := div_pos ha hl

theorem approved_width (a : ℝ) : physicalWidth a 10 = a / 10 := rfl

end
end S11AnalyticError
