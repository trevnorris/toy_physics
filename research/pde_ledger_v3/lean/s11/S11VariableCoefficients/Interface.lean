import S11VariableCoefficients.Common
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-! VC3: actual two-interval integration by parts with common test trace.
No multidimensional trace theorem or distribution product is postulated. -/
namespace S11VariableCoefficients
noncomputable section
open MeasureTheory Set

/-- Separately regular fluxes may have unequal interface values.
The same test function supplies the common trace at c. -/
theorem split_integration_by_parts {a c b : ℝ}
    {pm pp dm dp h dh : ℝ → ℝ}
    (hm : ∀ x ∈ uIcc a c, HasDerivAt pm (dm x) x)
    (hp : ∀ x ∈ uIcc c b, HasDerivAt pp (dp x) x)
    (hml : ∀ x ∈ uIcc a c, HasDerivAt h (dh x) x)
    (hpl : ∀ x ∈ uIcc c b, HasDerivAt h (dh x) x)
    (him : IntervalIntegrable dm volume a c) (hip : IntervalIntegrable dp volume c b)
    (hihm : IntervalIntegrable dh volume a c) (hihp : IntervalIntegrable dh volume c b) :
    (∫ x in a..c, pm x * dh x) + (∫ x in c..b, pp x * dh x) =
      pp b * h b - pm a * h a + (pm c - pp c) * h c -
        (∫ x in a..c, dm x * h x) - (∫ x in c..b, dp x * h x) := by
  rw [intervalIntegral.integral_mul_deriv_eq_deriv_mul hm hml him hihm,
    intervalIntegral.integral_mul_deriv_eq_deriv_mul hp hpl hip hihp]
  ring

theorem compact_endpoint_split {a c b : ℝ}
    {pm pp dm dp h dh : ℝ → ℝ}
    (hm : ∀ x ∈ uIcc a c, HasDerivAt pm (dm x) x)
    (hp : ∀ x ∈ uIcc c b, HasDerivAt pp (dp x) x)
    (hml : ∀ x ∈ uIcc a c, HasDerivAt h (dh x) x)
    (hpl : ∀ x ∈ uIcc c b, HasDerivAt h (dh x) x)
    (him : IntervalIntegrable dm volume a c) (hip : IntervalIntegrable dp volume c b)
    (hihm : IntervalIntegrable dh volume a c) (hihp : IntervalIntegrable dh volume c b)
    (ha : h a = 0) (hb : h b = 0) :
    (∫ x in a..c, pm x * dh x) + (∫ x in c..b, pp x * dh x) =
      (pm c - pp c) * h c -
        (∫ x in a..c, dm x * h x) - (∫ x in c..b, dp x * h x) := by
  rw [split_integration_by_parts hm hp hml hpl him hip hihm hihp, ha, hb]
  ring

end
end S11VariableCoefficients
