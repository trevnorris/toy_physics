import Mathlib.MeasureTheory.Integral.CircleIntegral
import Mathlib.Tactic

/-! NP2: actual circle integrals of finite Laurent operands. -/
namespace S11NonlinearPole
noncomputable section
open Complex Real Metric

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℂ E]

/-- Positive orientation, with the usual 1/(2πi) normalization. -/
def moment (R : ℝ) (f : ℂ → E) : E :=
  (2 * Real.pi * Complex.I : ℂ)⁻¹ • (∮ z in C(0, R), f z)

theorem normalization_ne_zero : (2 * Real.pi * Complex.I : ℂ) ≠ 0 := by
  exact mul_ne_zero (mul_ne_zero (by norm_num) (by exact_mod_cast Real.pi_ne_zero)) I_ne_zero

theorem circle_zpow (R : ℝ) (hR : 0 < R) (n : ℤ) :
    CircleIntegrable (fun z : ℂ => z ^ n) 0 R := by
  simpa using (circleIntegrable_sub_zpow_iff (c := (0 : ℂ)) (w := 0) (R := R)
    (n := n)).mpr (Or.inr (Or.inr (by simp [abs_of_pos hR, hR.ne])))

theorem moment_zpow (R : ℝ) (hR : 0 < R) (n : ℤ) :
    moment R (fun z : ℂ => z ^ n) = if n = -1 then 1 else 0 := by
  by_cases hn : n = -1
  · subst n
    simp only [moment, zpow_neg_one, smul_eq_mul]
    have hi := circleIntegral.integral_sub_center_inv (0 : ℂ) hR.ne'
    simp only [sub_zero] at hi
    rw [hi, inv_mul_cancel₀ normalization_ne_zero]
    simp
  · simp only [moment, if_neg hn]
    have hi := circleIntegral.integral_sub_zpow_of_ne hn (0 : ℂ) 0 R
    simpa using congrArg (fun w : ℂ => (2 * Real.pi * Complex.I : ℂ)⁻¹ • w) hi

theorem moment_congr {R : ℝ} (hR : 0 ≤ R) {f g : ℂ → E}
    (h : Set.EqOn f g (sphere 0 R)) : moment R f = moment R g := by
  unfold moment
  rw [circleIntegral.integral_congr hR h]

theorem moment_add {R : ℝ} {f g : ℂ → E}
    (hf : CircleIntegrable f 0 R) (hg : CircleIntegrable g 0 R) :
    moment R (fun z => f z + g z) = moment R f + moment R g := by
  simp only [moment, circleIntegral.integral_add hf hg, smul_add]

theorem moment_smul (R : ℝ) (c : ℂ) (f : ℂ → E) :
    moment R (fun z => c • f z) = c • moment R f := by
  simp only [moment, circleIntegral.integral_smul, smul_comm c]

variable [CompleteSpace E]

theorem moment_zpow_smul (R : ℝ) (hR : 0 < R) (n : ℤ) (v : E) :
    moment R (fun z : ℂ => z ^ n • v) = if n = -1 then v else 0 := by
  simp only [moment, circleIntegral.integral_smul_const, ← smul_assoc]
  change (moment R (fun z : ℂ => z ^ n)) • v = _
  rw [moment_zpow R hR]
  split_ifs <;> simp

theorem moment_finite_laurent {ι : Type*} (s : Finset ι) (R : ℝ) (hR : 0 < R)
    (n : ι → ℤ) (v : ι → E) :
    moment R (fun z => ∑ i ∈ s, z ^ n i • v i) =
      ∑ i ∈ s, if n i = -1 then v i else 0 := by
  unfold moment
  rw [circleIntegral.integral_fun_sum (s := s) (f := fun i z => z ^ n i • v i) (fun i _ => (circle_zpow R hR (n i)).smul_continuousOn (continuousOn_const (c := v i)))]
  rw [Finset.smul_sum]
  exact Finset.sum_congr rfl (fun i _ => moment_zpow_smul R hR (n i) (v i))

/-- Terms with nonnegative exponent have zero contour moment, proved by integration. -/
theorem moment_polynomial {ι : Type*} (s : Finset ι) (R : ℝ) (hR : 0 < R)
    (n : ι → ℕ) (v : ι → E) :
    moment R (fun z => ∑ i ∈ s, z ^ n i • v i) = 0 := by
  have h := moment_finite_laurent s R hR (fun i => (n i : ℤ)) v
  simpa using h

theorem nonzero_on_circle {R : ℝ} (hR : 0 < R) {z : ℂ} (hz : z ∈ sphere 0 R) :
    z ≠ 0 := by
  rintro rfl
  simp [hR.ne] at hz

end
end S11NonlinearPole
