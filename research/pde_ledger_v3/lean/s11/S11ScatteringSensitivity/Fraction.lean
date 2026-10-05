import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-! S3: real normalized observables with explicit absolute denominator margins. -/
namespace S11ScatteringSensitivity

theorem denominator_margin (j jh δ : ℝ) (herr : |jh - j| ≤ δ) :
    |j| - δ ≤ |jh| := by
  have ht : |j| ≤ |jh| + |jh - j| := by
    calc
      |j| = |jh + (j - jh)| := by congr 1; ring
      _ ≤ |jh| + |j - jh| := abs_add_le _ _
      _ = _ := by rw [abs_sub_comm j jh]
  linarith

theorem denominator_nonzero (j jh δ : ℝ) (herr : |jh - j| ≤ δ) (hδ : δ < |j|) :
    jh ≠ 0 := by
  have hm := denominator_margin j jh δ herr
  intro hzero
  rw [hzero, abs_zero] at hm
  linarith

theorem fraction_difference (n nh j jh : ℝ) (hj : j ≠ 0) (hjh : jh ≠ 0) :
    nh / jh - n / j = (nh - n) / jh + n * (j - jh) / (jh * j) := by
  field_simp
  ring

theorem fraction_error_bound (n nh j jh d : ℝ) (hd : 0 < d)
    (hj : d ≤ |j|) (hjh : d ≤ |jh|) :
    |nh / jh - n / j| ≤ |nh - n| / d + |n| * |jh - j| / d ^ 2 := by
  have hj0 : j ≠ 0 := abs_pos.mp (hd.trans_le hj)
  have hjh0 : jh ≠ 0 := abs_pos.mp (hd.trans_le hjh)
  have hp : d ^ 2 ≤ |jh| * |j| := by
    simpa only [pow_two] using mul_le_mul hjh hj hd.le (abs_nonneg jh)
  rw [fraction_difference n nh j jh hj0 hjh0]
  calc
    _ ≤ |(nh - n) / jh| + |n * (j - jh) / (jh * j)| := abs_add_le _ _
    _ = |nh - n| / |jh| + (|n| * |jh - j|) / (|jh| * |j|) := by
      rw [abs_div, abs_div, abs_mul, abs_mul, abs_sub_comm j jh]
    _ ≤ _ := add_le_add
      (div_le_div_of_nonneg_left (abs_nonneg _) hd hjh)
      (div_le_div_of_nonneg_left (mul_nonneg (abs_nonneg _) (abs_nonneg _)) (sq_pos_of_pos hd) hp)

theorem fraction_error_budget (n nh j jh d α δ : ℝ) (hd : 0 < d)
    (hj : d ≤ |j|) (hjh : d ≤ |jh|) (hn : |nh - n| ≤ α) (he : |jh - j| ≤ δ) :
    |nh / jh - n / j| ≤ α / d + |n| * δ / d ^ 2 := by
  exact (fraction_error_bound n nh j jh d hd hj hjh).trans (add_le_add
    (div_le_div_of_nonneg_right hn hd.le)
    (div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left he (abs_nonneg _)) (sq_nonneg d)))

theorem denominator_cases (j : ℝ) : j = 0 ∨ 0 < |j| := by
  by_cases h : j = 0
  · exact Or.inl h
  · exact Or.inr (abs_pos.mpr h)

end S11ScatteringSensitivity
