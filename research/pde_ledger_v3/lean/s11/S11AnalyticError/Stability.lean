import Mathlib.Analysis.Normed.Operator.Banach
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Tactic

/-! T2: actual bounded inverse existence and quantitative perturbation estimates.
The outgoing/graph-space realization of the physical operator is a premise. -/
namespace S11AnalyticError
noncomputable section

variable {𝕜 X Y Z : Type*} [NontriviallyNormedField 𝕜]
  [NormedAddCommGroup X] [NormedSpace 𝕜 X]
  [NormedAddCommGroup Y] [NormedSpace 𝕜 Y]
  [NormedAddCommGroup Z] [NormedSpace 𝕜 Z]

def perturbation (A : X ≃L[𝕜] Y) (E : X →L[𝕜] Y) : X →L[𝕜] Y :=
  (A : X →L[𝕜] Y) + E

def perturbEquiv [CompleteSpace X] (A : X ≃L[𝕜] Y) (E : X →L[𝕜] Y)
    (h : ‖(A.symm : Y →L[𝕜] X).comp E‖ < 1) : X ≃L[𝕜] Y :=
  (ContinuousLinearEquiv.ofUnit (Units.oneSub (-(A.symm : Y →L[𝕜] X).comp E)
    (by simpa only [norm_neg] using h))).trans A

theorem perturbEquiv_eq [CompleteSpace X] (A : X ≃L[𝕜] Y) (E : X →L[𝕜] Y)
    (h : ‖(A.symm : Y →L[𝕜] X).comp E‖ < 1) :
    (perturbEquiv A E h : X →L[𝕜] Y) = perturbation A E := by
  ext x
  change A (x - -(A.symm (E x))) = A x + E x
  simp

theorem relative_error_small {A : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1) :
    ‖(A.symm : Y →L[𝕜] X).comp E‖ < 1 := by
  exact lt_of_le_of_lt ((ContinuousLinearMap.opNorm_comp_le _ _).trans
    (mul_le_mul hκ hε (norm_nonneg _) ((norm_nonneg _).trans hκ))) h

theorem coercive_bound {A : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (x : X) :
    (1 - κ * ε) * ‖x‖ ≤ κ * ‖perturbation A E x‖ := by
  have hκ0 : 0 ≤ κ := (norm_nonneg _).trans hκ
  have hx : x = A.symm (perturbation A E x) - A.symm (E x) := by
    simp [perturbation]
  have hfirst : ‖A.symm (perturbation A E x)‖ ≤ κ * ‖perturbation A E x‖ :=
    ((A.symm : Y →L[𝕜] X).le_opNorm _).trans
      (mul_le_mul_of_nonneg_right hκ (norm_nonneg _))
  have hsecond : ‖A.symm (E x)‖ ≤ κ * (ε * ‖x‖) := by
    calc
      _ ≤ κ * ‖E x‖ := ((A.symm : Y →L[𝕜] X).le_opNorm _).trans
        (mul_le_mul_of_nonneg_right hκ (norm_nonneg _))
      _ ≤ _ := mul_le_mul_of_nonneg_left
        ((E.le_opNorm x).trans (mul_le_mul_of_nonneg_right hε (norm_nonneg _))) hκ0
  have hz : ‖x‖ ≤ κ * ‖perturbation A E x‖ + κ * (ε * ‖x‖) := by
    calc
      ‖x‖ = ‖A.symm (perturbation A E x) - A.symm (E x)‖ := congrArg norm hx
      _ ≤ _ := (norm_sub_le _ _).trans (add_le_add hfirst hsecond)
  nlinarith

theorem inverse_norm_bound {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hB : (B : X →L[𝕜] Y) = perturbation A E)
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1) :
    ‖(B.symm : Y →L[𝕜] X)‖ ≤ κ / (1 - κ * ε) := by
  have hd : 0 < 1 - κ * ε := sub_pos.mpr h
  apply ContinuousLinearMap.opNorm_le_bound _ (div_nonneg ((norm_nonneg _).trans hκ) hd.le)
  intro y
  have hc := coercive_bound hκ hε (B.symm y)
  rw [← hB] at hc
  simp only [ContinuousLinearEquiv.coe_apply, ContinuousLinearEquiv.apply_symm_apply] at hc
  rw [div_mul_eq_mul_div, le_div_iff₀ hd]
  simpa only [ContinuousLinearEquiv.coe_apply, mul_comm] using hc

theorem resolvent_identity {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y}
    (hB : (B : X →L[𝕜] Y) = perturbation A E) :
    (B.symm : Y →L[𝕜] X) - (A.symm : Y →L[𝕜] X) =
      -((B.symm : Y →L[𝕜] X).comp (E.comp (A.symm : Y →L[𝕜] X))) := by
  ext y
  change B.symm y - A.symm y = -B.symm (E (A.symm y))
  apply B.injective
  rw [map_sub, map_neg, B.apply_symm_apply, B.apply_symm_apply]
  have heq : B (A.symm y) = y + E (A.symm y) := by
    have he := congrArg (fun F : X →L[𝕜] Y => F (A.symm y)) hB
    simpa [perturbation] using he
  rw [heq]
  abel

theorem inverse_difference_bound {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hB : (B : X →L[𝕜] Y) = perturbation A E)
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1) :
    ‖(B.symm : Y →L[𝕜] X) - (A.symm : Y →L[𝕜] X)‖ ≤
      κ ^ 2 * ε / (1 - κ * ε) := by
  rw [resolvent_identity hB, norm_neg]
  calc
    _ ≤ ‖(B.symm : Y →L[𝕜] X)‖ * ‖E.comp (A.symm : Y →L[𝕜] X)‖ :=
      ContinuousLinearMap.opNorm_comp_le _ _
    _ ≤ ‖(B.symm : Y →L[𝕜] X)‖ * (‖E‖ * ‖(A.symm : Y →L[𝕜] X)‖) :=
      mul_le_mul_of_nonneg_left (ContinuousLinearMap.opNorm_comp_le _ _) (norm_nonneg _)
    _ ≤ (κ / (1 - κ * ε)) * (ε * κ) := by
      exact mul_le_mul (inverse_norm_bound hB hκ hε h)
        (mul_le_mul hε hκ (norm_nonneg _) ((norm_nonneg _).trans hε))
        (mul_nonneg (norm_nonneg _) (norm_nonneg _))
        (div_nonneg ((norm_nonneg _).trans hκ) (sub_pos.mpr h).le)
    _ = _ := by ring

theorem solution_error_identity {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y}
    (hB : (B : X →L[𝕜] Y) = perturbation A E) (f df : Y) :
    B.symm (f + df) - A.symm f = B.symm (df - E (A.symm f)) := by
  apply B.injective
  rw [map_sub, B.apply_symm_apply, B.apply_symm_apply]
  have heq : B (A.symm f) = f + E (A.symm f) := by
    have he := congrArg (fun F : X →L[𝕜] Y => F (A.symm f)) hB
    simpa [perturbation] using he
  rw [heq]
  abel

theorem solution_error_bound {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hB : (B : X →L[𝕜] Y) = perturbation A E)
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1)
    (f df : Y) :
    ‖B.symm (f + df) - A.symm f‖ ≤
      (κ / (1 - κ * ε)) * (‖df‖ + ε * ‖A.symm f‖) := by
  rw [solution_error_identity hB]
  calc
    _ ≤ ‖(B.symm : Y →L[𝕜] X)‖ * ‖df - E (A.symm f)‖ :=
      (B.symm : Y →L[𝕜] X).le_opNorm _
    _ ≤ (κ / (1 - κ * ε)) * (‖df‖ + ε * ‖A.symm f‖) := by
      apply mul_le_mul (inverse_norm_bound hB hκ hε h)
      · exact (norm_sub_le _ _).trans (add_le_add (le_refl _)
          ((E.le_opNorm _).trans (mul_le_mul_of_nonneg_right hε (norm_nonneg _))))
      · exact norm_nonneg _
      · exact div_nonneg ((norm_nonneg _).trans hκ) (sub_pos.mpr h).le

theorem observation_error_bound {A B : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hB : (B : X →L[𝕜] Y) = perturbation A E)
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1)
    (C : X →L[𝕜] Z) (f df : Y) :
    ‖C (B.symm (f + df)) - C (A.symm f)‖ ≤
      ‖C‖ * ((κ / (1 - κ * ε)) * (‖df‖ + ε * ‖A.symm f‖)) := by
  rw [← map_sub]
  exact (C.le_opNorm _).trans
    (mul_le_mul_of_nonneg_left (solution_error_bound hB hκ hε h f df) (norm_nonneg _))

theorem exists_controlled_inverse [CompleteSpace X]
    {A : X ≃L[𝕜] Y} {E : X →L[𝕜] Y} {κ ε : ℝ}
    (hκ : ‖(A.symm : Y →L[𝕜] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1) :
    ∃ B : X ≃L[𝕜] Y, (B : X →L[𝕜] Y) = perturbation A E ∧
      ‖(B.symm : Y →L[𝕜] X)‖ ≤ κ / (1 - κ * ε) ∧
      ‖(B.symm : Y →L[𝕜] X) - (A.symm : Y →L[𝕜] X)‖ ≤
        κ ^ 2 * ε / (1 - κ * ε) := by
  let B := perturbEquiv A E (relative_error_small hκ hε h)
  have hB : (B : X →L[𝕜] Y) = perturbation A E := perturbEquiv_eq _ _ _
  exact ⟨B, hB, inverse_norm_bound hB hκ hε h, inverse_difference_bound hB hκ hε h⟩

end
end S11AnalyticError
