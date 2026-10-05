import S11AnalyticError.Stability

/-! S1: a residual bound for a supplied invertible finite realization.
Norms, scaling maps and inverse bounds are explicit application premises. -/
namespace S11ScatteringSensitivity
noncomputable section

variable {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℂ X]
  [NormedAddCommGroup Y] [NormedSpace ℂ Y]

def residual (A : X →L[ℂ] Y) (b : Y) (x : X) : Y := A x - b

theorem residual_identity (A : X ≃L[ℂ] Y) (b : Y) (x : X) :
    x - A.symm b = A.symm (residual (A : X →L[ℂ] Y) b x) := by
  simp [residual]

theorem residual_bound (A : X ≃L[ℂ] Y) {κ : ℝ}
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ) (b : Y) (x : X) :
    ‖x - A.symm b‖ ≤ κ * ‖residual (A : X →L[ℂ] Y) b x‖ := by
  rw [residual_identity]
  exact ((A.symm : Y →L[ℂ] X).le_opNorm _).trans
    (mul_le_mul_of_nonneg_right hκ (norm_nonneg _))

theorem exact_residual (A : X ≃L[ℂ] Y) (b : Y) :
    residual (A : X →L[ℂ] Y) b (A.symm b) = 0 := by
  simp [residual]

theorem residual_zero_iff (A : X ≃L[ℂ] Y) (b : Y) (x : X) :
    residual (A : X →L[ℂ] Y) b x = 0 ↔ x = A.symm b := by
  constructor
  · intro h
    have he := residual_identity A b x
    rw [h, map_zero] at he
    exact sub_eq_zero.mp he
  · rintro rfl
    exact exact_residual A b

/-- R is the row scaling, D the unknown scaling; the new unknown is D x. -/
def balanced (A : X ≃L[ℂ] Y) (R : Y ≃L[ℂ] Y) (D : X ≃L[ℂ] X) : X ≃L[ℂ] Y :=
  D.symm.trans (A.trans R)

theorem balanced_apply (A : X ≃L[ℂ] Y) (R : Y ≃L[ℂ] Y) (D : X ≃L[ℂ] X)
    (z : X) : balanced A R D z = R (A (D.symm z)) := rfl

theorem balanced_residual (A : X ≃L[ℂ] Y) (R : Y ≃L[ℂ] Y) (D : X ≃L[ℂ] X)
    (b : Y) (x : X) :
    residual (balanced A R D : X →L[ℂ] Y) (R b) (D x) =
      R (residual (A : X →L[ℂ] Y) b x) := by
  simp [residual, balanced_apply, map_sub]

theorem unscale_error (D : X ≃L[ℂ] X) (z w : X) :
    ‖D.symm z - D.symm w‖ ≤ ‖(D.symm : X →L[ℂ] X)‖ * ‖z - w‖ := by
  rw [← map_sub]
  exact (D.symm : X →L[ℂ] X).le_opNorm _

/-- Reuse the reviewed T2 stability bound and add the actual numerical residual. -/
theorem perturbed_residual_bound {A B : X ≃L[ℂ] Y} {E : X →L[ℂ] Y} {κ ε : ℝ}
    (hB : (B : X →L[ℂ] Y) = S11AnalyticError.perturbation A E)
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1)
    (b db : Y) (x : X) :
    ‖x - A.symm b‖ ≤ (κ / (1 - κ * ε)) *
      (‖residual (B : X →L[ℂ] Y) (b + db) x‖ + ‖db‖ + ε * ‖A.symm b‖) := by
  have hr := residual_bound B (S11AnalyticError.inverse_norm_bound hB hκ hε h) (b + db) x
  have hs := S11AnalyticError.solution_error_bound hB hκ hε h b db
  calc
    _ = ‖(x - B.symm (b + db)) + (B.symm (b + db) - A.symm b)‖ := by
      rw [sub_add_sub_cancel]
    _ ≤ ‖x - B.symm (b + db)‖ + ‖B.symm (b + db) - A.symm b‖ := norm_add_le _ _
    _ ≤ _ := by nlinarith

end
end S11ScatteringSensitivity
