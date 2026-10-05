import S10Anisotropic.Certificate

/-! Q5 scaling for the extra branch, including its direction-dependent coefficient. -/
namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

theorem extraConeValue_scaling (e : Fin D) (sigma rho mu lambda : ℝ) (k : Vec D) :
    extraConeValue e sigma rho mu (lambda • k) = lambda ^ 2 * extraConeValue e sigma rho mu k := by
  simp only [extraConeValue, extraValue, extraNumerator, perpSq, normSq_smul,
    Pi.smul_apply, smul_eq_mul]
  ring

theorem extraConeValue_scaling_ratio {e : Fin D} {sigma rho mu : ℝ} {k : Vec D}
    (hs : 0 < sigma) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) (lambda : ℝ) :
    extraConeValue e sigma rho mu (lambda • k) / extraConeValue e sigma rho mu k = lambda ^ 2 := by
  rw [extraConeValue_scaling, mul_div_cancel_right₀]
  exact ne_of_gt (mul_pos (div_pos hmu hrho) (extraValue_pos hs hk))

end
end S10Anisotropic
