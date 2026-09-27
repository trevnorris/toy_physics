import S11ScatteringSensitivity.Residual
import S11ScatteringSensitivity.Current
import S11ScatteringSensitivity.Fraction

/-! Compose the actual finite-solve residual with the affine observation and
full current, without supplying a physical inverse or denominator certificate. -/
namespace S11ScatteringSensitivity
noncomputable section
open S11ScatteringFlux
variable {n X Y : Type*} [Fintype n]
  [NormedAddCommGroup X] [NormedSpace ℂ X]
  [NormedAddCommGroup Y] [NormedSpace ℂ Y]

def amplitudeBudget (κ : ℝ) (C : n → X →L[ℂ] ℂ) (A : X →L[ℂ] Y) (b : Y) (x : X) : ℝ :=
  observationBound C * (κ * ‖residual A b x‖)

def currentBudget (β δ : ℝ) (a : n → ℂ) : ℝ := β * (2 * mass a * δ + δ ^ 2)

theorem observed_residual_bound (A : X ≃L[ℂ] Y) {κ : ℝ}
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ)
    (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (b : Y) (x : X) :
    mass (amplitude C d x - amplitude C d (A.symm b)) ≤
      amplitudeBudget κ C (A : X →L[ℂ] Y) b x := by
  exact (amplitude_error_bound C d x (A.symm b)).trans
    (mul_le_mul_of_nonneg_left (residual_bound A hκ b x) (observation_bound_nonneg C))

theorem observed_perturbed_residual_bound {A B : X ≃L[ℂ] Y} {E : X →L[ℂ] Y} {κ ε : ℝ}
    (hB : (B : X →L[ℂ] Y) = S11AnalyticError.perturbation A E)
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ) (hε : ‖E‖ ≤ ε) (h : κ * ε < 1)
    (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (b db : Y) (x : X) :
    mass (amplitude C d x - amplitude C d (A.symm b)) ≤
      observationBound C * ((κ / (1 - κ * ε)) *
        (‖residual (B : X →L[ℂ] Y) (b + db) x‖ + ‖db‖ + ε * ‖A.symm b‖)) := by
  exact (amplitude_error_bound C d x (A.symm b)).trans
    (mul_le_mul_of_nonneg_left (perturbed_residual_bound hB hκ hε h b db x)
      (observation_bound_nonneg C))

theorem flux_distance_budget (J : Matrix n n ℂ) {β δ : ℝ} (hβ : 0 ≤ β)
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (a ah : n → ℂ) (he : mass (ah - a) ≤ δ) :
    |flux J ah - flux J a| ≤ currentBudget β δ a := by
  have ht := flux_change_budget J hβ hJ a (ah - a) he
  have hx : a + (ah - a) = ah := by abel
  rw [hx] at ht
  exact ht

theorem flux_residual_bound (A : X ≃L[ℂ] Y) {κ β : ℝ}
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ) (hβ : 0 ≤ β)
    (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (J : Matrix n n ℂ)
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (b : Y) (x : X) :
    |flux J (amplitude C d x) - flux J (amplitude C d (A.symm b))| ≤
      currentBudget β (amplitudeBudget κ C (A : X →L[ℂ] Y) b x)
        (amplitude C d (A.symm b)) :=
  flux_distance_budget J hβ hJ _ _ (observed_residual_bound A hκ C d b x)

theorem normalized_residual_bound (A : X ≃L[ℂ] Y) {κ β : ℝ}
    (hκ : ‖(A.symm : Y →L[ℂ] X)‖ ≤ κ) (hβ : 0 ≤ β)
    (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (J : Matrix n n ℂ)
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (b : Y) (x : X)
    (j jh margin η : ℝ) (hm : 0 < margin) (hj : margin ≤ |j|)
    (hjh : margin ≤ |jh|) (he : |jh - j| ≤ η) :
    |flux J (amplitude C d x) / jh - flux J (amplitude C d (A.symm b)) / j| ≤
      currentBudget β (amplitudeBudget κ C (A : X →L[ℂ] Y) b x)
        (amplitude C d (A.symm b)) / margin +
      |flux J (amplitude C d (A.symm b))| * η / margin ^ 2 :=
  fraction_error_budget _ _ j jh margin _ η hm hj hjh
    (flux_residual_bound A hκ hβ C d J hJ b x) he

end
end S11ScatteringSensitivity
