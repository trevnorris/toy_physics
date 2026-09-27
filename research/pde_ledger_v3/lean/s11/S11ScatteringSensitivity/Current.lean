import S11ScatteringFlux.Current
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Normed.Operator.Basic

/-! S2: explicit finite l1 amplitude mass and the full supplied current form.
Entrywise current bounds are not spectral norms or physical conservation. -/
namespace S11ScatteringSensitivity
noncomputable section
open S11ScatteringFlux
open scoped BigOperators
variable {n : Type*} [Fintype n]

def mass (a : n → ℂ) : ℝ := ∑ i, ‖a i‖

theorem mass_nonneg (a : n → ℂ) : 0 ≤ mass a := Finset.sum_nonneg fun _ _ => norm_nonneg _

theorem mass_zero : mass (0 : n → ℂ) = 0 := by simp [mass]

theorem mass_add_le (a e : n → ℂ) : mass (a + e) ≤ mass a + mass e := by
  simpa only [mass, Pi.add_apply, Finset.sum_add_distrib] using
    Finset.sum_le_sum (fun i (_ : i ∈ (Finset.univ : Finset n)) => norm_add_le (a i) (e i))

theorem pair_bound (J : Matrix n n ℂ) {β : ℝ}
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (a b : n → ℂ) :
    ‖pair J a b‖ ≤ β * mass a * mass b := by
  rw [pair_coordinates]
  calc
    _ ≤ ∑ j, ‖∑ i, star (a i) * J i j * b j‖ := norm_sum_le _ _
    _ ≤ ∑ j, ∑ i, ‖star (a i) * J i j * b j‖ :=
      Finset.sum_le_sum fun j _ => norm_sum_le _ _
    _ ≤ ∑ j, ∑ i, β * ‖a i‖ * ‖b j‖ := by
      apply Finset.sum_le_sum
      intro j _
      apply Finset.sum_le_sum
      intro i _
      simp only [norm_mul, Complex.star_def, Complex.norm_conj]
      have he := mul_le_mul_of_nonneg_right
        (mul_le_mul_of_nonneg_left (hJ i j) (norm_nonneg (a i))) (norm_nonneg (b j))
      simpa only [mul_comm, mul_left_comm, mul_assoc] using he
    _ = β * mass a * mass b := by
      simp only [mass, ← Finset.sum_mul, ← Finset.mul_sum]

theorem flux_bound (J : Matrix n n ℂ) {β : ℝ}
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (a : n → ℂ) :
    |flux J a| ≤ β * mass a ^ 2 := by
  exact (Complex.abs_re_le_norm _).trans (by simpa [pow_two, mul_assoc] using pair_bound J hJ a a)

theorem flux_change_bound (J : Matrix n n ℂ) {β : ℝ}
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (a e : n → ℂ) :
    |flux J (a + e) - flux J a| ≤ β * (2 * mass a * mass e + mass e ^ 2) := by
  have heq : flux J (a + e) - flux J a = flux J e + (pair J a e + pair J e a).re := by
    rw [flux_add]
    ring
  rw [heq]
  calc
    _ ≤ |flux J e| + |(pair J a e + pair J e a).re| := abs_add_le _ _
    _ ≤ β * mass e ^ 2 + ‖pair J a e + pair J e a‖ :=
      add_le_add (flux_bound J hJ e) (Complex.abs_re_le_norm _)
    _ ≤ β * mass e ^ 2 + (‖pair J a e‖ + ‖pair J e a‖) :=
      add_le_add (le_refl _) (norm_add_le _ _)
    _ ≤ β * mass e ^ 2 + (β * mass a * mass e + β * mass e * mass a) :=
      add_le_add (le_refl _) (add_le_add (pair_bound J hJ a e) (pair_bound J hJ e a))
    _ = _ := by ring

theorem flux_change_budget (J : Matrix n n ℂ) {β δ : ℝ} (hβ : 0 ≤ β)
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (a e : n → ℂ) (he : mass e ≤ δ) :
    |flux J (a + e) - flux J a| ≤ β * (2 * mass a * δ + δ ^ 2) := by
  have hm := mass_nonneg a
  have hn := mass_nonneg e
  have hp : mass e ^ 2 ≤ δ ^ 2 := pow_le_pow_left₀ hn he 2
  have hl : 2 * mass a * mass e ≤ 2 * mass a * δ :=
    mul_le_mul_of_nonneg_left he (mul_nonneg (le_of_lt zero_lt_two) hm)
  exact (flux_change_bound J hJ a e).trans
    (mul_le_mul_of_nonneg_left (add_le_add hl hp) hβ)

theorem current_change_bound (J H : Matrix n n ℂ) {γ : ℝ}
    (hH : ∀ i j, ‖H i j‖ ≤ γ) (a : n → ℂ) :
    |flux (J + H) a - flux J a| ≤ γ * mass a ^ 2 := by
  have heq : flux (J + H) a - flux J a = flux H a := by
    simp only [flux, pair_add_matrix, Complex.add_re, add_sub_cancel_left]
  rw [heq]
  exact flux_bound H hH a

theorem flux_and_current_change (J H : Matrix n n ℂ) {β γ : ℝ}
    (hJ : ∀ i j, ‖J i j‖ ≤ β) (hH : ∀ i j, ‖H i j‖ ≤ γ) (a e : n → ℂ) :
    |flux (J + H) (a + e) - flux J a| ≤
      β * (2 * mass a * mass e + mass e ^ 2) + γ * mass (a + e) ^ 2 := by
  calc
    _ = |(flux J (a + e) - flux J a) + (flux (J + H) (a + e) - flux J (a + e))| := by
      congr 1
      ring
    _ ≤ _ := (abs_add_le _ _).trans
      (add_le_add (flux_change_bound J hJ a e) (current_change_bound J H hH (a + e)))

variable {X : Type*} [NormedAddCommGroup X] [NormedSpace ℂ X]

def amplitude (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (x : X) : n → ℂ := fun i => C i x + d i
def observationBound (C : n → X →L[ℂ] ℂ) : ℝ := ∑ i, ‖C i‖

omit [Fintype n] in
theorem amplitude_difference (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (x y : X) :
    amplitude C d x - amplitude C d y = fun i => C i (x - y) := by
  funext i
  simp [amplitude, map_sub]

theorem observation_bound_nonneg (C : n → X →L[ℂ] ℂ) : 0 ≤ observationBound C :=
  Finset.sum_nonneg fun _ _ => norm_nonneg _

theorem amplitude_error_bound (C : n → X →L[ℂ] ℂ) (d : n → ℂ) (x y : X) :
    mass (amplitude C d x - amplitude C d y) ≤ observationBound C * ‖x - y‖ := by
  rw [amplitude_difference]
  exact (Finset.sum_le_sum fun i _ => (C i).le_opNorm (x - y)).trans_eq (by
    simp only [observationBound, Finset.sum_mul])

end
end S11ScatteringSensitivity
