import S11D4Odd.Calculus
import S11D4Odd.Variation

/-! D4B.2–D4B.3: an explicit current and vanishing bulk variation, retaining
the distinction between a nonzero density and its zero bulk response. -/
namespace S11D4Odd
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

theorem dualCurl_divergence_zero {u : Point 4 → Vec 4} (hu : SmoothField u)
    (x : Point 4) (j : Fin 4) :
    (∑ i : Fin 4, coordDeriv i.succ (fun y => dualCurl (fieldJet u y) i j) x) = 0 := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i k : Fin 5) (r : Fin 4) :
      coordDeriv i (coordDeriv k (fun y => u y r)) x =
      coordDeriv k (coordDeriv i (fun y => u y r)) x := partial_commute (hu r) i k x
  fin_cases j
  all_goals norm_num [Fin.sum_univ_four, dualCurl, antisym,
    Matrix.cons_val_two, Matrix.cons_val_three]
  all_goals simp (disch := fun_prop) only [partial_sub, partial_zero]
  all_goals simp only [fieldJet]
  all_goals simp only [hcomm]
  all_goals ring

theorem eulerLagrange_zero {u : Point 4 → Vec 4} (hu : SmoothField u)
    (beta : ℝ) (x : Point 4) : eulerLagrange beta u x = 0 := by
  ext i
  have hp : ∀ j : Fin 4,
      coordDeriv j.succ (fun y => -(beta / 2) * dualCurl (fieldJet u y) j i) x =
        -(beta / 2) * coordDeriv j.succ (fun y => dualCurl (fieldJet u y) j i) x :=
    fun j => partial_const_mul (smooth_dualCurl hu j i) _ _ _
  change -(∑ j : Fin 5, coordDeriv j (fun y => momentum beta (fieldJet u y) j i) x) = 0
  rw [Fin.sum_univ_succ]
  simp only [momentum_eq,
    Fin.cases_zero, Fin.cases_succ, partial_zero, zero_add]
  simp_rw [hp]
  rw [← Finset.mul_sum, dualCurl_divergence_zero hu]
  simp

/-- The factor 1/2 compensates the degree-two Euler contraction. -/
def currentJet (v : Vec 4) (J : Jet 4) (i : Fin 4) : ℝ :=
  (1 / 2 : ℝ) * ∑ j : Fin 4, v j * dualCurl J i j

def boundaryCurrent (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => currentJet (u x) (fieldJet u x) i

theorem boundary_identity {u : Point 4 → Vec 4} (hu : SmoothField u) (x : Point 4) :
    (∑ i : Fin 4, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) =
      orientationJet (fieldJet u x) := by
  have hU : ∀ j : Fin 4, ContDiff ℝ ∞ (fun x => u x j) := hu
  have hM : ∀ i j : Fin 4, ContDiff ℝ ∞ (fun x => dualCurl (fieldJet u x) i j) :=
    smooth_dualCurl hu
  have hp (i j : Fin 4) :
      coordDeriv i.succ (fun y => u y j * dualCurl (fieldJet u y) i j) x =
        fieldJet u x i.succ j * dualCurl (fieldJet u x) i j +
          u x j * coordDeriv i.succ (fun y => dualCurl (fieldJet u y) i j) x :=
    partial_mul (hu j) (hM i j) _ _
  have hd (i : Fin 4) : coordDeriv i.succ (fun y => boundaryCurrent u y i) x =
      (1 / 2 : ℝ) * ∑ j : Fin 4,
        (fieldJet u x i.succ j * dualCurl (fieldJet u x) i j +
          u x j * coordDeriv i.succ (fun y => dualCurl (fieldJet u y) i j) x) := by
    unfold boundaryCurrent currentJet
    rw [partial_const_mul (by fun_prop), partial_sum _ (fun j => (hu j).mul (hM i j))]
    simp_rw [hp]
  have hc : (∑ i : Fin 4, ∑ j : Fin 4,
      u x j * coordDeriv i.succ (fun y => dualCurl (fieldJet u y) i j) x) = 0 := by
    rw [Finset.sum_comm]
    simp_rw [← Finset.mul_sum, dualCurl_divergence_zero hu]
    simp
  simp_rw [hd, Finset.sum_add_distrib]
  rw [← Finset.mul_sum, Finset.sum_add_distrib, dualCurl_contraction, hc]
  ring

theorem lagrangian_is_divergence {u : Point 4 → Vec 4} (hu : SmoothField u)
    (beta : ℝ) (x : Point 4) :
    lagrangian beta (fieldJet u x) = -(beta / 2) *
      (∑ i : Fin 4, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) := by
  rw [boundary_identity hu]
  rfl

theorem firstVariation_zero {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    deriv (relativeAction beta u h) 0 = 0 := by
  rw [relativeAction_deriv_eq_eulerLagrange hu hh]
  simp [eulerLagrange_zero hu]

theorem every_background_stationary {u : Point 4 → Vec 4}
    (hu : SmoothField u) (beta : ℝ) : ActionStationary beta u :=
  fun _ hh => firstVariation_zero hu hh beta

def witnessJet : Jet 4 := Fin.cases (fun _ => 0)
  (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,1],![0,0,0,0]])

theorem nonzero_density : orientationJet witnessJet = 1 := by
  change ((1 : ℝ) - 0) * (1 - 0) - (0 - 0) * (0 - 0) + (0 - 0) * (0 - 0) = 1
  norm_num

theorem witness_lagrangian (beta : ℝ) : lagrangian beta witnessJet = -beta / 2 := by
  rw [lagrangian, nonzero_density]
  ring

theorem nonzero_momentum : momentum 2 witnessJet 1 1 = -1 := by
  rw [momentum_eq]
  change -(2 / 2 : ℝ) * (1 - 0) = -1
  norm_num

theorem current_normalization : currentJet ![0,2,0,0] witnessJet 0 = 1 := by
  unfold currentJet
  rw [Fin.sum_univ_four]
  change (1 / 2 : ℝ) * (0 * 0 + 2 * (1 - 0) + 0 * -(0 - 0) + 0 * (0 - 0)) = 1
  norm_num

end
end S11D4Odd
