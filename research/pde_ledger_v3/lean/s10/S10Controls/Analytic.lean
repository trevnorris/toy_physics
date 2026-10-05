import S10Controls.PlaneWave
import S10Pilot.Variation

/-! Analytic regularity of the momenta derived from each control action. -/

namespace S10Controls
open S10Pilot
noncomputable section
open MeasureTheory
open scoped ContDiff
variable {D : ℕ}

theorem smooth_momentum (form : Form) {u : Point D → Vec D} (hu : SmoothField u)
    (rho mu : ℝ) (j : Fin (D + 1)) (i : Fin D) :
    ContDiff ℝ ∞ (fun x => momentum form rho mu (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  cases j using Fin.cases
  · simp only [momentum_eq, Fin.cases_zero]
    fun_prop
  · cases form
    · simp only [momentum_eq, Fin.cases_succ, response]
      fun_prop
    · simp only [momentum_eq, Fin.cases_succ, response, divergence]
      split_ifs <;> fun_prop

theorem smooth_eulerLagrange (form : Form) {u : Point D → Vec D} (hu : SmoothField u)
    (rho mu : ℝ) (i : Fin D) :
    ContDiff ℝ ∞ (fun x => eulerLagrange form rho mu u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞
      (coordDeriv j (fun y => momentum form rho mu (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum form hu rho mu j i) j
  unfold eulerLagrange
  fun_prop

end
end S10Controls
