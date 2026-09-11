import S10Anisotropic.PlaneWave
import S10Pilot.Variation

/-! Analytic regularity of the momenta derived from the anisotropic action. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
open MeasureTheory
open scoped ContDiff
variable {D : ℕ}

theorem smooth_momentum (e : Fin D) (sigma : ℝ) {u : Point D → Vec D} (hu : SmoothField u)
    (rho mu : ℝ) (j : Fin (D + 1)) (i : Fin D) :
    ContDiff ℝ ∞ (fun x => momentum e sigma rho mu (fieldJet u x) j i) := by
  have hJ : ContDiff ℝ ∞ (fun x => fieldJet u x 0 e) := smooth_coordDeriv (hu e) 0
  have hp := S10Pilot.smooth_momentum hu rho mu j i
  simp only [momentum_eq]
  split_ifs <;> fun_prop

theorem smooth_eulerLagrange (e : Fin D) (sigma : ℝ) {u : Point D → Vec D} (hu : SmoothField u)
    (rho mu : ℝ) (i : Fin D) :
    ContDiff ℝ ∞ (fun x => eulerLagrange e sigma rho mu u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞
      (coordDeriv j (fun y => momentum e sigma rho mu (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum e sigma hu rho mu j i) j
  unfold eulerLagrange
  fun_prop

end
end S10Anisotropic
