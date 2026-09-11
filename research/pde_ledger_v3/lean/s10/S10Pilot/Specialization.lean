import S10Pilot.VariationalCertificate
import S9Pilot.VariationalCertificate

/-! Equality of the D=3 definitions with the original S9 pilot, including the
supplied density and the integrated stationarity predicate. -/

namespace S10Pilot
noncomputable section

theorem dot_three (a b : Vec 3) : dot a b = S9Pilot.dot a b := by
  simp [dot, S9Pilot.dot, Fin.sum_univ_succ]
  ring

theorem normSq_three (a : Vec 3) : normSq a = S9Pilot.normSq a := dot_three a a

theorem stiffness_three (J : Jet 3) : stiffness J = S9Pilot.normSq (S9Pilot.jetCurl J) := by
  simp [stiffness, antisym, S9Pilot.normSq, S9Pilot.dot, S9Pilot.jetCurl, Fin.sum_univ_succ]
  ring

theorem lagrangian_three (rho mu : ℝ) (J : Jet 3) :
    lagrangian rho mu J = S9Pilot.lagrangian rho mu J := by
  rw [lagrangian, normSq_three, stiffness_three]
  rfl

theorem modeJet_three (omega : ℝ) (k a : Vec 3) :
    modeJet omega k a = S9Pilot.modeJet omega k a := by
  ext j i
  fin_cases j <;> rfl

theorem modalAction_three (rho mu omega : ℝ) (k a : Vec 3) :
    modalAction rho mu omega k a = S9Pilot.modalAction rho mu omega k a := by
  simp only [modalAction, lagrangian_three, modeJet_three, S9Pilot.modalAction]

theorem modalOperator_three (rho mu omega : ℝ) (k a : Vec 3) :
    modalOperator rho mu omega k a = S9Pilot.modalOperator rho mu omega k a := by
  ext i
  simp only [modalOperator, normSq_three, dot_three, S9Pilot.modalOperator]

theorem coneValue_three (rho mu : ℝ) (k : Vec 3) :
    coneValue rho mu k = S9Pilot.coneValue rho mu k := by
  simp only [coneValue, normSq_three, S9Pilot.coneValue]

theorem transverseSpace_three (k : Vec 3) : transverseSpace k = S9Pilot.transverseSpace k := by
  ext a
  change dot k a = 0 ↔ S9Pilot.dot k a = 0
  rw [dot_three]

theorem longitudinalSpace_three (k : Vec 3) :
    longitudinalSpace k = S9Pilot.longitudinalSpace k := rfl

theorem phase_three (q x : Point 3) : phase q x = S9Pilot.phase q x := by
  simp [phase, dot, S9Pilot.phase, Fin.sum_univ_succ]
  ring

theorem waveCovector_three (omega : ℝ) (k : Vec 3) :
    waveCovector omega k = S9Pilot.waveCovector omega k := by
  ext j
  fin_cases j <;> rfl

theorem planeWave_three (omega : ℝ) (k a : Vec 3) :
    planeWave omega k a = S9Pilot.planeWave omega k a := by
  ext x i
  simp only [planeWave, phase_three, waveCovector_three, S9Pilot.planeWave]

theorem coordDeriv_three (j : Fin 4) (f : Point 3 → ℝ) (x : Point 3) :
    coordDeriv j f x = S9Pilot.coordDeriv j f x := rfl

theorem fieldJet_three (u : Point 3 → Vec 3) (x : Point 3) :
    fieldJet u x = S9Pilot.fieldJet u x := rfl

theorem momentum_three (rho mu : ℝ) (J : Jet 3) (j : Fin 4) (i : Fin 3) :
    momentum rho mu J j i = S9Pilot.momentum rho mu J j i := by
  unfold momentum S9Pilot.momentum
  congr 1
  funext s
  exact lagrangian_three rho mu (J + s • basisJet j i)

theorem eulerLagrange_three (rho mu : ℝ) (u : Point 3 → Vec 3) (x : Point 3) :
    eulerLagrange rho mu u x = S9Pilot.eulerLagrange rho mu u x := by
  ext i
  change -(∑ j : Fin 4, coordDeriv j (fun y => momentum rho mu (fieldJet u y) j i) x) =
    -(∑ j : Fin 4, S9Pilot.coordDeriv j (fun y =>
      S9Pilot.momentum rho mu (S9Pilot.fieldJet u y) j i) x)
  congr 1
  apply Finset.sum_congr rfl
  intro j _
  rw [coordDeriv_three]
  congr 1
  funext y
  exact momentum_three rho mu (fieldJet u y) j i

theorem relativeAction_three (rho mu s : ℝ) (u h : Point 3 → Vec 3) :
    relativeAction rho mu u h s = S9Pilot.relativeAction rho mu u h s := by
  simp only [relativeAction, fieldJet_three, lagrangian_three, S9Pilot.relativeAction]

theorem actionStationary_three (rho mu : ℝ) (u : Point 3 → Vec 3) :
    ActionStationary rho mu u ↔ S9Pilot.ActionStationary rho mu u := by
  have heq : ∀ h, relativeAction rho mu u h = S9Pilot.relativeAction rho mu u h :=
    fun h => funext (fun s => relativeAction_three rho mu s u h)
  simp only [ActionStationary, heq, S9Pilot.ActionStationary]
  rfl

end
end S10Pilot
