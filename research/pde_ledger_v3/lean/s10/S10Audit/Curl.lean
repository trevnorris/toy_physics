import S10Audit.Packages
import Physlib.Mathematics.LeviCivita.Basic

/-! Q7: construct the ordinary curl by epsilon contraction, then compare
the actual package stiffness operand, independently of its action coefficient. -/

namespace S10Audit
open S10Pilot
noncomputable section
set_option backward.isDefEq.respectTransparency false

def epsilonCurl (J : Jet 3) : Vec 3 := fun i =>
  ∑ j, ∑ k, (leviCivitaSymbol ![i, j, k] : ℝ) * J j.succ k

theorem epsilonCurl_eq_jetCurl (J : Jet 3) : epsilonCurl J = S9Pilot.jetCurl J := by
  ext i
  fin_cases i <;>
    simp only [epsilonCurl, S9Pilot.jetCurl, Fin.isValue, Fin.zero_eta,
      Fin.mk_one, Fin.reduceFinMk, Fin.sum_univ_three,
      leviCivitaSymbol_eq_det, Matrix.det_fin_three, KroneckerDelta.kroneckerDelta,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons] <;> norm_num <;>
    try simp only [show (Fin.succ (2 : Fin 3)) = (3 : Fin 4) from rfl]
  all_goals ring

def curlSquared (J : Jet 3) : ℝ := normSq (epsilonCurl J)

theorem curlSquared_eq_stiffness (J : Jet 3) : curlSquared J = S10Pilot.stiffness J := by
  rw [curlSquared, epsilonCurl_eq_jetCurl, normSq_three, stiffness_three]

def q7Difference (p : Package) (J : Jet 3) : ℝ := packageStiffness p J - curlSquared J

theorem main_q7 (J : Jet 3) : q7Difference .main J = 0 := by
  simp [q7Difference, packageStiffness, curlSquared_eq_stiffness]

theorem signFlip_q7 (J : Jet 3) : q7Difference .signFlip J = 0 := main_q7 J

theorem anisotropic_q7 (J : Jet 3) : q7Difference .anisotropic J = 0 := main_q7 J

theorem coefficientScale_q7 (J : Jet 3) : q7Difference .coefficientScale J = 0 := main_q7 J

theorem fullGradient_q7 (J : Jet 3) :
    q7Difference .fullGradient J = ∑ i, ∑ j, J i.succ j * J j.succ i := by
  simp [q7Difference, packageStiffness, curlSquared_eq_stiffness,
    S10Controls.stiffness, S10Controls.fullGradientStiffness,
    S10Pilot.stiffness, antisym, Fin.sum_univ_succ]
  ring

theorem divergenceOnly_q7 (J : Jet 3) :
    q7Difference .divergenceOnly J = (∑ i : Fin 3, J i.succ i) ^ 2 -
      (∑ i : Fin 3, ∑ j : Fin 3, J i.succ j ^ 2) + ∑ i, ∑ j, J i.succ j * J j.succ i := by
  simp [q7Difference, packageStiffness, curlSquared_eq_stiffness,
    S10Controls.stiffness, S10Controls.divergenceOnlyStiffness, S10Controls.divergence,
    S10Pilot.stiffness, antisym, Fin.sum_univ_succ]
  ring

/-- A symmetric diagonal gradient has zero curl but nonzero control stiffness. -/
theorem q7_control_counterexample :
    let J : Jet 3 := ![![0, 0, 0], ![1, 0, 0], ![0, 0, 0], ![0, 0, 0]]
    curlSquared J = 0 ∧ q7Difference .fullGradient J = 1 ∧
      q7Difference .divergenceOnly J = 1 := by
  norm_num [curlSquared_eq_stiffness, q7Difference, packageStiffness,
    S10Controls.stiffness, S10Controls.fullGradientStiffness,
    S10Controls.divergenceOnlyStiffness, S10Controls.divergence,
    S10Pilot.stiffness, antisym, Fin.sum_univ_succ]

end
end S10Audit
