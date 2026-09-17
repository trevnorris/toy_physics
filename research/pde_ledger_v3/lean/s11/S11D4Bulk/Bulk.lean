import S11D4Bulk.Variation
import S11D4Odd.Boundary
import S11Homogeneous.Action

/-! D4C.2–D4C.3: exact bulk response and both boundary currents. -/
namespace S11D4Bulk
noncomputable section
open S10Pilot
open S11D4Odd (partial_add partial_sub partial_neg partial_const_mul partial_mul
  partial_zero partial_commute)
open scoped ContDiff
namespace Even

def gradDiv (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => ∑ j : Fin 4, coordDeriv i.succ (coordDeriv j.succ (fun y => u y j)) x

def laplacian (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => ∑ j : Fin 4, coordDeriv j.succ (coordDeriv j.succ (fun y => u y i)) x

theorem gradDiv_eq {u : Point 4 → Vec 4} (hu : SmoothField u)
    (x : Point 4) (i : Fin 4) :
    gradDiv u x i = coordDeriv i.succ (fun y => divergence (fieldJet u y)) x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  simp only [gradDiv, divergence, Fin.sum_univ_four]
  simp (disch := fun_prop) only [partial_add]
  rfl

theorem eulerLagrange_eq {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (x : Point 4) :
    eulerLagrange v u x = (v 0 + v 1) • gradDiv u x + v 2 • laplacian u x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i j : Fin 5) (r : Fin 4) :
      coordDeriv i (coordDeriv j (fun y => u y r)) x =
      coordDeriv j (coordDeriv i (fun y => u y r)) x := partial_commute (hu r) i j x
  ext i
  fin_cases i
  all_goals simp only [eulerLagrange, momentum_eq, Fin.sum_univ_succ,
    Fin.sum_univ_zero, Fin.cases_zero, Fin.cases_succ, divergence,
    Pi.smul_apply, Pi.add_apply, smul_eq_mul, gradDiv, laplacian]
  all_goals norm_num [Fin.succ] at *
  all_goals simp (disch := fun_prop) only [partial_neg, partial_add, partial_const_mul, partial_zero]
  all_goals simp only [fieldJet]
  · erw [hcomm 2 1 1, hcomm 3 1 2, hcomm 4 1 3]
    ring!
  · erw [hcomm 1 2 0, hcomm 3 2 2, hcomm 4 2 3]
    ring!
  · erw [hcomm 1 3 0, hcomm 2 3 1, hcomm 4 3 3]
    ring!
  · erw [hcomm 1 4 0, hcomm 2 4 1, hcomm 3 4 2]
    ring!

def boundaryCurrent (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => ∑ j : Fin 4,
    (u x i * fieldJet u x j.succ j - u x j * fieldJet u x j.succ i)

theorem boundary_identity {u : Point 4 → Vec 4} (hu : SmoothField u) (x : Point 4) :
    (∑ i : Fin 4, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) =
      divergence (fieldJet u x) ^ 2 - transposePair (fieldJet u x) := by
  have hU : ∀ i : Fin 4, ContDiff ℝ ∞ (fun x => u x i) := hu
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i j : Fin 5) (r : Fin 4) :
      coordDeriv i (coordDeriv j (fun y => u y r)) x =
      coordDeriv j (coordDeriv i (fun y => u y r)) x := partial_commute (hu r) i j x
  simp only [boundaryCurrent, divergence, transposePair, Fin.sum_univ_four]
  simp (disch := fun_prop) only [partial_add, partial_sub, partial_mul]
  simp only [fieldJet]
  simp only [hcomm]
  ring

def modalOperator (v : Vec 4) (k a : Vec 4) : Vec 4 :=
  (-v 2 * normSq k) • a + (-(v 0 + v 1) * dot k a) • k

theorem homogeneous_operator (v : Vec 4) (omega : ℝ) (k a : Vec 4) :
    modalOperator v k a = S11Homogeneous.modalOperator 0 (v 2)
      (v 0 + v 1 + v 2) omega k a := by
  ext i
  simp [modalOperator, S11Homogeneous.modalOperator]

theorem momentum_contraction (v : Vec 4) (omega : ℝ) (k a : Vec 4) (i : Fin 4) :
    (∑ j : Fin 5, waveCovector omega k j * momentum v (modeJet (-omega) k a) j i) =
      modalOperator v k a i := by
  simp only [momentum_eq, waveCovector, modeJet, divergence, Fin.sum_univ_succ,
    Fin.sum_univ_zero, Fin.cases_zero, Fin.cases_succ]
  fin_cases i <;> simp [modalOperator, normSq, dot, Fin.sum_univ_four] <;> ring

theorem eulerLagrange_planeWave (v : Vec 4) (omega : ℝ) (k a : Vec 4) (x : Point 4) :
    eulerLagrange v (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator v k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have he : ∀ j : Fin 5,
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum v (modeJet (-omega) k a) j i) =
      (fun y => -momentum v (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [he, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction v omega]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

theorem modal_longitudinal (v : Vec 4) (k : Vec 4) :
    modalOperator v k k = (-(v 0 + v 1 + v 2) * normSq k) • k := by
  ext i
  simp [modalOperator, normSq]
  ring

theorem modal_transverse (v : Vec 4) (k a : Vec 4) (ha : dot k a = 0) :
    modalOperator v k a = (-v 2 * normSq k) • a := by
  simp [modalOperator, ha]

theorem modal_zero_wavevector (v a : Vec 4) : modalOperator v 0 a = 0 := by
  simp [modalOperator, normSq, dot]

end Even

abbrev gradDiv := Even.gradDiv
abbrev laplacian := Even.laplacian
abbrev boundaryCurrent := Even.boundaryCurrent
abbrev modalOperator := Even.modalOperator

theorem eulerLagrange_eq {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (x : Point 4) :
    eulerLagrange v u x = (v 0 + v 1) • gradDiv u x + v 2 • laplacian u x := by
  rw [eulerLagrange_even hu, Even.eulerLagrange_eq hu]

theorem null_density_is_divergence {u : Point 4 → Vec 4} (hu : SmoothField u)
    (t beta : ℝ) (x : Point 4) :
    lagrangian ![t,-t,0,beta] (fieldJet u x) =
      -t / 2 * (∑ i : Fin 4, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) -
      beta / 2 * (∑ i : Fin 4, coordDeriv i.succ
        (fun y => S11D4Odd.boundaryCurrent u y i) x) := by
  rw [Even.boundary_identity hu, S11D4Odd.boundary_identity hu]
  simp only [lagrangian, Even.lagrangian, S11D4Odd.lagrangian, Matrix.cons_val]
  ring

theorem homogeneous_operator (v : Vec 4) (omega : ℝ) (k a : Vec 4) :
    modalOperator v k a = S11Homogeneous.modalOperator 0 (v 2)
      (v 0 + v 1 + v 2) omega k a := Even.homogeneous_operator v omega k a

theorem eulerLagrange_planeWave (v : Vec 4) (omega : ℝ) (k a : Vec 4) (x : Point 4) :
    eulerLagrange v (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator v k a := by
  rw [eulerLagrange_even (smooth_planeWave _ _ _), Even.eulerLagrange_planeWave]

theorem modal_longitudinal (v : Vec 4) (k : Vec 4) :
    modalOperator v k k = (-(v 0 + v 1 + v 2) * normSq k) • k :=
  Even.modal_longitudinal v k

theorem modal_transverse (v : Vec 4) (k a : Vec 4) (ha : dot k a = 0) :
    modalOperator v k a = (-v 2 * normSq k) • a := Even.modal_transverse v k a ha

theorem modal_zero_wavevector (v a : Vec 4) : modalOperator v 0 a = 0 :=
  Even.modal_zero_wavevector v a

end
end S11D4Bulk
