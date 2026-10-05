import S11D3Bulk.Calculus
import S11D3Bulk.Variation
import S11Homogeneous.Action

/-! K2–K3: the actual local operator, divergence current and physical map. -/
namespace S11D3Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

def gradDiv (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun i => ∑ j : Fin 3, coordDeriv i.succ (coordDeriv j.succ (fun y => u y j)) x

def laplacian (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun i => ∑ j : Fin 3, coordDeriv j.succ (coordDeriv j.succ (fun y => u y i)) x

theorem gradDiv_eq {u : Point 3 → Vec 3} (hu : SmoothField u)
    (x : Point 3) (i : Fin 3) :
    gradDiv u x i = coordDeriv i.succ (fun y => divergence (fieldJet u y)) x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  simp only [gradDiv, divergence, Fin.sum_univ_three]
  simp (disch := fun_prop) only [partial_add]
  rfl

theorem eulerLagrange_eq {u : Point 3 → Vec 3} (hu : SmoothField u)
    (v : Vec 3) (x : Point 3) :
    eulerLagrange v u x = (v 0 + v 1) • gradDiv u x + v 2 • laplacian u x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i j : Fin 4) (r : Fin 3) :
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
  · erw [hcomm 2 1 1, hcomm 3 1 2]
    ring!
  · erw [hcomm 1 2 0, hcomm 3 2 2]
    ring!
  · erw [hcomm 1 3 0, hcomm 2 3 1]
    ring!

def boundaryCurrent (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun i => ∑ j : Fin 3,
    (u x i * fieldJet u x j.succ j - u x j * fieldJet u x j.succ i)

theorem boundary_identity {u : Point 3 → Vec 3} (hu : SmoothField u) (x : Point 3) :
    (∑ i : Fin 3, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) =
      divergence (fieldJet u x) ^ 2 - transposePair (fieldJet u x) := by
  have hU : ∀ i : Fin 3, ContDiff ℝ ∞ (fun x => u x i) := hu
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i j : Fin 4) (r : Fin 3) :
      coordDeriv i (coordDeriv j (fun y => u y r)) x =
      coordDeriv j (coordDeriv i (fun y => u y r)) x := partial_commute (hu r) i j x
  simp only [boundaryCurrent, divergence, transposePair, Fin.sum_univ_three]
  simp (disch := fun_prop) only [partial_add, partial_sub, partial_mul]
  simp only [fieldJet]
  simp only [hcomm]
  ring

theorem null_density_is_divergence {u : Point 3 → Vec 3} (hu : SmoothField u)
    (a : ℝ) (x : Point 3) :
    lagrangian ![a,-a,0] (fieldJet u x) =
      -a / 2 * (∑ i : Fin 3, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) := by
  rw [boundary_identity hu]
  simp [lagrangian]
  ring

def modalOperator (v : Vec 3) (k a : Vec 3) : Vec 3 :=
  (-v 2 * normSq k) • a + (-(v 0 + v 1) * dot k a) • k

theorem homogeneous_operator (v : Vec 3) (omega : ℝ) (k a : Vec 3) :
    modalOperator v k a = S11Homogeneous.modalOperator 0 (v 2)
      (v 0 + v 1 + v 2) omega k a := by
  ext i
  simp [modalOperator, S11Homogeneous.modalOperator]

theorem momentum_contraction (v : Vec 3) (omega : ℝ) (k a : Vec 3) (i : Fin 3) :
    (∑ j : Fin 4, waveCovector omega k j * momentum v (modeJet (-omega) k a) j i) =
      modalOperator v k a i := by
  simp only [momentum_eq, waveCovector, modeJet, divergence, Fin.sum_univ_succ,
    Fin.sum_univ_zero, Fin.cases_zero, Fin.cases_succ]
  fin_cases i <;> simp [modalOperator, normSq, dot, Fin.sum_univ_three] <;> ring

theorem eulerLagrange_planeWave (v : Vec 3) (omega : ℝ) (k a : Vec 3) (x : Point 3) :
    eulerLagrange v (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator v k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have he : ∀ j : Fin 4,
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

theorem modal_longitudinal (v : Vec 3) (k : Vec 3) :
    modalOperator v k k = (-(v 0 + v 1 + v 2) * normSq k) • k := by
  ext i
  simp [modalOperator, normSq]
  ring

theorem modal_transverse (v : Vec 3) (k a : Vec 3) (ha : dot k a = 0) :
    modalOperator v k a = (-v 2 * normSq k) • a := by
  simp [modalOperator, ha]

theorem modal_zero_wavevector (v a : Vec 3) : modalOperator v 0 a = 0 := by
  simp [modalOperator, normSq, dot]

end
end S11D3Bulk
