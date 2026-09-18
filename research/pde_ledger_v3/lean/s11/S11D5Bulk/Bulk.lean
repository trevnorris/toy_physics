import S11D4Odd.Calculus
import S11D5Bulk.Variation
import S11Homogeneous.Action

/-! D5B.2–D5B.3: the actual local operator, divergence current and physical map. -/
namespace S11D5Bulk
noncomputable section
open S10Pilot
open S11D4Odd (partial_add partial_sub partial_neg partial_const_mul partial_mul
  partial_zero partial_commute)
open scoped ContDiff

def gradDiv (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => ∑ j : Fin 5, coordDeriv i.succ (coordDeriv j.succ (fun y => u y j)) x

def laplacian (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => ∑ j : Fin 5, coordDeriv j.succ (coordDeriv j.succ (fun y => u y i)) x

theorem gradDiv_eq {u : Point 5 → Vec 5} (hu : SmoothField u)
    (x : Point 5) (i : Fin 5) :
    gradDiv u x i = coordDeriv i.succ (fun y => divergence (fieldJet u y)) x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  simp only [gradDiv, divergence, S11D5Invariants.sum_five]
  simp (disch := fun_prop) only [partial_add]
  rfl

theorem eulerLagrange_eq {u : Point 5 → Vec 5} (hu : SmoothField u)
    (v : Coeff) (x : Point 5) :
    eulerLagrange v u x = (v 0 + v 1) • gradDiv u x + v 2 • laplacian u x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hdiv : ContDiff ℝ ∞ (fun y => divergence (fieldJet u y)) := by
    unfold divergence
    fun_prop
  have hp (r i : Fin 5) :
      coordDeriv r.succ (fun y => momentum v (fieldJet u y) r.succ i) x =
        -(v 0 * (if r = i then coordDeriv i.succ
            (fun y => divergence (fieldJet u y)) x else 0) +
          v 1 * coordDeriv i.succ (coordDeriv r.succ (fun y => u y r)) x +
          v 2 * coordDeriv r.succ (coordDeriv r.succ (fun y => u y i)) x) := by
    simp only [momentum_eq, Fin.cases_succ]
    by_cases h : r = i
    · subst r
      simp only [if_true]
      simp (disch := fun_prop) only [partial_neg, partial_add, partial_const_mul]
      rfl
    · simp only [h, if_false, mul_zero, zero_add]
      simp (disch := fun_prop) only [partial_neg, partial_add, partial_const_mul]
      simp only [fieldJet]
      rw [partial_commute (hu r) r.succ i.succ]
  ext i
  change -(∑ j : Fin 6, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x) = _
  rw [Fin.sum_univ_succ]
  simp only [momentum_time, partial_zero, zero_add, hp,
    Finset.sum_neg_distrib, neg_neg, Finset.sum_add_distrib, ← Finset.mul_sum]
  simp only [Finset.sum_ite_eq', Finset.mem_univ, if_true]
  rw [← gradDiv_eq hu]
  change v 0 * gradDiv u x i + v 1 * gradDiv u x i + v 2 * laplacian u x i = _
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

def boundaryCurrent (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => ∑ j : Fin 5,
    (u x i * fieldJet u x j.succ j - u x j * fieldJet u x j.succ i)

theorem boundary_identity {u : Point 5 → Vec 5} (hu : SmoothField u) (x : Point 5) :
    (∑ i : Fin 5, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) =
      divergence (fieldJet u x) ^ 2 - transposePair (fieldJet u x) := by
  have hU : ∀ i : Fin 5, ContDiff ℝ ∞ (fun x => u x i) := hu
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hcomm (i j : Fin 6) (r : Fin 5) :
      coordDeriv i (coordDeriv j (fun y => u y r)) x =
      coordDeriv j (coordDeriv i (fun y => u y r)) x := partial_commute (hu r) i j x
  simp only [boundaryCurrent, divergence, transposePair, S11D5Invariants.sum_five]
  simp (disch := fun_prop) only [partial_add, partial_sub, partial_mul]
  simp only [fieldJet]
  simp only [hcomm]
  ring

theorem null_density_is_divergence {u : Point 5 → Vec 5} (hu : SmoothField u)
    (a : ℝ) (x : Point 5) :
    lagrangian ![a,-a,0] (fieldJet u x) =
      -a / 2 * (∑ i : Fin 5, coordDeriv i.succ (fun y => boundaryCurrent u y i) x) := by
  rw [boundary_identity hu]
  simp [lagrangian]
  ring

def modalOperator (v : Coeff) (k a : Vec 5) : Vec 5 :=
  (-v 2 * normSq k) • a + (-(v 0 + v 1) * dot k a) • k

theorem homogeneous_operator (v : Coeff) (omega : ℝ) (k a : Vec 5) :
    modalOperator v k a = S11Homogeneous.modalOperator 0 (v 2)
      (v 0 + v 1 + v 2) omega k a := by
  ext i
  simp [modalOperator, S11Homogeneous.modalOperator]

theorem momentum_contraction (v : Coeff) (omega : ℝ) (k a : Vec 5) (i : Fin 5) :
    (∑ j : Fin 6, waveCovector omega k j * momentum v (modeJet (-omega) k a) j i) =
      modalOperator v k a i := by
  simp only [momentum_eq, waveCovector, modeJet, divergence, Fin.sum_univ_succ,
    Fin.sum_univ_zero, Fin.cases_zero, Fin.cases_succ]
  fin_cases i <;> simp [modalOperator, normSq, dot, S11D5Invariants.sum_five] <;> ring

theorem eulerLagrange_planeWave (v : Coeff) (omega : ℝ) (k a : Vec 5) (x : Point 5) :
    eulerLagrange v (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator v k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have he : ∀ j : Fin 6,
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

theorem modal_longitudinal (v : Coeff) (k : Vec 5) :
    modalOperator v k k = (-(v 0 + v 1 + v 2) * normSq k) • k := by
  ext i
  simp [modalOperator, normSq]
  ring

theorem modal_transverse (v : Coeff) (k a : Vec 5) (ha : dot k a = 0) :
    modalOperator v k a = (-v 2 * normSq k) • a := by
  simp [modalOperator, ha]

theorem modal_zero_wavevector (v : Coeff) (a : Vec 5) : modalOperator v 0 a = 0 := by
  simp [modalOperator, normSq, dot]

end
end S11D5Bulk
