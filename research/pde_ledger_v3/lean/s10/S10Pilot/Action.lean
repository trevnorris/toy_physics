import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.LinearAlgebra.Matrix.DotProduct
import Mathlib.Tactic

/-! S10's supplied antisymmetric-derivative action in arbitrary spatial dimension.
The double sum, including its factor 1/2, is the action input. Reduced forms and
the modal operator below are conclusions checked by Lean. -/

namespace S10Pilot
noncomputable section

abbrev Vec (D : ℕ) := Fin D → ℝ
abbrev Point (D : ℕ) := Fin (D + 1) → ℝ
abbrev Jet (D : ℕ) := Fin (D + 1) → Vec D

variable {D : ℕ}

def dot (a b : Vec D) : ℝ := ∑ i, a i * b i
def normSq (a : Vec D) : ℝ := dot a a
def antisym (J : Jet D) (i j : Fin D) : ℝ := J i.succ j - J j.succ i
def stiffness (J : Jet D) : ℝ := (1 / 2 : ℝ) * ∑ i, ∑ j, antisym J i j ^ 2
def lagrangian (rho mu : ℝ) (J : Jet D) : ℝ :=
  rho / 2 * normSq (J 0) - mu / 2 * stiffness J

theorem dot_self_pos {a : Vec D} (ha : a ≠ 0) : 0 < normSq a := by
  have hn : normSq a ≠ 0 := fun h => ha (dotProduct_self_eq_zero.mp h)
  apply lt_of_le_of_ne ?_ hn.symm
  change 0 ≤ ∑ i : Fin D, a i * a i
  exact Finset.sum_nonneg (fun i _ => mul_self_nonneg (a i))

theorem dot_comm (a b : Vec D) : dot a b = dot b a := by
  simp [dot, mul_comm]

theorem dot_add_right (a b c : Vec D) : dot a (b + c) = dot a b + dot a c := by
  simp [dot, mul_add, Finset.sum_add_distrib]

theorem dot_add_left (a b c : Vec D) : dot (a + b) c = dot a c + dot b c := by
  rw [dot_comm, dot_add_right, dot_comm c a, dot_comm c b]

theorem dot_smul_right (a b : Vec D) (s : ℝ) : dot a (s • b) = s * dot a b := by
  simp [dot, Finset.mul_sum, mul_left_comm]

theorem dot_smul_left (a b : Vec D) (s : ℝ) : dot (s • a) b = s * dot a b := by
  rw [dot_comm, dot_smul_right, dot_comm b a]

theorem normSq_smul (a : Vec D) (s : ℝ) : normSq (s • a) = s ^ 2 * normSq a := by
  simp [normSq, dot_smul_left, dot_smul_right]
  ring

@[simp] theorem dot_zero (a : Vec D) : dot a 0 = 0 := by simp [dot]
@[simp] theorem zero_dot (a : Vec D) : dot 0 a = 0 := by simp [dot]

theorem sum_quadratic {ι : Type*} [Fintype ι] (a b : ι → ℝ) (s : ℝ) :
    (∑ i, (a i + s * b i) ^ 2) =
      (∑ i, a i ^ 2) + s * (2 * ∑ i, a i * b i) + s ^ 2 * ∑ i, b i ^ 2 := by
  simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem antisym_add_smul (J H : Jet D) (s : ℝ) (i j : Fin D) :
    antisym (J + s • H) i j = antisym J i j + s * antisym H i j := by
  simp [antisym]
  ring

theorem antisym_smul (J : Jet D) (s : ℝ) (i j : Fin D) :
    antisym (s • J) i j = s * antisym J i j := by
  simp [antisym, mul_sub]

theorem lagrangian_smul (rho mu s : ℝ) (J : Jet D) :
    lagrangian rho mu (s • J) = s ^ 2 * lagrangian rho mu J := by
  simp only [lagrangian, Pi.smul_apply, normSq_smul, stiffness, antisym_smul, mul_pow,
    ← Finset.mul_sum]
  ring

def variationDensity (rho mu : ℝ) (J H : Jet D) : ℝ :=
  rho * dot (J 0) (H 0) - mu / 2 * ∑ i, ∑ j, antisym J i j * antisym H i j

theorem lagrangian_increment (rho mu s : ℝ) (J H : Jet D) :
    lagrangian rho mu (J + s • H) = lagrangian rho mu J +
      s * variationDensity rho mu J H + s ^ 2 * lagrangian rho mu H := by
  unfold lagrangian variationDensity normSq dot stiffness
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, antisym_add_smul, ← sq]
  simp_rw [sum_quadratic]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  ring

theorem lagrangian_variation (rho mu : ℝ) (J H : Jet D) :
    HasDerivAt (fun s : ℝ => lagrangian rho mu (J + s • H))
      (variationDensity rho mu J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian rho mu J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity rho mu J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian rho mu H))
  convert! h using 1
  · funext s
    exact lagrangian_increment rho mu s J H
  · simp

def modeJet (omega : ℝ) (k a : Vec D) : Jet D :=
  Fin.cases (omega • a) (fun i => k i • a)

def modalAction (rho mu omega : ℝ) (k a : Vec D) : ℝ :=
  lagrangian rho mu (modeJet omega k a)

/-- Candidate formula, identified with the differentiated action below. -/
def modalOperator (rho mu omega : ℝ) (k a : Vec D) : Vec D :=
  fun i => rho * omega ^ 2 * a i - mu * (normSq k * a i - k i * dot k a)

theorem modalOperator_eq (rho mu omega : ℝ) (k a : Vec D) :
    modalOperator rho mu omega k a =
      (rho * omega ^ 2 - mu * normSq k) • a + (mu * dot k a) • k := by
  ext i
  simp [modalOperator]
  ring

theorem antisym_modeJet (omega : ℝ) (k a : Vec D) (i j : Fin D) :
    antisym (modeJet omega k a) i j = k i * a j - k j * a i := by
  simp [antisym, modeJet]

/-- The contraction identity is proved from the unreduced double sum in every D. -/
theorem antisym_mode_pair (omega : ℝ) (k a b : Vec D) :
    (∑ i, ∑ j, antisym (modeJet omega k a) i j * antisym (modeJet omega k b) i j) =
      2 * (normSq k * dot a b - dot k a * dot k b) := by
  simp only [antisym_modeJet, normSq, dot]
  have hexpand : ∀ i j : Fin D,
      (k i * a j - k j * a i) * (k i * b j - k j * b i) =
      k i * k i * (a j * b j) - (k i * b i) * (k j * a j) -
      (k i * a i) * (k j * b j) + (a i * b i) * (k j * k j) := by
    intro i j
    ring
  simp_rw [hexpand, Finset.sum_add_distrib, Finset.sum_sub_distrib,
    ← Finset.mul_sum, ← Finset.sum_mul]
  ring

theorem modalAction_eq (rho mu omega : ℝ) (k a : Vec D) :
    modalAction rho mu omega k a =
      rho / 2 * omega ^ 2 * normSq a -
      mu / 2 * (normSq k * normSq a - dot k a ^ 2) := by
  unfold modalAction lagrangian stiffness
  simp only [sq, antisym_mode_pair, modeJet, Fin.cases_zero]
  rw [normSq, dot_smul_left, dot_smul_right]
  unfold normSq
  ring

theorem modal_variationDensity (rho mu omega : ℝ) (k a b : Vec D) :
    variationDensity rho mu (modeJet omega k a) (modeJet omega k b) =
      dot (modalOperator rho mu omega k a) b := by
  have hm : dot (modalOperator rho mu omega k a) b =
      (rho * omega ^ 2 - mu * normSq k) * dot a b + mu * dot k a * dot k b := by
    rw [modalOperator_eq, dot_add_left, dot_smul_left, dot_smul_left]
  rw [hm]
  unfold variationDensity
  rw [antisym_mode_pair]
  simp only [modeJet, Fin.cases_zero, dot_smul_left, dot_smul_right]
  ring

theorem modeJet_add_smul (omega s : ℝ) (k a b : Vec D) :
    modeJet omega k (a + s • b) = modeJet omega k a + s • modeJet omega k b := by
  ext j i
  refine Fin.cases ?_ (fun j => ?_) j <;> simp [modeJet] <;> ring

theorem modalAction_variation (rho mu omega : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => modalAction rho mu omega k (a + s • b))
      (dot (modalOperator rho mu omega k a) b) 0 := by
  have h := lagrangian_variation rho mu (modeJet omega k a) (modeJet omega k b)
  simp_rw [modal_variationDensity] at h
  simpa only [modalAction, modeJet_add_smul] using h

def ModalStationary (rho mu omega : ℝ) (k a : Vec D) : Prop :=
  ∀ b : Vec D, deriv (fun s : ℝ => modalAction rho mu omega k (a + s • b)) 0 = 0

theorem modal_stationary_iff (rho mu omega : ℝ) (k a : Vec D) :
    ModalStationary rho mu omega k a ↔ modalOperator rho mu omega k a = 0 := by
  constructor
  · intro h
    have hs := h (modalOperator rho mu omega k a)
    rw [(modalAction_variation _ _ _ _ _ _).deriv] at hs
    by_contra hn
    exact (ne_of_gt (dot_self_pos hn)) hs
  · intro h b
    rw [(modalAction_variation _ _ _ _ _ _).deriv, h, zero_dot]

end
end S10Pilot
