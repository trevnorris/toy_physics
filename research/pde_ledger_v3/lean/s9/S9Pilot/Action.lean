import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Tactic

/-!
S9's supplied local Lagrangian, with independent first-derivative coordinates.
This file proves finite-dimensional variation identities. It does not assume an
Euler–Lagrange equation, a spectrum, or a CAS-produced matrix as an axiom.
-/

namespace S9Pilot
noncomputable section

abbrev Vec := Fin 3 → ℝ
abbrev Point := Fin 4 → ℝ
abbrev Jet := Fin 4 → Vec

def dot (a b : Vec) : ℝ := a 0 * b 0 + a 1 * b 1 + a 2 * b 2
def normSq (a : Vec) : ℝ := dot a a
/-- First row is time derivative; rows 1,2,3 are x,y,z derivatives. -/
def jetCurl (J : Jet) : Vec := ![J 2 2 - J 3 1, J 3 0 - J 1 2, J 1 1 - J 2 0]

/-- The postulated S9 density: rho/2 |u_t|² - mu/2 |curl u|². -/
def lagrangian (rho mu : ℝ) (J : Jet) : ℝ :=
  rho / 2 * normSq (J 0) - mu / 2 * normSq (jetCurl J)

theorem dot_self_pos {a : Vec} (ha : a ≠ 0) : 0 < normSq a := by
  have hx := sq_nonneg (a 0)
  have hy := sq_nonneg (a 1)
  have hz := sq_nonneg (a 2)
  have hne : a 0 ≠ 0 ∨ a 1 ≠ 0 ∨ a 2 ≠ 0 := by
    by_contra h
    push Not at h
    apply ha
    ext i
    fin_cases i <;> simp_all
  rcases hne with h | h | h
  all_goals
    have hp := sq_pos_of_ne_zero h
    unfold normSq dot
    nlinarith

theorem lagrangian_variation (rho mu : ℝ) (J H : Jet) :
    HasDerivAt (fun s : ℝ => lagrangian rho mu (J + s • H))
      (rho * dot (J 0) (H 0) - mu * dot (jetCurl J) (jetCurl H)) 0 := by
  let C := rho * dot (J 0) (H 0) - mu * dot (jetCurl J) (jetCurl H)
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian rho mu J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const C)).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian rho mu H))
  convert! h using 1
  · funext s
    simp [C, lagrangian, normSq, dot, jetCurl]
    ring
  · simp [C]

/-- Modal derivative amplitudes, prior to any variation or spectral solve. -/
def modeJet (omega : ℝ) (k a : Vec) : Jet :=
  ![omega • a, k 0 • a, k 1 • a, k 2 • a]

def modalAction (rho mu omega : ℝ) (k a : Vec) : ℝ :=
  lagrangian rho mu (modeJet omega k a)

/-- Candidate written in coordinate form, certified by the variation theorem below. -/
def modalOperator (rho mu omega : ℝ) (k a : Vec) : Vec :=
  fun i => rho * omega ^ 2 * a i - mu * (normSq k * a i - k i * dot k a)

theorem modalAction_variation (rho mu omega : ℝ) (k a b : Vec) :
    HasDerivAt (fun s : ℝ => modalAction rho mu omega k (a + s • b))
      (dot (modalOperator rho mu omega k a) b) 0 := by
  have h := lagrangian_variation rho mu (modeJet omega k a) (modeJet omega k b)
  have heq : (fun s : ℝ => modalAction rho mu omega k (a + s • b)) =
      (fun s : ℝ => lagrangian rho mu (modeJet omega k a + s • modeJet omega k b)) := by
    funext s
    unfold modalAction
    congr 1
    ext j i
    fin_cases j <;> simp [modeJet] <;> ring
  rw [heq]
  convert! h using 1
  simp [dot, normSq, jetCurl, modeJet, modalOperator]
  ring

/-- Stationarity of the modal restriction, expressed using the actual derivative. -/
def ModalStationary (rho mu omega : ℝ) (k a : Vec) : Prop :=
  ∀ b : Vec, deriv (fun s : ℝ => modalAction rho mu omega k (a + s • b)) 0 = 0

theorem modal_stationary_iff (rho mu omega : ℝ) (k a : Vec) :
    ModalStationary rho mu omega k a ↔ modalOperator rho mu omega k a = 0 := by
  constructor
  · intro h
    have hself := h (modalOperator rho mu omega k a)
    rw [(modalAction_variation _ _ _ _ _ _).deriv] at hself
    by_contra hne
    exact (ne_of_gt (dot_self_pos hne)) hself
  · intro h b
    rw [(modalAction_variation _ _ _ _ _ _).deriv, h]
    simp [dot]

end
end S9Pilot
