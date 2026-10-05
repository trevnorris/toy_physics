import S10Pilot.PlaneWave

/-! The two stiffness-form controls, each entering through its own supplied
spatial derivative density. Shared coordinates and calculus come from S10Pilot. -/

namespace S10Controls
open S10Pilot
noncomputable section
variable {D : ℕ}

inductive Form where
  | fullGradient
  | divergenceOnly
  deriving DecidableEq

def fullGradientStiffness (J : Jet D) : ℝ := ∑ i : Fin D, ∑ j : Fin D, J i.succ j ^ 2
def divergence (J : Jet D) : ℝ := ∑ i : Fin D, J i.succ i
def divergenceOnlyStiffness (J : Jet D) : ℝ := divergence J ^ 2

def stiffness (form : Form) (J : Jet D) : ℝ := match form with
  | .fullGradient => fullGradientStiffness J
  | .divergenceOnly => divergenceOnlyStiffness J

def lagrangian (form : Form) (rho mu : ℝ) (J : Jet D) : ℝ :=
  rho / 2 * normSq (J 0) - mu / 2 * stiffness form J

def stiffnessPair (form : Form) (J H : Jet D) : ℝ := match form with
  | .fullGradient => ∑ i : Fin D, ∑ j : Fin D, J i.succ j * H i.succ j
  | .divergenceOnly => divergence J * divergence H

def variationDensity (form : Form) (rho mu : ℝ) (J H : Jet D) : ℝ :=
  rho * dot (J 0) (H 0) - mu * stiffnessPair form J H

theorem divergence_add_smul (J H : Jet D) (s : ℝ) :
    divergence (J + s • H) = divergence J + s * divergence H := by
  simp [divergence, Finset.sum_add_distrib, Finset.mul_sum]

theorem divergence_smul (J : Jet D) (s : ℝ) : divergence (s • J) = s * divergence J := by
  simp [divergence, Finset.mul_sum]

theorem lagrangian_smul (form : Form) (rho mu s : ℝ) (J : Jet D) :
    lagrangian form rho mu (s • J) = s ^ 2 * lagrangian form rho mu J := by
  cases form
  · simp only [lagrangian, stiffness, fullGradientStiffness, Pi.smul_apply,
      smul_eq_mul, normSq_smul, mul_pow, ← Finset.mul_sum]
    ring
  · simp only [lagrangian, stiffness, divergenceOnlyStiffness, Pi.smul_apply,
      normSq_smul, divergence_smul]
    ring

theorem lagrangian_increment (form : Form) (rho mu s : ℝ) (J H : Jet D) :
    lagrangian form rho mu (J + s • H) = lagrangian form rho mu J +
      s * variationDensity form rho mu J H + s ^ 2 * lagrangian form rho mu H := by
  cases form
  · unfold lagrangian variationDensity stiffness stiffnessPair fullGradientStiffness normSq dot
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, ← sq]
    simp_rw [sum_quadratic]
    simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
    ring
  · unfold lagrangian variationDensity stiffness stiffnessPair divergenceOnlyStiffness normSq dot
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, divergence_add_smul, ← sq]
    rw [sum_quadratic]
    ring

theorem lagrangian_variation (form : Form) (rho mu : ℝ) (J H : Jet D) :
    HasDerivAt (fun s : ℝ => lagrangian form rho mu (J + s • H))
      (variationDensity form rho mu J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian form rho mu J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity form rho mu J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian form rho mu H))
  convert! h using 1
  · funext s
    exact lagrangian_increment form rho mu s J H
  · simp

theorem divergence_modeJet (omega : ℝ) (k a : Vec D) :
    divergence (modeJet omega k a) = dot k a := by simp [divergence, modeJet, dot]

theorem fullGradient_mode_pair (omega : ℝ) (k a b : Vec D) :
    stiffnessPair .fullGradient (modeJet omega k a) (modeJet omega k b) =
      normSq k * dot a b := by
  simp only [stiffnessPair, modeJet, Fin.cases_succ, Pi.smul_apply, smul_eq_mul]
  have he : ∀ i j : Fin D, (k i * a j) * (k i * b j) = (k i * k i) * (a j * b j) := by
    intro i j
    ring
  simp_rw [he, ← Finset.mul_sum, ← Finset.sum_mul]
  rfl

theorem stiffnessPair_self (form : Form) (J : Jet D) :
    stiffnessPair form J J = stiffness form J := by
  cases form <;> simp [stiffnessPair, stiffness, fullGradientStiffness, divergenceOnlyStiffness, sq]

def modalAction (form : Form) (rho mu omega : ℝ) (k a : Vec D) : ℝ :=
  lagrangian form rho mu (modeJet omega k a)

/-- Candidate operators; identification with the action's derivative is proved below. -/
def modalOperator (form : Form) (rho mu omega : ℝ) (k a : Vec D) : Vec D := match form with
  | .fullGradient => (rho * omega ^ 2 - mu * normSq k) • a
  | .divergenceOnly => (rho * omega ^ 2) • a + (-mu * dot k a) • k

theorem modalAction_eq (form : Form) (rho mu omega : ℝ) (k a : Vec D) :
    modalAction form rho mu omega k a = rho / 2 * omega ^ 2 * normSq a - mu / 2 *
      (match form with
       | .fullGradient => normSq k * normSq a
       | .divergenceOnly => dot k a ^ 2) := by
  unfold modalAction lagrangian
  rw [← stiffnessPair_self]
  cases form
  · rw [fullGradient_mode_pair]
    simp only [modeJet, Fin.cases_zero, normSq_smul]
    unfold normSq
    ring
  · simp only [stiffnessPair, divergence_modeJet, modeJet, Fin.cases_zero, normSq_smul]
    ring

theorem modal_variationDensity (form : Form) (rho mu omega : ℝ) (k a b : Vec D) :
    variationDensity form rho mu (modeJet omega k a) (modeJet omega k b) =
      dot (modalOperator form rho mu omega k a) b := by
  unfold variationDensity
  cases form
  · rw [fullGradient_mode_pair]
    simp only [modeJet, Fin.cases_zero, dot_smul_left, dot_smul_right, modalOperator]
    ring
  · simp only [stiffnessPair, divergence_modeJet, modeJet, Fin.cases_zero,
      modalOperator, dot_add_left, dot_smul_left, dot_smul_right]
    ring

theorem modalAction_variation (form : Form) (rho mu omega : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => modalAction form rho mu omega k (a + s • b))
      (dot (modalOperator form rho mu omega k a) b) 0 := by
  have h := lagrangian_variation form rho mu (modeJet omega k a) (modeJet omega k b)
  simp_rw [modal_variationDensity] at h
  simpa only [modalAction, modeJet_add_smul] using h

def ModalStationary (form : Form) (rho mu omega : ℝ) (k a : Vec D) : Prop :=
  ∀ b : Vec D, deriv (fun s : ℝ => modalAction form rho mu omega k (a + s • b)) 0 = 0

theorem modal_stationary_iff (form : Form) (rho mu omega : ℝ) (k a : Vec D) :
    ModalStationary form rho mu omega k a ↔ modalOperator form rho mu omega k a = 0 := by
  constructor
  · intro h
    have hs := h (modalOperator form rho mu omega k a)
    rw [(modalAction_variation _ _ _ _ _ _ _).deriv] at hs
    by_contra hn
    exact (ne_of_gt (dot_self_pos hn)) hs
  · intro h b
    rw [(modalAction_variation _ _ _ _ _ _ _).deriv, h, zero_dot]

end
end S10Controls
