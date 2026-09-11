import S10Pilot.PlaneWave

/-! The supplied one-axis inertia control. The distinguished index is arbitrary;
the other D-1 inertia coefficients remain exactly one. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

def unit (e : Fin D) : Vec D := Pi.single e 1

def kinetic (e : Fin D) (sigma : ℝ) (v : Vec D) : ℝ :=
  ∑ i, (if i = e then sigma else 1) * v i ^ 2

def lagrangian (e : Fin D) (sigma rho mu : ℝ) (J : Jet D) : ℝ :=
  rho / 2 * kinetic e sigma (J 0) - mu / 2 * S10Pilot.stiffness J

theorem kinetic_eq (e : Fin D) (sigma : ℝ) (v : Vec D) :
    kinetic e sigma v = normSq v + (sigma - 1) * v e ^ 2 := by
  have hi : ∀ i, (if i = e then sigma else 1) * v i ^ 2 =
      v i * v i + (if i = e then (sigma - 1) * v e ^ 2 else 0) := by
    intro i
    by_cases h : i = e <;> simp [h] <;> ring
  simp only [kinetic, hi, Finset.sum_add_distrib]
  simp [normSq, dot]

theorem lagrangian_eq (e : Fin D) (sigma rho mu : ℝ) (J : Jet D) :
    lagrangian e sigma rho mu J = S10Pilot.lagrangian rho mu J +
      rho / 2 * (sigma - 1) * J 0 e ^ 2 := by
  rw [lagrangian, kinetic_eq, S10Pilot.lagrangian]
  ring

def variationDensity (e : Fin D) (sigma rho mu : ℝ) (J H : Jet D) : ℝ :=
  S10Pilot.variationDensity rho mu J H + rho * (sigma - 1) * J 0 e * H 0 e

theorem lagrangian_smul (e : Fin D) (sigma rho mu s : ℝ) (J : Jet D) :
    lagrangian e sigma rho mu (s • J) = s ^ 2 * lagrangian e sigma rho mu J := by
  simp only [lagrangian_eq, S10Pilot.lagrangian_smul, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_increment (e : Fin D) (sigma rho mu s : ℝ) (J H : Jet D) :
    lagrangian e sigma rho mu (J + s • H) = lagrangian e sigma rho mu J +
      s * variationDensity e sigma rho mu J H + s ^ 2 * lagrangian e sigma rho mu H := by
  simp only [lagrangian_eq, S10Pilot.lagrangian_increment, variationDensity,
    Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (e : Fin D) (sigma rho mu : ℝ) (J H : Jet D) :
    HasDerivAt (fun s : ℝ => lagrangian e sigma rho mu (J + s • H))
      (variationDensity e sigma rho mu J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian e sigma rho mu J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity e sigma rho mu J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian e sigma rho mu H))
  convert! h using 1
  · funext s
    exact lagrangian_increment e sigma rho mu s J H
  · simp

theorem dot_unit_left (e : Fin D) (a : Vec D) : dot (unit e) a = a e := by
  simp [dot, unit, Pi.single_apply, ite_mul]

theorem dot_unit_right (e : Fin D) (a : Vec D) : dot a (unit e) = a e := by
  rw [dot_comm, dot_unit_left]

def modalAction (e : Fin D) (sigma rho mu omega : ℝ) (k a : Vec D) : ℝ :=
  lagrangian e sigma rho mu (modeJet omega k a)

def modalOperator (e : Fin D) (sigma rho mu omega : ℝ) (k a : Vec D) : Vec D :=
  S10Pilot.modalOperator rho mu omega k a + (rho * omega ^ 2 * (sigma - 1) * a e) • unit e

theorem modalAction_eq (e : Fin D) (sigma rho mu omega : ℝ) (k a : Vec D) :
    modalAction e sigma rho mu omega k a =
      rho / 2 * omega ^ 2 * (normSq a + (sigma - 1) * a e ^ 2) -
      mu / 2 * (normSq k * normSq a - dot k a ^ 2) := by
  unfold modalAction
  rw [lagrangian_eq]
  change S10Pilot.modalAction rho mu omega k a + _ = _
  rw [S10Pilot.modalAction_eq]
  simp only [modeJet, Fin.cases_zero, Pi.smul_apply, smul_eq_mul]
  ring

theorem modal_variationDensity (e : Fin D) (sigma rho mu omega : ℝ) (k a b : Vec D) :
    variationDensity e sigma rho mu (modeJet omega k a) (modeJet omega k b) =
      dot (modalOperator e sigma rho mu omega k a) b := by
  simp only [variationDensity, S10Pilot.modal_variationDensity, modalOperator,
    dot_add_left, dot_smul_left, dot_unit_left, modeJet, Fin.cases_zero,
    Pi.smul_apply, smul_eq_mul]
  ring

theorem modalAction_variation (e : Fin D) (sigma rho mu omega : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => modalAction e sigma rho mu omega k (a + s • b))
      (dot (modalOperator e sigma rho mu omega k a) b) 0 := by
  have h := lagrangian_variation e sigma rho mu (modeJet omega k a) (modeJet omega k b)
  simp_rw [modal_variationDensity] at h
  simpa only [modalAction, modeJet_add_smul] using h

def ModalStationary (e : Fin D) (sigma rho mu omega : ℝ) (k a : Vec D) : Prop :=
  ∀ b : Vec D, deriv (fun s : ℝ => modalAction e sigma rho mu omega k (a + s • b)) 0 = 0

theorem modal_stationary_iff (e : Fin D) (sigma rho mu omega : ℝ) (k a : Vec D) :
    ModalStationary e sigma rho mu omega k a ↔ modalOperator e sigma rho mu omega k a = 0 := by
  constructor
  · intro h
    have hs := h (modalOperator e sigma rho mu omega k a)
    rw [(modalAction_variation _ _ _ _ _ _ _ _).deriv] at hs
    by_contra hn
    exact (ne_of_gt (dot_self_pos hn)) hs
  · intro h b
    rw [(modalAction_variation _ _ _ _ _ _ _ _).deriv, h, zero_dot]

end
end S10Anisotropic
