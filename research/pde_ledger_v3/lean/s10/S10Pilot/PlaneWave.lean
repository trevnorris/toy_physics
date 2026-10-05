import S10Pilot.Action

/-! Actual real plane-wave derivatives and the action-derived local PDE in every D. -/

namespace S10Pilot
noncomputable section

variable {D : ℕ}

def axis (j : Fin (D + 1)) : Point D := Pi.single j 1
def basisJet (j : Fin (D + 1)) (i : Fin D) : Jet D := Pi.single j (Pi.single i 1)

theorem basisJet_apply (j a : Fin (D + 1)) (i b : Fin D) :
    basisJet j i a b = if a = j then (if b = i then 1 else 0) else 0 := by
  simp [basisJet, Pi.single_apply, ite_apply]

def phase (q x : Point D) : ℝ := dot q x
def coordDeriv (j : Fin (D + 1)) (f : Point D → ℝ) (x : Point D) : ℝ :=
  deriv (fun s : ℝ => f (x + s • axis j)) 0
def fieldJet (u : Point D → Vec D) (x : Point D) : Jet D :=
  fun j i => coordDeriv j (fun x => u x i) x
def waveCovector (omega : ℝ) (k : Vec D) : Point D := Fin.cases (-omega) k
def planeWave (omega : ℝ) (k a : Vec D) (x : Point D) : Vec D :=
  fun i => a i * Real.cos (phase (waveCovector omega k) x)

theorem phase_shift (q x : Point D) (j : Fin (D + 1)) (s : ℝ) :
    phase q (x + s • axis j) = phase q x + s * q j := by
  unfold phase
  rw [dot_add_right, dot_smul_right]
  simp [dot, axis, Pi.single_apply]

theorem partial_const_cos (q x : Point D) (j : Fin (D + 1)) (c : ℝ) :
    coordDeriv j (fun y => c * Real.cos (phase q y)) x =
      -c * Real.sin (phase q x) * q j := by
  unfold coordDeriv
  simp_rw [phase_shift]
  have h := (((hasDerivAt_const (0 : ℝ) (phase q x)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (q j))).cos).const_mul c
  convert! h.deriv using 1
  simp
  ring

theorem partial_const_sin (q x : Point D) (j : Fin (D + 1)) (c : ℝ) :
    coordDeriv j (fun y => c * Real.sin (phase q y)) x =
      c * Real.cos (phase q x) * q j := by
  unfold coordDeriv
  simp_rw [phase_shift]
  have h := (((hasDerivAt_const (0 : ℝ) (phase q x)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (q j))).sin).const_mul c
  convert! h.deriv using 1
  simp
  ring

theorem fieldJet_planeWave (omega : ℝ) (k a : Vec D) (x : Point D) :
    fieldJet (planeWave omega k a) x =
      (-Real.sin (phase (waveCovector omega k) x)) • modeJet (-omega) k a := by
  ext j i
  simp only [fieldJet, planeWave, partial_const_cos]
  refine Fin.cases ?_ (fun j => ?_) j <;> simp [waveCovector, modeJet] <;> ring

/-- Exact pointwise density identity; no infinite-volume action integral is taken. -/
theorem planeWave_lagrangian (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    lagrangian rho mu (fieldJet (planeWave omega k a) x) =
      Real.sin (phase (waveCovector omega k) x) ^ 2 * modalAction rho mu omega k a := by
  rw [fieldJet_planeWave]
  rw [lagrangian_smul]
  change (-Real.sin (phase (waveCovector omega k) x)) ^ 2 * modalAction rho mu (-omega) k a = _
  simp [modalAction_eq]

/-- Conjugate derivative obtained from the density by an actual derivative. -/
def momentum (rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) : ℝ :=
  deriv (fun s : ℝ => lagrangian rho mu (J + s • basisJet j i)) 0

theorem momentum_eq (rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) :
    momentum rho mu J j i =
      Fin.cases (rho * J 0 i) (fun m => -mu * antisym J m i) j := by
  unfold momentum
  rw [(lagrangian_variation rho mu J (basisJet j i)).deriv]
  refine Fin.cases ?_ (fun m => ?_) j
  · simp [variationDensity, dot, antisym, basisJet_apply]
  · have hm : (0 : Fin (D + 1)) ≠ m.succ := (Fin.succ_ne_zero m).symm
    simp [variationDensity, dot, antisym, basisJet_apply, mul_sub,
      Finset.sum_sub_distrib, mul_ite, hm]
    ring

theorem momentum_smul (rho mu s : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) :
    momentum rho mu (s • J) j i = s * momentum rho mu J j i := by
  rw [momentum_eq, momentum_eq]
  refine Fin.cases ?_ (fun j => ?_) j <;> simp [antisym] <;> ring

/-- First-order, field-independent density: dL/du - sum_j d_j(dL/dJ_ji). -/
def eulerLagrange (rho mu : ℝ) (u : Point D → Vec D) (x : Point D) : Vec D :=
  fun i => -∑ j : Fin (D + 1), coordDeriv j (fun y => momentum rho mu (fieldJet u y) j i) x

theorem momentum_contraction (rho mu omega : ℝ) (k a : Vec D) (i : Fin D) :
    (∑ j : Fin (D + 1), waveCovector omega k j * momentum rho mu (modeJet (-omega) k a) j i) =
      modalOperator rho mu omega k a i := by
  simp only [Fin.sum_univ_succ, waveCovector, Fin.cases_zero, Fin.cases_succ,
    momentum_eq, modeJet, antisym, Pi.smul_apply, smul_eq_mul]
  simp only [mul_sub, mul_neg, neg_mul, Finset.sum_sub_distrib, Finset.sum_neg_distrib]
  have h₁ : (∑ j : Fin D, k j * (mu * (k j * a i))) = mu * normSq k * a i := by
    simp only [normSq, dot, Finset.mul_sum, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro j _
    ring
  have h₂ : (∑ j : Fin D, k j * (mu * (k i * a j))) = mu * k i * dot k a := by
    simp only [dot, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro j _
    ring
  rw [h₁, h₂]
  simp only [modalOperator]
  ring

/-- The independently differentiated position-space expression yields the same matrix. -/
theorem eulerLagrange_planeWave (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    eulerLagrange rho mu (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator rho mu omega k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have hrearrange : ∀ j : Fin (D + 1),
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum rho mu (modeJet (-omega) k a) j i) =
      (fun y => -momentum rho mu (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [hrearrange, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- A global plane wave solves the local EL expression iff its amplitude is stationary. -/
theorem planeWave_solves_iff (rho mu omega : ℝ) (k a : Vec D) :
    (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
      ModalStationary rho mu omega k a := by
  rw [modal_stationary_iff]
  constructor
  · intro h
    simpa [eulerLagrange_planeWave, phase, dot] using h 0
  · intro h x
    simp [eulerLagrange_planeWave, h]

end
end S10Pilot
