import S10Controls.Action

/-! Coordinate variation and actual plane-wave evaluation for the two control actions. -/

namespace S10Controls
open S10Pilot
noncomputable section
variable {D : ℕ}

/-- Exact pointwise density identity; no infinite-volume action integral is taken. -/
theorem planeWave_lagrangian (form : Form) (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    lagrangian form rho mu (fieldJet (planeWave omega k a) x) =
      Real.sin (phase (waveCovector omega k) x) ^ 2 * modalAction form rho mu omega k a := by
  rw [fieldJet_planeWave]
  rw [lagrangian_smul form]
  change (-Real.sin (phase (waveCovector omega k) x)) ^ 2 * modalAction form rho mu (-omega) k a = _
  simp [modalAction_eq form]

def response (form : Form) (J : Jet D) (j i : Fin D) : ℝ := match form with
  | .fullGradient => J j.succ i
  | .divergenceOnly => if j = i then divergence J else 0

/-- Momentum is an actual derivative of the selected density. -/
def momentum (form : Form) (rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) : ℝ :=
  deriv (fun s : ℝ => lagrangian form rho mu (J + s • basisJet j i)) 0

theorem momentum_eq (form : Form) (rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) :
    momentum form rho mu J j i =
      Fin.cases (rho * J 0 i) (fun m => -mu * response form J m i) j := by
  unfold momentum
  rw [(lagrangian_variation form rho mu J (basisJet j i)).deriv]
  cases form
  · refine Fin.cases ?_ (fun m => ?_) j
    · simp [variationDensity, stiffnessPair, dot, basisJet_apply]
    · have hm : (0 : Fin (D + 1)) ≠ m.succ := (Fin.succ_ne_zero m).symm
      simp [variationDensity, stiffnessPair, dot, basisJet_apply, response, hm, mul_ite]
  · refine Fin.cases ?_ (fun m => ?_) j
    · simp [variationDensity, stiffnessPair, divergence, dot, basisJet_apply]
    · have hm : (0 : Fin (D + 1)) ≠ m.succ := (Fin.succ_ne_zero m).symm
      by_cases hi : m = i
      · subst i
        simp [variationDensity, stiffnessPair, divergence, dot, basisJet_apply, response, hm]
      · simp [variationDensity, stiffnessPair, divergence, dot, basisJet_apply, response, hm, hi]

theorem momentum_smul (form : Form) (rho mu s : ℝ) (J : Jet D)
    (j : Fin (D + 1)) (i : Fin D) :
    momentum form rho mu (s • J) j i = s * momentum form rho mu J j i := by
  rw [momentum_eq, momentum_eq]
  refine Fin.cases ?_ (fun j => ?_) j
  · simp
    ring
  · cases form
    · simp [response]
      ring
    · by_cases hi : j = i
      · simp [response, divergence_smul, hi]
        ring
      · simp [response, hi]

def eulerLagrange (form : Form) (rho mu : ℝ) (u : Point D → Vec D) (x : Point D) : Vec D :=
  fun i => -∑ j : Fin (D + 1), coordDeriv j (fun y => momentum form rho mu (fieldJet u y) j i) x

theorem momentum_contraction (form : Form) (rho mu omega : ℝ) (k a : Vec D) (i : Fin D) :
    (∑ j : Fin (D + 1), waveCovector omega k j *
      momentum form rho mu (modeJet (-omega) k a) j i) =
      modalOperator form rho mu omega k a i := by
  cases form
  · simp only [Fin.sum_univ_succ, waveCovector, momentum_eq, response, modeJet,
      Fin.cases_zero, Fin.cases_succ, Pi.smul_apply, smul_eq_mul]
    have hs : (∑ j : Fin D, k j * (-mu * (k j * a i))) = -mu * normSq k * a i := by
      simp only [normSq, dot, Finset.mul_sum, Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro j _
      ring
    rw [hs]
    simp only [modalOperator, Pi.smul_apply, smul_eq_mul]
    ring
  · simp [Fin.sum_univ_succ, waveCovector, momentum_eq, response, divergence_modeJet,
      modeJet, modalOperator, mul_ite]
    ring

/-- The independently differentiated position-space expression yields the same matrix. -/
theorem eulerLagrange_planeWave (form : Form) (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    eulerLagrange form rho mu (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator form rho mu omega k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul form]
  have hrearrange : ∀ j : Fin (D + 1),
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum form rho mu (modeJet (-omega) k a) j i) =
      (fun y => -momentum form rho mu (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [hrearrange, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction form]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- A global plane wave solves the local EL expression iff its amplitude is stationary. -/
theorem planeWave_solves_iff (form : Form) (rho mu omega : ℝ) (k a : Vec D) :
    (∀ x, eulerLagrange form rho mu (planeWave omega k a) x = 0) ↔
      ModalStationary form rho mu omega k a := by
  rw [modal_stationary_iff form]
  constructor
  · intro h
    simpa [eulerLagrange_planeWave form, phase, dot] using h 0
  · intro h x
    simp [eulerLagrange_planeWave form, h]

end
end S10Controls
