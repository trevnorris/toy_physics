import S10Anisotropic.Action

/-! Coordinate variation and actual plane-wave evaluation for the one-axis inertia action. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

/-- Exact pointwise density identity; no infinite-volume action integral is taken. -/
theorem planeWave_lagrangian (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    lagrangian e sigma rho mu (fieldJet (planeWave omega k a) x) =
      Real.sin (phase (waveCovector omega k) x) ^ 2 * modalAction e sigma rho mu omega k a := by
  rw [fieldJet_planeWave]
  rw [lagrangian_smul e sigma]
  change (-Real.sin (phase (waveCovector omega k) x)) ^ 2 * modalAction e sigma rho mu (-omega) k a = _
  simp [modalAction_eq e sigma]

/-- Momentum is an actual derivative of the anisotropic density. -/
def momentum (e : Fin D) (sigma rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) : ℝ :=
  deriv (fun s : ℝ => lagrangian e sigma rho mu (J + s • basisJet j i)) 0

theorem momentum_eq (e : Fin D) (sigma rho mu : ℝ) (J : Jet D) (j : Fin (D + 1)) (i : Fin D) :
    momentum e sigma rho mu J j i = S10Pilot.momentum rho mu J j i +
      (if j = 0 then rho * (sigma - 1) * J 0 e * (if e = i then 1 else 0) else 0) := by
  unfold momentum S10Pilot.momentum
  rw [(lagrangian_variation e sigma rho mu J (basisJet j i)).deriv,
    (S10Pilot.lagrangian_variation rho mu J (basisJet j i)).deriv]
  simp only [variationDensity, basisJet_apply]
  by_cases h : j = 0 <;> simp [h, eq_comm]

theorem momentum_smul (e : Fin D) (sigma rho mu s : ℝ) (J : Jet D)
    (j : Fin (D + 1)) (i : Fin D) :
    momentum e sigma rho mu (s • J) j i = s * momentum e sigma rho mu J j i := by
  rw [momentum_eq, momentum_eq, S10Pilot.momentum_smul]
  by_cases h : j = 0
  · by_cases hi : e = i
    · simp [h, hi]
      ring
    · simp [h, hi]
  · simp [h]

def eulerLagrange (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (u : Point D → Vec D) (x : Point D) : Vec D :=
  fun i => -∑ j : Fin (D + 1), coordDeriv j (fun y => momentum e sigma rho mu (fieldJet u y) j i) x

theorem momentum_contraction (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (i : Fin D) :
    (∑ j : Fin (D + 1), waveCovector omega k j *
      momentum e sigma rho mu (modeJet (-omega) k a) j i) =
      modalOperator e sigma rho mu omega k a i := by
  simp only [momentum_eq, mul_add, Finset.sum_add_distrib]
  rw [S10Pilot.momentum_contraction]
  simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  by_cases h : i = e
  · subst i
    simp [modalOperator, unit, waveCovector, modeJet]
    ring
  · simp [modalOperator, unit, h, Ne.symm h]

/-- The independently differentiated position-space expression yields the same matrix. -/
theorem eulerLagrange_planeWave (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    eulerLagrange e sigma rho mu (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator e sigma rho mu omega k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul e sigma]
  have hrearrange : ∀ j : Fin (D + 1),
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum e sigma rho mu (modeJet (-omega) k a) j i) =
      (fun y => -momentum e sigma rho mu (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [hrearrange, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction e sigma]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- A global plane wave solves the local EL expression iff its amplitude is stationary. -/
theorem planeWave_solves_iff (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) :
    (∀ x, eulerLagrange e sigma rho mu (planeWave omega k a) x = 0) ↔
      ModalStationary e sigma rho mu omega k a := by
  rw [modal_stationary_iff e sigma]
  constructor
  · intro h
    simpa [eulerLagrange_planeWave e sigma, phase, dot] using h 0
  · intro h x
    simp [eulerLagrange_planeWave e sigma, h]

end
end S10Anisotropic
