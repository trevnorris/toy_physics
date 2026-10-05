import S10Anisotropic.PlaneWave
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-! S10's phase average over 0..2π, with the normalization kept explicit.
No time period involving division by omega is used. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
variable {D : ℕ}

def phaseDensity (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (phi : ℝ) : ℝ :=
  lagrangian e sigma rho mu ((-Real.sin phi) • modeJet (-omega) k a)

theorem phaseDensity_eq (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (phi : ℝ) :
    phaseDensity e sigma rho mu omega k a phi = Real.sin phi ^ 2 * modalAction e sigma rho mu omega k a := by
  unfold phaseDensity
  rw [lagrangian_smul e sigma]
  change (-Real.sin phi) ^ 2 * modalAction e sigma rho mu (-omega) k a = _
  simp [modalAction_eq e sigma]

theorem phaseDensity_field (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) (x : Point D) :
    lagrangian e sigma rho mu (fieldJet (planeWave omega k a) x) =
      phaseDensity e sigma rho mu omega k a (phase (waveCovector omega k) x) := by
  rw [planeWave_lagrangian e sigma, phaseDensity_eq e sigma]

def phaseAverage (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) : ℝ :=
  (1 / (2 * Real.pi)) * ∫ phi in (0 : ℝ)..(2 * Real.pi), phaseDensity e sigma rho mu omega k a phi

theorem phaseAverage_eq (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) :
    phaseAverage e sigma rho mu omega k a = (1 / 2 : ℝ) * modalAction e sigma rho mu omega k a := by
  unfold phaseAverage
  simp_rw [phaseDensity_eq e sigma]
  rw [intervalIntegral.integral_mul_const, integral_sin_sq]
  simp only [Real.sin_zero, zero_mul, Real.sin_two_pi, sub_zero, zero_add]
  field_simp

theorem phaseAverage_variation (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => phaseAverage e sigma rho mu omega k (a + s • b))
      ((1 / 2 : ℝ) * dot (modalOperator e sigma rho mu omega k a) b) 0 := by
  simpa only [phaseAverage_eq e sigma] using (modalAction_variation e sigma rho mu omega k a b).const_mul (1 / 2)

theorem phaseAverage_stationary_iff (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) :
    (∀ b : Vec D, deriv (fun s : ℝ => phaseAverage e sigma rho mu omega k (a + s • b)) 0 = 0) ↔
      ModalStationary e sigma rho mu omega k a := by
  unfold ModalStationary
  simp_rw [(phaseAverage_variation e sigma rho mu omega k a _).deriv,
    (modalAction_variation e sigma rho mu omega k a _).deriv]
  simp

end
end S10Anisotropic
