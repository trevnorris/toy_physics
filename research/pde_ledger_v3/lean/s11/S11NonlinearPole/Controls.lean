import S11NonlinearPole.Modal
import S11NonlinearPole.Jordan

/-! NP4: admissible modal data and noncommuting/response witnesses. -/
namespace S11NonlinearPole
noncomputable section
open Matrix
open scoped Matrix.Norms.Elementwise

def scaledPair : PairingData ℂ ℂ ℂ ℂ where
  derivative := (2 : ℂ) • LinearMap.id
  right := LinearMap.id
  left := LinearMap.id
  pairing := LinearEquiv.smulOfNeZero ℂ ℂ 2 (by norm_num)
  pairing_eq := by intro k; rfl

theorem scaled_residue : scaledPair.residueCandidate 1 = 1/2 := by
  norm_num [PairingData.residueCandidate, scaledPair, LinearEquiv.smulOfNeZero,
    LinearEquiv.smulOfUnit, Units.smul_def]

theorem scaled_projection : scaledPair.fieldProjection 1 = 1 := by
  exact scaledPair.fieldProjection_fixes 1

def fullPair : PairingData ℂ (Fin 2 → ℂ) (Fin 2 → ℂ) (Fin 2 → ℂ) where
  derivative := LinearMap.id
  right := LinearMap.id
  left := LinearMap.id
  pairing := LinearEquiv.refl ℂ _
  pairing_eq := by intro k; rfl

theorem full_projection_rank : Module.finrank ℂ (LinearMap.range fullPair.fieldProjection) = 2 := by
  rw [fullPair.finrank_fieldProjection]
  simp

theorem zero_pairing_not_invertible : ¬ Function.Injective (fun x : ℂ => (0 : ℂ) * x) := by
  intro h
  have he := h (a₁ := 0) (a₂ := 1) (by simp)
  norm_num at he

/-- Reversing the coefficient/derivative order changes a measured entry. -/
theorem ordered_log_entry :
    moment 1 (fun z => doublePrincipal jordanN 0 z * (jordanN.transpose + z • 0)) 0 0 = 1 := by
  rw [double_log_moment _ _ _ _ 1 (by norm_num)]
  simp only [Matrix.add_apply, Matrix.mul_apply]
  norm_num [jordanN, Fin.sum_univ_two]

theorem reversed_log_entry : (jordanN.transpose * jordanN) 0 0 = 0 := by
  norm_num [jordanN, Matrix.mul_apply, Fin.sum_univ_two]

/-- Freezing both maps at zero misses the actual response residue 11. -/
theorem frozen_response_residue : (2 : ℂ) * moment 1 squareInverse * 1 = 0 := by
  rw [square_residue 1 (by norm_num)]
  ring

theorem simple_residue : moment 1 (fun z : ℂ => z⁻¹) = 1 := by
  simpa using moment_zpow 1 (by norm_num) (-1)

end
end S11NonlinearPole
