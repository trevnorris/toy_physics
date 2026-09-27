import S11ScatteringFlux.Balance
import Mathlib.Tactic.NormNum
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Algebra.BigOperators.Fin

/-! F4: exact nonzero and exceptional-domain witnesses. -/
namespace S11ScatteringFlux
noncomputable section
open Matrix

def coherent : Matrix (Fin 2) (Fin 2) ℂ := !![1, 1; 1, 1]
def diagonalOnly : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, 1]
def both : Fin 2 → ℂ := ![1, 1]
def oneCurrent : Matrix (Fin 1) (Fin 1) ℂ := fun _ _ => 1
def phaseAmplitude : Fin 1 → ℂ := fun _ => Complex.I
def realAmplitude : Fin 1 → ℂ := fun _ => 1
def doubleBasis : Matrix (Fin 1) (Fin 1) ℂ := fun _ _ => 2

theorem coherent_flux : flux coherent both = 4 := by
  norm_num [flux, pair, coherent, both, Matrix.mulVec, dotProduct, Fin.sum_univ_two]

theorem diagonal_flux : flux diagonalOnly both = 2 := by
  norm_num [flux, pair, diagonalOnly, both, Matrix.mulVec, dotProduct, Fin.sum_univ_two]

theorem phase_flux : flux oneCurrent phaseAmplitude = 1 := by
  norm_num [flux, pair, oneCurrent, phaseAmplitude, Matrix.mulVec, dotProduct]

theorem basis_flux : flux (pullback oneCurrent doubleBasis) realAmplitude = 4 := by
  norm_num [flux, pair, pullback, oneCurrent, doubleBasis, realAmplitude,
    Matrix.mul_apply, Matrix.conjTranspose_apply, Matrix.mulVec, dotProduct]

/-- An indefinite supplied current may have zero flux at a nonzero amplitude. -/
theorem null_flux_nonzero_amplitude :
    flux (!![1, 0; 0, -1] : Matrix (Fin 2) (Fin 2) ℂ) both = 0 ∧ both ≠ 0 := by
  constructor
  · norm_num [flux, pair, both, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  · intro h
    have := congrFun h 0
    norm_num [both] at this

theorem nonhermitian_imaginary_witness :
    (pair (fun _ _ : Fin 1 => Complex.I) realAmplitude realAmplitude).im = 1 := by
  norm_num [pair, realAmplitude, Matrix.mulVec, dotProduct]

theorem empty_flux (J : Matrix (Fin 0) (Fin 0) ℂ) (x : Fin 0 → ℂ) : flux J x = 0 := by
  simp [flux, pair, dotProduct]

theorem signed_incident_witness : incident .left (-3) = -3 := by norm_num [incident_left]
theorem oriented_outgoing_witness : outgoing 3 1 = -2 := by norm_num [outgoing_eq]
theorem zero_denominator_witness : fraction 7 0 = none := by norm_num [fraction]
theorem positive_fraction_witness : fraction 2 4 = some (1 / 2) := by norm_num [fraction]
theorem negative_fraction_witness : fraction 2 (-4) = some (-1 / 2) := by norm_num [fraction]
theorem fraction_above_one_witness : fraction 4 2 = some 2 := by norm_num [fraction]
theorem nonconservative_witness : (1 : ℝ) / 4 + 1 / 4 = 1 - 2 / 4 := by norm_num

end
end S11ScatteringFlux
