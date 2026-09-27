import S11ScatteringBookkeeping.Observables
import S11ScatteringBookkeeping.Quotient
import Mathlib.Tactic.NormNum
import Mathlib.Algebra.BigOperators.Fin

namespace S11ScatteringBookkeeping
noncomputable section
open Matrix S11ScatteringFlux

def scalarA (z : ℂ) : Fin 1 → ℂ := fun _ => z
def scalarB (z : ℂ) : Matrix (Fin 1) (Fin 1) ℂ := fun _ _ => z

theorem path_witness : rectangle (scalarA 1) (scalarA 2) (scalarA 3) (scalarA 4) 1 2 0 = 17 := by
  norm_num [rectangle, scalarA, smul_eq_mul]

theorem delta_witness : delta (scalarA 2) (scalarA 3) (scalarA 4) 1 2 0 = 16 := by
  norm_num [delta, scalarA, smul_eq_mul]

theorem second_witness : q2 (scalarB 1) (scalarB 4) (scalarB 5)
    (scalarA 1) (scalarA 2) (scalarA 3) = 31 := by
  norm_num [q2, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem first_witness : q1 (scalarB 1) (scalarB 4) (scalarB 5)
    (scalarA 1) (scalarA 2) (scalarA 3) = 8 := by
  norm_num [q1, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem higher_witness : q3 (scalarB 1) (scalarB 4) (scalarB 5)
    (scalarA 1) (scalarA 2) (scalarA 3) = 72 := by
  norm_num [q3, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem induced_witness : q2 (scalarB 1) (scalarB 4) (scalarB 5)
    0 (scalarA 2) (scalarA 3) = 4 := by
  rw [induced_coefficient]
  norm_num [pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem subtracted_witness :
    flux (scalarB 1) (scalarA 1 + scalarA 1) - flux (scalarB 1) (scalarA 1) = 3 := by
  norm_num [flux, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem parent_second_witness :
    q2 (scalarB 1) 0 0 (scalarA 1) 0 (scalarA 3) = 6 := by
  norm_num [q2, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem epsilon_witness : flux (scalarB 1) ((2 : ℂ) • scalarA 1) = 4 := by
  calc
    _ = (2 : ℝ) ^ 2 * flux (scalarB 1) (scalarA 1) := by
      simpa only [Complex.ofReal_ofNat] using
        flux_epsilon_squared (scalarB 1) (scalarA 1) 2
    _ = 4 := by norm_num [flux, pair, scalarB, scalarA, Matrix.mulVec, dotProduct]

theorem zero_model_two_solutions :
    (0 : ℝ) * 0 = 0 ∧ (0 : ℝ) * 1 = 0 ∧ (0 : ℝ) ≠ 1 := by norm_num

end
end S11ScatteringBookkeeping
