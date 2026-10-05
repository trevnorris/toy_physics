import S11ScatteringBookkeeping.Quadratic

/-! P4: total, baseline subtraction and coherent amplitude subtraction differ. -/
namespace S11ScatteringBookkeeping
noncomputable section
open Matrix S11ScatteringFlux
variable {n : Type*} [Fintype n]

theorem subtracted_total (B : Matrix n n ℂ) (baseline induced : n → ℂ) :
    flux B (baseline + induced) - flux B baseline =
      flux B induced + (pair B baseline induced + pair B induced baseline).re := by
  rw [flux_add]
  ring

theorem subtracted_eq_induced_iff (B : Matrix n n ℂ) (baseline induced : n → ℂ) :
    flux B (baseline + induced) - flux B baseline = flux B induced ↔
      (pair B baseline induced + pair B induced baseline).re = 0 := by
  rw [subtracted_total]
  constructor <;> intro h <;> linarith

theorem baseline_zero_total (B : Matrix n n ℂ) (induced : n → ℂ) :
    flux B (0 + induced) = flux B induced := by rw [zero_add]

theorem induced_coefficient (B0 B1 B2 : Matrix n n ℂ) (a1 a2 : n → ℂ) :
    q2 B0 B1 B2 0 a1 a2 = pair B0 a1 a1 := by
  have hz : star (0 : n → ℂ) = 0 := funext fun _ => map_zero (starRingEnd ℂ)
  simp [q2, pair, hz]

theorem induced_low_coefficients (B0 B1 B2 : Matrix n n ℂ) (a1 a2 : n → ℂ) :
    q0 B0 B1 B2 0 a1 a2 = 0 ∧ q1 B0 B1 B2 0 a1 a2 = 0 := by
  have hz : star (0 : n → ℂ) = 0 := funext fun _ => map_zero (starRingEnd ℂ)
  simp [q0, q1, pair, hz]

end
end S11ScatteringBookkeeping
