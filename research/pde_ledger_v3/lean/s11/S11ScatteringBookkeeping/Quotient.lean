import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Ring

/-! P3: formal coefficients of the supplied rational model, with explicit domain. -/
namespace S11ScatteringBookkeeping
noncomputable section
variable {K : Type*} [Field K]

def c0 (n0 j0 : K) : K := n0 / j0
def c1 (n0 n1 j0 j1 : K) : K := (n1 - j1 * c0 n0 j0) / j0
def c2 (n0 n1 n2 j0 j1 j2 : K) : K :=
  (n2 - j1 * c1 n0 n1 j0 j1 - j2 * c0 n0 j0) / j0

theorem quotient_equations (n0 n1 n2 j0 j1 j2 : K) (h : j0 ≠ 0) :
    j0 * c0 n0 j0 = n0 ∧
    j0 * c1 n0 n1 j0 j1 + j1 * c0 n0 j0 = n1 ∧
    j0 * c2 n0 n1 n2 j0 j1 j2 + j1 * c1 n0 n1 j0 j1 + j2 * c0 n0 j0 = n2 := by
  constructor
  · exact mul_div_cancel₀ n0 h
  constructor
  · change j0 * ((n1 - j1 * c0 n0 j0) / j0) + j1 * c0 n0 j0 = n1
    rw [mul_div_cancel₀ _ h]
    ring
  · change j0 * ((n2 - j1 * c1 n0 n1 j0 j1 - j2 * c0 n0 j0) / j0) +
        j1 * c1 n0 n1 j0 j1 + j2 * c0 n0 j0 = n2
    rw [mul_div_cancel₀ _ h]
    ring

/-- Uniqueness uses the nonzero leading denominator; no analytic neighborhood is inferred. -/
theorem quotient_unique (n0 n1 n2 j0 j1 j2 u0 u1 u2 : K) (h : j0 ≠ 0)
    (h0 : j0 * u0 = n0) (h1 : j0 * u1 + j1 * u0 = n1)
    (h2 : j0 * u2 + j1 * u1 + j2 * u0 = n2) :
    u0 = c0 n0 j0 ∧ u1 = c1 n0 n1 j0 j1 ∧ u2 = c2 n0 n1 n2 j0 j1 j2 := by
  have e := quotient_equations n0 n1 n2 j0 j1 j2 h
  have e0 : u0 = c0 n0 j0 := mul_left_cancel₀ h (h0.trans e.1.symm)
  have e1 : u1 = c1 n0 n1 j0 j1 := by
    apply mul_left_cancel₀ h
    exact add_right_cancel (by simpa only [e0] using h1.trans e.2.1.symm)
  have e2 : u2 = c2 n0 n1 n2 j0 j1 j2 := by
    apply mul_left_cancel₀ h
    exact add_right_cancel (add_right_cancel
      (by simpa only [e0, e1] using h2.trans e.2.2.symm))
  exact ⟨e0, e1, e2⟩

theorem quotient_residual (n0 n1 n2 j0 j1 j2 t : K) (h : j0 ≠ 0) :
    (j0 + t * j1 + t ^ 2 * j2) *
      (c0 n0 j0 + t * c1 n0 n1 j0 j1 + t ^ 2 * c2 n0 n1 n2 j0 j1 j2) -
      (n0 + t * n1 + t ^ 2 * n2) =
    t ^ 3 * (j1 * c2 n0 n1 n2 j0 j1 j2 + j2 * c1 n0 n1 j0 j1) +
      t ^ 4 * (j2 * c2 n0 n1 n2 j0 j1 j2) := by
  have e := quotient_equations n0 n1 n2 j0 j1 j2 h
  calc
    _ = (j0 * c0 n0 j0 - n0) +
      t * (j0 * c1 n0 n1 j0 j1 + j1 * c0 n0 j0 - n1) +
      t ^ 2 * (j0 * c2 n0 n1 n2 j0 j1 j2 + j1 * c1 n0 n1 j0 j1 +
        j2 * c0 n0 j0 - n2) +
      t ^ 3 * (j1 * c2 n0 n1 n2 j0 j1 j2 + j2 * c1 n0 n1 j0 j1) +
      t ^ 4 * (j2 * c2 n0 n1 n2 j0 j1 j2) := by ring
    _ = _ := by rw [e.1, e.2.1, e.2.2]; ring

theorem leading_denominator_cases (j0 : K) : j0 = 0 ∨ j0 ≠ 0 := eq_or_ne j0 0

theorem zero_leading_obstruction (n0 : K) (h : n0 ≠ 0) (u0 : K) :
    (0 : K) * u0 ≠ n0 := by simpa only [zero_mul, ne_eq, eq_comm] using h

theorem scaled_denominator_nonzero (eps den : K) (he : eps ≠ 0) (hd : den ≠ 0) :
    eps ^ 2 * den ≠ 0 := mul_ne_zero (pow_ne_zero 2 he) hd

theorem epsilon_cancels (eps num den : K) (he : eps ≠ 0) :
    (eps ^ 2 * num) / (eps ^ 2 * den) = num / den :=
  mul_div_mul_left num den (pow_ne_zero 2 he)

end
end S11ScatteringBookkeeping
