import S11ScatteringFlux.Current

/-! P1: a supplied four-term rectangle and its specified real path. -/
namespace S11ScatteringBookkeeping
noncomputable section
open Matrix S11ScatteringFlux
variable {n : Type*}

def amplitude (a0 a1 a2 : n → ℂ) (t : ℝ) : n → ℂ :=
  a0 + (t : ℂ) • a1 + ((t ^ 2 : ℝ) : ℂ) • a2

def rectangle (a00 a10 a01 a11 : n → ℂ) (eta sigma : ℝ) : n → ℂ :=
  a00 + (eta : ℂ) • a10 + (sigma : ℂ) • a01 + ((eta * sigma : ℝ) : ℂ) • a11

def delta (a10 a01 a11 : n → ℂ) (eta sigma : ℝ) : n → ℂ :=
  (eta : ℂ) • a10 + ((sigma : ℂ) • a01 + ((eta * sigma : ℝ) : ℂ) • a11)

theorem rectangle_decomposition (a00 a10 a01 a11 : n → ℂ) (eta sigma : ℝ) :
    rectangle a00 a10 a01 a11 eta sigma = a00 + delta a10 a01 a11 eta sigma := by
  simp only [rectangle, delta, add_assoc]

theorem rectangle_path (a00 a10 a01 a11 : n → ℂ) (t r : ℝ) :
    rectangle a00 a10 a01 a11 t (r * t) =
      amplitude a00 (a10 + (r : ℂ) • a01) ((r : ℂ) • a11) t := by
  ext i
  simp only [rectangle, amplitude, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    Complex.ofReal_mul, Complex.ofReal_pow]
  ring

theorem baseline_free (a10 a01 a11 : n → ℂ) (eta sigma : ℝ) :
    rectangle 0 a10 a01 a11 eta sigma = delta a10 a01 a11 eta sigma := by
  rw [rectangle_decomposition, zero_add]

/-- The map to one real path can have a kernel even on the retained rectangle. -/
theorem path_does_not_identify_rectangle (v : n → ℂ) (t r : ℝ) :
    rectangle 0 (-((r : ℂ) • v)) v 0 t (r * t) = 0 := by
  ext i
  simp only [rectangle, Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply,
    smul_eq_mul, Complex.ofReal_mul]
  ring

end
end S11ScatteringBookkeeping
