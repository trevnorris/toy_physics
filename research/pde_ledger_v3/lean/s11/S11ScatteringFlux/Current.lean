import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-! F1–F2: the entire supplied finite complex current form, with no diagonal
or positivity assumption. Physical validity of that form is an application premise. -/
namespace S11ScatteringFlux
noncomputable section
open Matrix
variable {n m k : Type*} [Fintype n] [Fintype m] [Fintype k]

def pair (J : Matrix n n ℂ) (x y : n → ℂ) : ℂ := star x ⬝ᵥ (J *ᵥ y)
def flux (J : Matrix n n ℂ) (x : n → ℂ) : ℝ := (pair J x x).re
def pullback (J : Matrix n n ℂ) (C : Matrix n m ℂ) : Matrix m m ℂ := Cᴴ * J * C

theorem pair_coordinates (J : Matrix n n ℂ) (x y : n → ℂ) :
    pair J x y = ∑ j, ∑ i, star (x i) * J i j * y j :=
  dot_mulVec_eq_sum_sum _ _ _

theorem pair_add_left (J : Matrix n n ℂ) (x y z : n → ℂ) :
    pair J (x + y) z = pair J x z + pair J y z := by
  have h : star (x + y) = star x + star y := funext fun i => star_add (x i) (y i)
  simp only [pair, h, add_dotProduct]

theorem pair_add_right (J : Matrix n n ℂ) (x y z : n → ℂ) :
    pair J x (y + z) = pair J x y + pair J x z := by
  simp only [pair, mulVec_add, dotProduct_add]

theorem pair_add_matrix (J K : Matrix n n ℂ) (x y : n → ℂ) :
    pair (J + K) x y = pair J x y + pair K x y := by
  simp only [pair, add_mulVec, dotProduct_add]

theorem pair_conjTranspose (J : Matrix n n ℂ) (x y : n → ℂ) :
    pair Jᴴ x y = star (pair J y x) := by
  unfold pair
  rw [dotProduct_mulVec, ← star_mulVec]
  exact star_dotProduct _ _

theorem pair_hermitian (J : Matrix n n ℂ) (hJ : Jᴴ = J) (x y : n → ℂ) :
    pair J x y = star (pair J y x) := by
  rw [← pair_conjTranspose, hJ]

theorem hermitian_flux_real (J : Matrix n n ℂ) (hJ : Jᴴ = J) (x : n → ℂ) :
    (pair J x x).im = 0 := by
  have h := congrArg Complex.im (pair_hermitian J hJ x x)
  simp only [Complex.star_def, Complex.conj_im] at h
  linarith

theorem flux_add (J : Matrix n n ℂ) (x y : n → ℂ) :
    flux J (x + y) = flux J x + flux J y + (pair J x y + pair J y x).re := by
  simp only [flux, pair_add_left, pair_add_right, Complex.add_re]
  ring

theorem flux_add_iff_cross_zero (J : Matrix n n ℂ) (x y : n → ℂ) :
    flux J (x + y) = flux J x + flux J y ↔ (pair J x y + pair J y x).re = 0 := by
  rw [flux_add]
  constructor <;> intro h <;> linarith

theorem pair_pullback (J : Matrix n n ℂ) (C : Matrix n m ℂ) (x y : m → ℂ) :
    pair (pullback J C) x y = pair J (C *ᵥ x) (C *ᵥ y) := by
  simp only [pair, pullback, star_mulVec, dotProduct_mulVec, vecMul_vecMul]

theorem flux_pullback (J : Matrix n n ℂ) (C : Matrix n m ℂ) (x : m → ℂ) :
    flux (pullback J C) x = flux J (C *ᵥ x) := by
  exact congrArg Complex.re (pair_pullback J C x x)

omit [Fintype k] in
theorem pullback_comp (J : Matrix n n ℂ) (C : Matrix n m ℂ) (D : Matrix m k ℂ) :
    pullback (pullback J C) D = pullback J (C * D) := by
  simp only [pullback, conjTranspose_mul, Matrix.mul_assoc]

/-- C_out converts new outgoing coordinates to old ones; D_out is its supplied
right inverse. Input coordinates change through C_in. No channel selection is proved. -/
theorem scattering_flux_covariant [DecidableEq n]
    (J : Matrix n n ℂ) (S : Matrix n m ℂ)
    (Cout Dout : Matrix n n ℂ) (Cin : Matrix m k ℂ)
    (hInv : Cout * Dout = 1) (x : k → ℂ) :
    flux (pullback J Cout) ((Dout * S * Cin) *ᵥ x) = flux J (S *ᵥ (Cin *ᵥ x)) := by
  rw [flux_pullback]
  congr 1
  simp only [mulVec_mulVec, ← Matrix.mul_assoc, hInv, Matrix.one_mul]

end
end S11ScatteringFlux
