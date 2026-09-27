import S11ScatteringBookkeeping.Path

/-! P2: exact contraction of supplied degree-two amplitude and current polynomials.
All higher powers remain in an explicit remainder. No parent-theory Taylor claim. -/
namespace S11ScatteringBookkeeping
noncomputable section
open Matrix S11ScatteringFlux
variable {n : Type*} [Fintype n]

def current (B0 B1 B2 : Matrix n n ℂ) (t : ℝ) : Matrix n n ℂ :=
  B0 + (t : ℂ) • B1 + ((t ^ 2 : ℝ) : ℂ) • B2

theorem pair_smul_left (B : Matrix n n ℂ) (x y : n → ℂ) (s : ℂ) :
    pair B (s • x) y = star s * pair B x y := by
  have hs : star (s • x) = star s • star x :=
    funext fun i => map_mul (starRingEnd ℂ) s (x i)
  simp only [pair, hs, smul_dotProduct, smul_eq_mul]

theorem pair_smul_right (B : Matrix n n ℂ) (x y : n → ℂ) (s : ℂ) :
    pair B x (s • y) = s * pair B x y := by
  simp only [pair, mulVec_smul, dotProduct_smul, smul_eq_mul]

theorem pair_smul_matrix (B : Matrix n n ℂ) (x y : n → ℂ) (s : ℂ) :
    pair (s • B) x y = s * pair B x y := by
  simp only [pair, smul_mulVec, dotProduct_smul, smul_eq_mul]


def q0 (B0 _B1 _B2 : Matrix n n ℂ) (a0 _a1 _a2 : n → ℂ) : ℂ :=
  pair B0 a0 a0

def q1 (B0 B1 _B2 : Matrix n n ℂ) (a0 a1 _a2 : n → ℂ) : ℂ :=
  pair B0 a0 a1 +
    pair B1 a0 a0 +
    pair B0 a1 a0

def q2 (B0 B1 B2 : Matrix n n ℂ) (a0 a1 a2 : n → ℂ) : ℂ :=
  pair B0 a0 a2 +
    pair B1 a0 a1 +
    pair B2 a0 a0 +
    pair B0 a1 a1 +
    pair B1 a1 a0 +
    pair B0 a2 a0

def q3 (B0 B1 B2 : Matrix n n ℂ) (a0 a1 a2 : n → ℂ) : ℂ :=
  pair B1 a0 a2 +
    pair B2 a0 a1 +
    pair B0 a1 a2 +
    pair B1 a1 a1 +
    pair B2 a1 a0 +
    pair B0 a2 a1 +
    pair B1 a2 a0

def q4 (B0 B1 B2 : Matrix n n ℂ) (a0 a1 a2 : n → ℂ) : ℂ :=
  pair B2 a0 a2 +
    pair B1 a1 a2 +
    pair B2 a1 a1 +
    pair B0 a2 a2 +
    pair B1 a2 a1 +
    pair B2 a2 a0

def q5 (_B0 B1 B2 : Matrix n n ℂ) (_a0 a1 a2 : n → ℂ) : ℂ :=
  pair B2 a1 a2 +
    pair B1 a2 a2 +
    pair B2 a2 a1

def q6 (_B0 _B1 B2 : Matrix n n ℂ) (_a0 _a1 a2 : n → ℂ) : ℂ :=
  pair B2 a2 a2

theorem quadratic_exact (B0 B1 B2 : Matrix n n ℂ) (a0 a1 a2 : n → ℂ) (t : ℝ) :
    pair (current B0 B1 B2 t) (amplitude a0 a1 a2 t) (amplitude a0 a1 a2 t) =
      q0 B0 B1 B2 a0 a1 a2 +
      ((t ^ 1 : ℝ) : ℂ) * q1 B0 B1 B2 a0 a1 a2 +
      ((t ^ 2 : ℝ) : ℂ) * q2 B0 B1 B2 a0 a1 a2 +
      ((t ^ 3 : ℝ) : ℂ) * q3 B0 B1 B2 a0 a1 a2 +
      ((t ^ 4 : ℝ) : ℂ) * q4 B0 B1 B2 a0 a1 a2 +
      ((t ^ 5 : ℝ) : ℂ) * q5 B0 B1 B2 a0 a1 a2 +
      ((t ^ 6 : ℝ) : ℂ) * q6 B0 B1 B2 a0 a1 a2 := by
  simp only [current, amplitude, pair_add_matrix, pair_add_left, pair_add_right,
    pair_smul_left, pair_smul_right, pair_smul_matrix,
    q0, q1, q2, q3, q4, q5, q6, Complex.star_def, map_pow, Complex.conj_ofReal,
    Complex.ofReal_pow]
  ring

def remainder (B0 B1 B2 : Matrix n n ℂ) (a0 a1 a2 : n → ℂ) (t : ℝ) : ℂ :=
  q3 B0 B1 B2 a0 a1 a2 + (t : ℂ) * q4 B0 B1 B2 a0 a1 a2 +
    ((t ^ 2 : ℝ) : ℂ) * q5 B0 B1 B2 a0 a1 a2 + ((t ^ 3 : ℝ) : ℂ) * q6 B0 B1 B2 a0 a1 a2

theorem retained_with_remainder (B0 B1 B2 : Matrix n n ℂ)
    (a0 a1 a2 : n → ℂ) (t : ℝ) :
    pair (current B0 B1 B2 t) (amplitude a0 a1 a2 t) (amplitude a0 a1 a2 t) =
      q0 B0 B1 B2 a0 a1 a2 + (t : ℂ) * q1 B0 B1 B2 a0 a1 a2 +
        ((t ^ 2 : ℝ) : ℂ) * q2 B0 B1 B2 a0 a1 a2 + ((t ^ 3 : ℝ) : ℂ) * remainder B0 B1 B2 a0 a1 a2 t := by
  rw [quadratic_exact]
  simp only [remainder, Complex.ofReal_pow]
  ring

theorem real_flux_expansion (B0 B1 B2 : Matrix n n ℂ)
    (a0 a1 a2 : n → ℂ) (t : ℝ) :
    flux (current B0 B1 B2 t) (amplitude a0 a1 a2 t) =
      (q0 B0 B1 B2 a0 a1 a2).re + t * (q1 B0 B1 B2 a0 a1 a2).re +
        t ^ 2 * (q2 B0 B1 B2 a0 a1 a2).re +
          t ^ 3 * (remainder B0 B1 B2 a0 a1 a2 t).re := by
  unfold flux
  rw [retained_with_remainder]
  simp only [Complex.add_re, Complex.re_ofReal_mul]

theorem flux_epsilon_squared (B : Matrix n n ℂ) (a : n → ℂ) (eps : ℝ) :
    flux B ((eps : ℂ) • a) = eps ^ 2 * flux B a := by
  simp only [flux, pair_smul_left, pair_smul_right, Complex.star_def, Complex.conj_ofReal,
    Complex.re_ofReal_mul]
  ring

end
end S11ScatteringBookkeeping
