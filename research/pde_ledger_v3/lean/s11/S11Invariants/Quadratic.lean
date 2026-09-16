import Mathlib.LinearAlgebra.QuadraticForm.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.Tactic

/-! Actual quadratic forms on all real 2x2 matrices. The ten-coefficient
presentation below is proved exhaustive using the associated bilinear form. -/

namespace S11Invariants
noncomputable section

abbrev Mat := Matrix (Fin 2) (Fin 2) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 4 → ℝ
abbrev Coefficients := Fin 10 → ℝ

def coordinates (G : Mat) : Coordinates :=
  ![G 0 0 + G 1 1, G 0 1 - G 1 0, G 0 0 - G 1 1, G 0 1 + G 1 0]

def decode (v : Coordinates) : Mat :=
  ![![(v 0 + v 2) / 2, (v 1 + v 3) / 2],
    ![(v 3 - v 1) / 2, (v 0 - v 2) / 2]]

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode] <;> ring

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 4) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates] <;> ring
  map_smul' r G := by fin_cases i <;> simp [coordinates] <;> ring

def frame : Fin 4 → Mat :=
  ![decode ![1,0,0,0], decode ![0,1,0,0], decode ![0,0,1,0], decode ![0,0,0,1]]

theorem frame_expansion (G : Mat) :
    G = coordinates G 0 • frame 0 + coordinates G 1 • frame 1 +
        coordinates G 2 • frame 2 + coordinates G 3 • frame 3 := by
  ext i j
  fin_cases i <;> fin_cases j
  · change G 0 0 = (G 0 0 + G 1 1) * ((1+0)/2) + (G 0 1-G 1 0) * ((0+0)/2) +
      (G 0 0-G 1 1) * ((0+1)/2) + (G 0 1+G 1 0) * ((0+0)/2)
    ring
  · change G 0 1 = (G 0 0 + G 1 1) * ((0+0)/2) + (G 0 1-G 1 0) * ((1+0)/2) +
      (G 0 0-G 1 1) * ((0+0)/2) + (G 0 1+G 1 0) * ((0+1)/2)
    ring
  · change G 1 0 = (G 0 0 + G 1 1) * ((0-0)/2) + (G 0 1-G 1 0) * ((0-1)/2) +
      (G 0 0-G 1 1) * ((0-0)/2) + (G 0 1+G 1 0) * ((1-0)/2)
    ring
  · change G 1 1 = (G 0 0 + G 1 1) * ((1-0)/2) + (G 0 1-G 1 0) * ((0-0)/2) +
      (G 0 0-G 1 1) * ((0-1)/2) + (G 0 1+G 1 0) * ((0-0)/2)
    ring

def polynomial (c : Coefficients) (v : Coordinates) : ℝ :=
  c 0 * v 0 ^ 2 + c 1 * v 0 * v 1 + c 2 * v 0 * v 2 + c 3 * v 0 * v 3 +
  c 4 * v 1 ^ 2 + c 5 * v 1 * v 2 + c 6 * v 1 * v 3 +
  c 7 * v 2 ^ 2 + c 8 * v 2 * v 3 + c 9 * v 3 ^ 2

def monomial (i j : Fin 4) : Quad :=
  QuadraticMap.linMulLin (coordinateMap i) (coordinateMap j)

def polynomialForm (c : Coefficients) : Quad :=
  c 0 • monomial 0 0 + c 1 • monomial 0 1 + c 2 • monomial 0 2 + c 3 • monomial 0 3 +
  c 4 • monomial 1 1 + c 5 • monomial 1 2 + c 6 • monomial 1 3 +
  c 7 • monomial 2 2 + c 8 • monomial 2 3 + c 9 • monomial 3 3

theorem polynomialForm_apply (c : Coefficients) (G : Mat) :
    polynomialForm c G = polynomial c (coordinates G) := by
  simp [polynomialForm, monomial, coordinateMap, polynomial, pow_two]
  ring

theorem quadratic_representation (Q : Quad) :
    ∃ c : Coefficients, ∀ G, Q G = polynomial c (coordinates G) := by
  let b := Q.associated
  let c : Coefficients :=
    ![b (frame 0) (frame 0), b (frame 0) (frame 1) + b (frame 1) (frame 0),
      b (frame 0) (frame 2) + b (frame 2) (frame 0),
      b (frame 0) (frame 3) + b (frame 3) (frame 0), b (frame 1) (frame 1),
      b (frame 1) (frame 2) + b (frame 2) (frame 1),
      b (frame 1) (frame 3) + b (frame 3) (frame 1), b (frame 2) (frame 2),
      b (frame 2) (frame 3) + b (frame 3) (frame 2), b (frame 3) (frame 3)]
  refine ⟨c, fun G => ?_⟩
  have hb : Q G = b G G := (Q.associated_eq_self_apply ℝ G).symm
  rw [hb, congrArg (fun Z => b Z Z) (frame_expansion G)]
  simp only [polynomial, c, Matrix.cons_val]
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
  ring

theorem polynomial_injective {c d : Coefficients}
    (h : ∀ v, polynomial c v = polynomial d v) : c = d := by
  have h0 := h ![1,0,0,0]
  have h1 := h ![1,1,0,0]
  have h2 := h ![1,0,1,0]
  have h3 := h ![1,0,0,1]
  have h4 := h ![0,1,0,0]
  have h5 := h ![0,1,1,0]
  have h6 := h ![0,1,0,1]
  have h7 := h ![0,0,1,0]
  have h8 := h ![0,0,1,1]
  have h9 := h ![0,0,0,1]
  norm_num [polynomial, Matrix.cons_val_two, Matrix.cons_val_three] at h0 h1 h2 h3 h4 h5 h6 h7 h8 h9
  ext i
  fin_cases i
  · change c 0 = d 0
    linarith
  · change c 1 = d 1
    linarith
  · change c 2 = d 2
    linarith
  · change c 3 = d 3
    linarith
  · change c 4 = d 4
    linarith
  · change c 5 = d 5
    linarith
  · change c 6 = d 6
    linarith
  · change c 7 = d 7
    linarith
  · change c 8 = d 8
    linarith
  · change c 9 = d 9
    linarith

theorem polynomialForm_injective : Function.Injective polynomialForm := by
  intro c d h
  apply polynomial_injective
  intro v
  have hv := congrArg (fun Q : Quad => Q (decode v)) h
  simpa only [polynomialForm_apply, coordinates_decode] using hv

theorem polynomialForm_surjective : Function.Surjective polynomialForm := by
  intro Q
  obtain ⟨c, hc⟩ := quadratic_representation Q
  refine ⟨c, ?_⟩
  ext G
  rw [polynomialForm_apply, hc]

/-- Images of basis polynomials stored as rows act on coefficient columns by
the transpose. This identity is independent of dimension or a CAS row layout. -/
theorem coefficient_action {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ)
    (c m : Fin n → ℝ) :
    (∑ i, c i * (∑ j, A i j * m j)) =
      ∑ j, (A.transpose.mulVec c) j * m j := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
  simp_rw [Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  congr 1
  funext j
  apply Finset.sum_congr rfl
  intro i _
  ring

end
end S11Invariants
