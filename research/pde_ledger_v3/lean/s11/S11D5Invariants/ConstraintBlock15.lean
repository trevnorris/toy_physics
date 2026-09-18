import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock15 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 38 + (1) * c 43 + (-1) * c 259 + (1) * c 304 = 0) ∧
    ((-1) * c 38 + (-1) * c 43 + (1) * c 259 + (-1) * c 304 = 0) ∧
    ((-1) * c 46 + (1) * c 47 + (-1) * c 319 + (1) * c 322 = 0) ∧
    ((-1) * c 46 + (-1) * c 47 + (1) * c 319 + (-1) * c 322 = 0) ∧
    ((-1) * c 49 + (-1) * c 51 + (1) * c 72 + (1) * c 73 = 0) ∧
    ((-1) * c 49 + (-1) * c 57 + (1) * c 72 + (1) * c 84 + (-1) * c 205 + (1) * c 270 = 0) ∧
    ((-1) * c 49 + (-1) * c 67 + (1) * c 72 + (1) * c 89 = 0) ∧
    ((-1) * c 100 + (1) * c 105 + (-1) * c 205 + (1) * c 270 = 0) ∧
    ((-1) * c 205 + (-1) * c 215 + (1) * c 270 + (1) * c 275 = 0) ∧
    ((-1) * c 234 + (-1) * c 236 + (-1) * c 259 + (1) * c 297 + (1) * c 298 + (1) * c 304 = 0) ∧
    ((-1) * c 234 + (-1) * c 244 + (1) * c 297 + (1) * c 302 + (-1) * c 319 + (1) * c 322 = 0) ∧
    ((-1) * c 234 + (-1) * c 246 + (1) * c 297 + (1) * c 303 = 0) ∧
    ((-1) * c 259 + (-1) * c 267 + (1) * c 304 + (1) * c 308 + (-1) * c 319 + (1) * c 322 = 0) ∧
    ((-1) * c 259 + (-1) * c 269 + (1) * c 304 + (1) * c 309 = 0) ∧
    ((-1) * c 319 + (-1) * c 321 + (1) * c 322 + (1) * c 323 = 0) ∧
    ((-1) * c 3 + (1) * c 4 + (-1) * c 72 + (1) * c 94 = 0) ∧
    ((-1) * c 15 + (1) * c 20 + (-1) * c 270 + (1) * c 310 = 0) ∧
    ((-1) * c 18 + (1) * c 24 + (-1) * c 297 + (1) * c 324 = 0) ∧
    ((-1) * c 37 + (1) * c 38 + (-1) * c 247 + (1) * c 259 = 0) ∧
    ((-1) * c 42 + (1) * c 48 + (-1) * c 297 + (1) * c 324 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
