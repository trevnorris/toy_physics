import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock10 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 270 + (-1) * c 279 + (1) * c 280 + (1) * c 288 = 0) ∧
    ((1) * c 270 + (-1) * c 272 + (-1) * c 280 + (-1) * c 281 = 0) ∧
    ((1) * c 270 + (-1) * c 273 + (-1) * c 280 + (-1) * c 282 = 0) ∧
    ((1) * c 270 + (-1) * c 274 + (-1) * c 280 + (-1) * c 283 = 0) ∧
    ((1) * c 270 + (-1) * c 277 + (-1) * c 280 + (-1) * c 286 = 0) ∧
    ((1) * c 270 + (-1) * c 278 + (-1) * c 280 + (-1) * c 287 = 0) ∧
    ((1) * c 270 + (-1) * c 279 + (-1) * c 280 + (-1) * c 288 = 0) ∧
    ((-1) * c 292 + (1) * c 293 + (-1) * c 310 + (1) * c 315 = 0) ∧
    ((-1) * c 292 + (-1) * c 293 + (1) * c 310 + (-1) * c 315 = 0) ∧
    ((-1) * c 299 + (1) * c 300 + (-1) * c 310 + (1) * c 315 = 0) ∧
    ((-1) * c 299 + (-1) * c 300 + (1) * c 310 + (-1) * c 315 = 0) ∧
    ((-1) * c 305 + (1) * c 306 + (-1) * c 310 + (1) * c 315 = 0) ∧
    ((-1) * c 305 + (-1) * c 306 + (1) * c 310 + (-1) * c 315 = 0) ∧
    ((-2) * c 311 = 0) ∧
    ((-1) * c 310 + (-1) * c 312 + (1) * c 315 + (1) * c 316 = 0) ∧
    ((-1) * c 310 + (-1) * c 313 + (1) * c 315 + (1) * c 317 = 0) ∧
    ((-1) * c 310 + (-1) * c 314 + (1) * c 315 + (1) * c 318 = 0) ∧
    ((1) * c 310 + (-1) * c 312 + (-1) * c 315 + (-1) * c 316 = 0) ∧
    ((1) * c 310 + (-1) * c 313 + (-1) * c 315 + (-1) * c 317 = 0) ∧
    ((1) * c 310 + (-1) * c 314 + (-1) * c 315 + (-1) * c 318 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
