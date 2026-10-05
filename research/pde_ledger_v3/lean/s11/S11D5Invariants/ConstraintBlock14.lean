import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock14 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 189 + (-1) * c 201 + (1) * c 259 + (1) * c 267 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 189 + (-1) * c 203 + (1) * c 259 + (1) * c 268 = 0) ∧
    ((-1) * c 189 + (-1) * c 204 + (1) * c 259 + (1) * c 269 = 0) ∧
    ((-1) * c 276 + (1) * c 277 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 280 + (-1) * c 282 + (1) * c 289 + (1) * c 290 = 0) ∧
    ((-1) * c 280 + (-1) * c 283 + (1) * c 289 + (1) * c 291 = 0) ∧
    ((-1) * c 280 + (-1) * c 285 + (1) * c 289 + (1) * c 294 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 280 + (-1) * c 287 + (1) * c 289 + (1) * c 295 = 0) ∧
    ((-1) * c 280 + (-1) * c 288 + (1) * c 289 + (1) * c 296 = 0) ∧
    ((-1) * c 300 + (1) * c 301 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 306 + (1) * c 307 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 315 + (-1) * c 317 + (1) * c 319 + (1) * c 320 = 0) ∧
    ((-1) * c 315 + (-1) * c 318 + (1) * c 319 + (1) * c 321 = 0) ∧
    ((-1) * c 2 + (1) * c 3 + (-1) * c 49 + (1) * c 72 = 0) ∧
    ((-1) * c 10 + (1) * c 15 + (-1) * c 205 + (1) * c 270 = 0) ∧
    ((-1) * c 12 + (1) * c 18 + (-1) * c 234 + (1) * c 297 = 0) ∧
    ((-1) * c 14 + (1) * c 19 + (-1) * c 259 + (1) * c 304 = 0) ∧
    ((-1) * c 22 + (1) * c 23 + (-1) * c 319 + (1) * c 322 = 0) ∧
    ((-1) * c 36 + (1) * c 42 + (-1) * c 234 + (1) * c 297 = 0) ∧
    ((-1) * c 37 + (-1) * c 41 + (-1) * c 247 + (1) * c 289 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
