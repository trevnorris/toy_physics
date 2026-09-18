import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock11 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 1 + (1) * c 2 + (-1) * c 25 + (1) * c 49 = 0) ∧
    ((-1) * c 1 + (-1) * c 2 + (1) * c 25 + (-1) * c 49 = 0) ∧
    ((-1) * c 5 + (1) * c 10 + (-1) * c 115 + (1) * c 205 = 0) ∧
    ((-1) * c 6 + (1) * c 12 + (-1) * c 135 + (1) * c 234 = 0) ∧
    ((-1) * c 7 + (-1) * c 11 + (-1) * c 154 + (1) * c 220 = 0) ∧
    ((-1) * c 8 + (1) * c 13 + (-1) * c 172 + (1) * c 247 = 0) ∧
    ((-1) * c 9 + (1) * c 14 + (-1) * c 189 + (1) * c 259 = 0) ∧
    ((-1) * c 16 + (1) * c 17 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 21 + (1) * c 22 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 25 + (-1) * c 27 + (1) * c 49 + (1) * c 50 = 0) ∧
    ((-1) * c 25 + (-1) * c 28 + (1) * c 49 + (1) * c 51 = 0) ∧
    ((-1) * c 25 + (-1) * c 29 + (1) * c 49 + (1) * c 57 + (-1) * c 115 + (1) * c 205 = 0) ∧
    ((-1) * c 25 + (-1) * c 30 + (1) * c 49 + (1) * c 59 + (-1) * c 135 + (1) * c 234 = 0) ∧
    ((-1) * c 25 + (-1) * c 31 + (1) * c 49 + (-1) * c 58 + (-1) * c 154 + (1) * c 220 = 0) ∧
    ((-1) * c 25 + (-1) * c 32 + (1) * c 49 + (1) * c 60 + (-1) * c 172 + (1) * c 247 = 0) ∧
    ((-1) * c 25 + (-1) * c 33 + (1) * c 49 + (1) * c 61 + (-1) * c 189 + (1) * c 259 = 0) ∧
    ((-1) * c 25 + (-1) * c 36 + (1) * c 49 + (1) * c 53 + (1) * c 135 + (-1) * c 234 = 0) ∧
    ((-1) * c 25 + (-1) * c 37 + (1) * c 49 + (-1) * c 55 + (1) * c 172 + (-1) * c 247 = 0) ∧
    ((-1) * c 25 + (-1) * c 38 + (1) * c 49 + (-1) * c 56 + (1) * c 189 + (-1) * c 259 = 0) ∧
    ((-1) * c 25 + (-1) * c 39 + (1) * c 49 + (1) * c 62 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
