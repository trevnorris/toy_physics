import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock12 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 25 + (-1) * c 40 + (1) * c 49 + (1) * c 64 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 25 + (-1) * c 41 + (1) * c 49 + (-1) * c 63 + (1) * c 280 + (-1) * c 289 = 0) ∧
    ((-1) * c 25 + (-1) * c 42 + (1) * c 49 + (1) * c 65 = 0) ∧
    ((-1) * c 25 + (-1) * c 43 + (1) * c 49 + (1) * c 66 = 0) ∧
    ((-1) * c 25 + (-1) * c 44 + (1) * c 49 + (1) * c 67 = 0) ∧
    ((-1) * c 25 + (-1) * c 45 + (1) * c 49 + (1) * c 69 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 25 + (-1) * c 46 + (1) * c 49 + (-1) * c 68 + (1) * c 315 + (-1) * c 319 = 0) ∧
    ((-1) * c 25 + (-1) * c 47 + (1) * c 49 + (1) * c 70 = 0) ∧
    ((-1) * c 25 + (-1) * c 48 + (1) * c 49 + (1) * c 71 = 0) ∧
    ((-1) * c 74 + (1) * c 79 + (-1) * c 115 + (1) * c 205 = 0) ∧
    ((-1) * c 76 + (-1) * c 80 + (-1) * c 154 + (1) * c 220 = 0) ∧
    ((-1) * c 78 + (1) * c 83 + (-1) * c 189 + (1) * c 259 = 0) ∧
    ((-1) * c 85 + (1) * c 86 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 90 + (1) * c 91 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 95 + (1) * c 100 + (-1) * c 115 + (1) * c 205 = 0) ∧
    ((-1) * c 97 + (-1) * c 101 + (-1) * c 154 + (1) * c 220 = 0) ∧
    ((-1) * c 106 + (1) * c 107 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 111 + (1) * c 112 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 115 + (-1) * c 125 + (1) * c 205 + (1) * c 210 = 0) ∧
    ((-1) * c 115 + (-1) * c 127 + (1) * c 205 + (-1) * c 211 + (1) * c 280 + (-1) * c 289 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
