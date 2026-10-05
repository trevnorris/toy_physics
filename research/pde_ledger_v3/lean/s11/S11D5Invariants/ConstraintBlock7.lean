import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock7 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((1) * c 72 + (-1) * c 82 + (-1) * c 172 + (-1) * c 177 = 0) ∧
    ((1) * c 72 + (-1) * c 83 + (-1) * c 172 + (-1) * c 178 = 0) ∧
    ((1) * c 72 + (-1) * c 86 + (-1) * c 172 + (-1) * c 181 = 0) ∧
    ((1) * c 72 + (-1) * c 87 + (-1) * c 172 + (-1) * c 182 = 0) ∧
    ((1) * c 72 + (-1) * c 88 + (-1) * c 172 + (-1) * c 183 = 0) ∧
    ((1) * c 72 + (-1) * c 91 + (-1) * c 172 + (-1) * c 186 = 0) ∧
    ((1) * c 72 + (-1) * c 92 + (-1) * c 172 + (-1) * c 187 = 0) ∧
    ((1) * c 72 + (-1) * c 93 + (-1) * c 172 + (-1) * c 188 = 0) ∧
    ((1) * c 94 + (-1) * c 102 + (-1) * c 189 + (-1) * c 192 = 0) ∧
    ((1) * c 94 + (-1) * c 103 + (-1) * c 189 + (-1) * c 193 = 0) ∧
    ((1) * c 94 + (-1) * c 104 + (-1) * c 189 + (-1) * c 194 = 0) ∧
    ((1) * c 94 + (-1) * c 107 + (-1) * c 189 + (-1) * c 197 = 0) ∧
    ((1) * c 94 + (-1) * c 108 + (-1) * c 189 + (-1) * c 198 = 0) ∧
    ((1) * c 94 + (-1) * c 109 + (-1) * c 189 + (-1) * c 199 = 0) ∧
    ((1) * c 94 + (-1) * c 112 + (-1) * c 189 + (-1) * c 202 = 0) ∧
    ((1) * c 94 + (-1) * c 113 + (-1) * c 189 + (-1) * c 203 = 0) ∧
    ((1) * c 94 + (-1) * c 114 + (-1) * c 189 + (-1) * c 204 = 0) ∧
    ((-2) * c 206 = 0) ∧
    ((-1) * c 205 + (-1) * c 207 + (1) * c 220 + (1) * c 221 = 0) ∧
    ((-1) * c 205 + (-1) * c 208 + (1) * c 220 + (1) * c 222 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
