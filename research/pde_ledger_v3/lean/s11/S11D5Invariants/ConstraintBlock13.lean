import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock13 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 115 + (-1) * c 130 + (1) * c 205 + (1) * c 215 = 0) ∧
    ((-1) * c 115 + (-1) * c 132 + (1) * c 205 + (-1) * c 216 + (1) * c 315 + (-1) * c 319 = 0) ∧
    ((-1) * c 135 + (-1) * c 137 + (-1) * c 172 + (1) * c 234 + (1) * c 235 + (1) * c 247 = 0) ∧
    ((-1) * c 135 + (-1) * c 138 + (-1) * c 189 + (1) * c 234 + (1) * c 236 + (1) * c 259 = 0) ∧
    ((-1) * c 135 + (-1) * c 145 + (1) * c 234 + (1) * c 239 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 135 + (-1) * c 147 + (1) * c 234 + (1) * c 240 = 0) ∧
    ((-1) * c 135 + (-1) * c 148 + (1) * c 234 + (1) * c 241 = 0) ∧
    ((-1) * c 135 + (-1) * c 150 + (1) * c 234 + (1) * c 244 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 135 + (-1) * c 152 + (1) * c 234 + (1) * c 245 = 0) ∧
    ((-1) * c 135 + (-1) * c 153 + (1) * c 234 + (1) * c 246 = 0) ∧
    ((-1) * c 172 + (-1) * c 173 + (-1) * c 189 + (1) * c 247 + (1) * c 248 + (1) * c 259 = 0) ∧
    ((-1) * c 172 + (-1) * c 180 + (1) * c 247 + (1) * c 251 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 172 + (-1) * c 182 + (1) * c 247 + (1) * c 252 = 0) ∧
    ((-1) * c 172 + (-1) * c 183 + (1) * c 247 + (1) * c 253 = 0) ∧
    ((-1) * c 172 + (-1) * c 185 + (1) * c 247 + (1) * c 256 + (-1) * c 315 + (1) * c 319 = 0) ∧
    ((-1) * c 172 + (-1) * c 187 + (1) * c 247 + (1) * c 257 = 0) ∧
    ((-1) * c 172 + (-1) * c 188 + (1) * c 247 + (1) * c 258 = 0) ∧
    ((-1) * c 189 + (-1) * c 196 + (1) * c 259 + (1) * c 262 + (-1) * c 280 + (1) * c 289 = 0) ∧
    ((-1) * c 189 + (-1) * c 198 + (1) * c 259 + (1) * c 263 = 0) ∧
    ((-1) * c 189 + (-1) * c 199 + (1) * c 259 + (1) * c 264 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_orthogonal (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
