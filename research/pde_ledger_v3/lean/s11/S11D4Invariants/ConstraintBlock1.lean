import S11D4Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit proper rotations.
This is the same rational certificate, split only to limit compilation resources. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock1 {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 31 + (-1) * c 36 + (1) * c 45 + (-1) * c 48 + (1) * c 81 + (-1) * c 91 = 0) ∧
    ((-1) * c 31 + (-1) * c 37 + (1) * c 81 + (1) * c 84 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((-1) * c 31 + (-1) * c 38 + (1) * c 81 + (-1) * c 83 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((-1) * c 31 + (-1) * c 39 + (1) * c 81 + (1) * c 85 = 0) ∧
    ((-1) * c 31 + (-1) * c 40 + (1) * c 81 + (1) * c 86 = 0) ∧
    ((-1) * c 31 + (-1) * c 41 + (1) * c 81 + (1) * c 88 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 31 + (-1) * c 42 + (1) * c 81 + (-1) * c 87 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 31 + (-1) * c 43 + (1) * c 81 + (1) * c 89 = 0) ∧
    ((-1) * c 31 + (-1) * c 44 + (1) * c 81 + (1) * c 90 = 0) ∧
    ((-1) * c 45 + (1) * c 91 = 0) ∧
    ((1) * c 16 + (-1) * c 22 + (-1) * c 45 + (-1) * c 46 + (-1) * c 58 + (1) * c 91 = 0) ∧
    ((1) * c 0 + (1) * c 7 + (-1) * c 45 + (-1) * c 47 + (-1) * c 70 + (1) * c 91 = 0) ∧
    ((-2) * c 49 = 0) ∧
    ((-1) * c 45 + (-1) * c 50 + (1) * c 91 + (1) * c 93 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((-1) * c 45 + (-1) * c 51 + (1) * c 91 + (-1) * c 92 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((-1) * c 45 + (-1) * c 52 + (1) * c 91 + (1) * c 94 = 0) ∧
    ((-1) * c 45 + (-1) * c 53 + (1) * c 91 + (1) * c 95 = 0) ∧
    ((-1) * c 45 + (-1) * c 54 + (1) * c 91 + (1) * c 97 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 45 + (-1) * c 55 + (1) * c 91 + (-1) * c 96 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 45 + (-1) * c 56 + (1) * c 91 + (1) * c 98 = 0) ∧
    ((-1) * c 45 + (-1) * c 57 + (1) * c 91 + (1) * c 99 = 0) ∧
    ((1) * c 16 + (1) * c 17 + (1) * c 31 + (-1) * c 58 + (-1) * c 60 + (-1) * c 81 = 0) ∧
    ((1) * c 16 + (1) * c 18 + (1) * c 45 + (-1) * c 58 + (-1) * c 61 + (-1) * c 91 = 0) ∧
    ((1) * c 16 + (-1) * c 24 + (-1) * c 58 + (-1) * c 62 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((1) * c 16 + (1) * c 23 + (-1) * c 58 + (-1) * c 63 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((1) * c 16 + (-1) * c 28 + (-1) * c 58 + (-1) * c 66 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((1) * c 16 + (1) * c 27 + (-1) * c 58 + (-1) * c 67 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((1) * c 0 + (-1) * c 2 + (1) * c 31 + (-1) * c 70 + (-1) * c 71 + (-1) * c 81 = 0) ∧
    ((1) * c 0 + (-1) * c 3 + (1) * c 45 + (-1) * c 70 + (-1) * c 72 + (-1) * c 91 = 0) ∧
    ((1) * c 0 + (1) * c 9 + (-1) * c 70 + (-1) * c 73 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((1) * c 0 + (-1) * c 8 + (-1) * c 70 + (-1) * c 74 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((1) * c 0 + (1) * c 13 + (-1) * c 70 + (-1) * c 77 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((1) * c 0 + (-1) * c 12 + (-1) * c 70 + (-1) * c 78 + (1) * c 126 + (-1) * c 130 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D4Invariants
