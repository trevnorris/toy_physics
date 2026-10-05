import S11D4Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit proper rotations.
This is the same rational certificate, split only to limit compilation resources. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock0 {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 0 + (1) * c 70 = 0) ∧
    ((-1) * c 0 + (-1) * c 1 + (-1) * c 16 + (1) * c 58 + (-1) * c 59 + (1) * c 70 = 0) ∧
    ((-1) * c 0 + (-1) * c 2 + (-1) * c 31 + (1) * c 70 + (1) * c 71 + (1) * c 81 = 0) ∧
    ((-1) * c 0 + (-1) * c 3 + (-1) * c 45 + (1) * c 70 + (1) * c 72 + (1) * c 91 = 0) ∧
    ((-1) * c 0 + (-1) * c 4 + (1) * c 16 + (-1) * c 20 + (-1) * c 58 + (1) * c 70 = 0) ∧
    ((-1) * c 0 + (-1) * c 6 + (1) * c 31 + (-1) * c 34 + (1) * c 70 + (-1) * c 81 = 0) ∧
    ((-1) * c 0 + (-1) * c 7 + (1) * c 45 + (-1) * c 47 + (1) * c 70 + (-1) * c 91 = 0) ∧
    ((-1) * c 0 + (-1) * c 8 + (1) * c 70 + (1) * c 74 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((-1) * c 0 + (-1) * c 9 + (1) * c 70 + (-1) * c 73 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((-1) * c 0 + (-1) * c 10 + (1) * c 70 + (1) * c 75 = 0) ∧
    ((-1) * c 0 + (-1) * c 11 + (1) * c 70 + (1) * c 76 = 0) ∧
    ((-1) * c 0 + (-1) * c 12 + (1) * c 70 + (1) * c 78 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 0 + (-1) * c 13 + (1) * c 70 + (-1) * c 77 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 0 + (-1) * c 14 + (1) * c 70 + (1) * c 79 = 0) ∧
    ((-1) * c 0 + (-1) * c 15 + (1) * c 70 + (1) * c 80 = 0) ∧
    ((-1) * c 16 + (1) * c 58 = 0) ∧
    ((-1) * c 16 + (-1) * c 17 + (-1) * c 31 + (1) * c 58 + (-1) * c 60 + (1) * c 81 = 0) ∧
    ((-1) * c 16 + (-1) * c 18 + (-1) * c 45 + (1) * c 58 + (-1) * c 61 + (1) * c 91 = 0) ∧
    ((-1) * c 16 + (-1) * c 21 + (1) * c 31 + (1) * c 33 + (1) * c 58 + (-1) * c 81 = 0) ∧
    ((-1) * c 16 + (-1) * c 22 + (1) * c 45 + (1) * c 46 + (1) * c 58 + (-1) * c 91 = 0) ∧
    ((-1) * c 16 + (-1) * c 23 + (1) * c 58 + (-1) * c 63 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((-1) * c 16 + (-1) * c 24 + (1) * c 58 + (1) * c 62 + (1) * c 100 + (-1) * c 108 = 0) ∧
    ((-1) * c 16 + (-1) * c 25 + (1) * c 58 + (-1) * c 64 = 0) ∧
    ((-1) * c 16 + (-1) * c 26 + (1) * c 58 + (-1) * c 65 = 0) ∧
    ((-1) * c 16 + (-1) * c 27 + (1) * c 58 + (-1) * c 67 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 16 + (-1) * c 28 + (1) * c 58 + (1) * c 66 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 16 + (-1) * c 29 + (1) * c 58 + (-1) * c 68 = 0) ∧
    ((-1) * c 16 + (-1) * c 30 + (1) * c 58 + (-1) * c 69 = 0) ∧
    ((-1) * c 31 + (1) * c 81 = 0) ∧
    ((-1) * c 31 + (-1) * c 32 + (-1) * c 45 + (1) * c 81 + (1) * c 82 + (1) * c 91 = 0) ∧
    ((1) * c 16 + (-1) * c 21 + (-1) * c 31 + (-1) * c 33 + (-1) * c 58 + (1) * c 81 = 0) ∧
    ((1) * c 0 + (1) * c 6 + (-1) * c 31 + (-1) * c 34 + (-1) * c 70 + (1) * c 81 = 0) ∧
    ((-2) * c 35 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D4Invariants
