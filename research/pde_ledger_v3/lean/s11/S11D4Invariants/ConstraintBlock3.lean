import S11D4Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit proper rotations.
This is the same rational certificate, split only to limit compilation resources. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock3 {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 5 + (1) * c 10 + (-1) * c 70 + (1) * c 115 = 0) ∧
    ((-1) * c 6 + (-1) * c 9 + (-1) * c 81 + (1) * c 108 = 0) ∧
    ((-1) * c 7 + (1) * c 11 + (-1) * c 91 + (1) * c 121 = 0) ∧
    ((-1) * c 13 + (1) * c 14 + (-1) * c 130 + (1) * c 133 = 0) ∧
    ((-1) * c 16 + (-1) * c 18 + (1) * c 31 + (1) * c 32 = 0) ∧
    ((-1) * c 16 + (-1) * c 19 + (1) * c 31 + (1) * c 37 + (-1) * c 58 + (1) * c 100 = 0) ∧
    ((-1) * c 16 + (-1) * c 20 + (1) * c 31 + (1) * c 39 + (-1) * c 70 + (1) * c 115 = 0) ∧
    ((-1) * c 16 + (-1) * c 21 + (1) * c 31 + (-1) * c 38 + (-1) * c 81 + (1) * c 108 = 0) ∧
    ((-1) * c 16 + (-1) * c 22 + (1) * c 31 + (1) * c 40 + (-1) * c 91 + (1) * c 121 = 0) ∧
    ((-1) * c 16 + (-1) * c 25 + (1) * c 31 + (1) * c 34 + (1) * c 70 + (-1) * c 115 = 0) ∧
    ((-1) * c 16 + (-1) * c 26 + (1) * c 31 + (-1) * c 36 + (1) * c 91 + (-1) * c 121 = 0) ∧
    ((-1) * c 16 + (-1) * c 27 + (1) * c 31 + (1) * c 41 = 0) ∧
    ((-1) * c 16 + (-1) * c 28 + (1) * c 31 + (1) * c 43 + (-1) * c 130 + (1) * c 133 = 0) ∧
    ((-1) * c 16 + (-1) * c 29 + (1) * c 31 + (-1) * c 42 + (1) * c 130 + (-1) * c 133 = 0) ∧
    ((-1) * c 16 + (-1) * c 30 + (1) * c 31 + (1) * c 44 = 0) ∧
    ((-1) * c 46 + (1) * c 50 + (-1) * c 58 + (1) * c 100 = 0) ∧
    ((-1) * c 48 + (-1) * c 51 + (-1) * c 81 + (1) * c 108 = 0) ∧
    ((-1) * c 55 + (1) * c 56 + (-1) * c 130 + (1) * c 133 = 0) ∧
    ((-1) * c 58 + (-1) * c 66 + (1) * c 100 + (1) * c 104 = 0) ∧
    ((-1) * c 58 + (-1) * c 68 + (1) * c 100 + (-1) * c 105 + (1) * c 130 + (-1) * c 133 = 0) ∧
    ((-1) * c 70 + (-1) * c 72 + (-1) * c 91 + (1) * c 115 + (1) * c 116 + (1) * c 121 = 0) ∧
    ((-1) * c 70 + (-1) * c 78 + (1) * c 115 + (1) * c 119 + (-1) * c 130 + (1) * c 133 = 0) ∧
    ((-1) * c 70 + (-1) * c 80 + (1) * c 115 + (1) * c 120 = 0) ∧
    ((-1) * c 91 + (-1) * c 97 + (1) * c 121 + (1) * c 124 + (-1) * c 130 + (1) * c 133 = 0) ∧
    ((-1) * c 91 + (-1) * c 99 + (1) * c 121 + (1) * c 125 = 0) ∧
    ((-1) * c 130 + (-1) * c 132 + (1) * c 133 + (1) * c 134 = 0) ∧
    ((-1) * c 2 + (1) * c 3 + (-1) * c 31 + (1) * c 45 = 0) ∧
    ((-1) * c 8 + (1) * c 12 + (-1) * c 100 + (1) * c 126 = 0) ∧
    ((-1) * c 10 + (1) * c 15 + (-1) * c 115 + (1) * c 135 = 0) ∧
    ((-1) * c 25 + (1) * c 30 + (-1) * c 115 + (1) * c 135 = 0) ∧
    ((-1) * c 26 + (-1) * c 29 + (-1) * c 121 + (1) * c 133 = 0) ∧
    ((-1) * c 31 + (-1) * c 37 + (1) * c 45 + (1) * c 54 + (-1) * c 100 + (1) * c 126 = 0) ∧
    (((-544/625)) * c 0 + ((108/625)) * c 1 + ((108/625)) * c 4 + ((144/625)) * c 5 + ((144/625)) * c 16 + ((144/625)) * c 19 + ((192/625)) * c 20 + ((144/625)) * c 58 + ((192/625)) * c 59 + ((256/625)) * c 70 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationZW 0 1) (rotationZW_proper (by norm_num)) (decode ![0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationZW, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY (3/5) (4/5)) (rotationXY_proper (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D4Invariants
