import S11D4Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit proper rotations.
This is the same rational certificate, split only to limit compilation resources. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock2 {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((1) * c 31 + (-1) * c 38 + (-1) * c 81 + (-1) * c 83 + (-1) * c 100 + (1) * c 108 = 0) ∧
    ((1) * c 31 + (-1) * c 39 + (-1) * c 81 + (-1) * c 85 = 0) ∧
    ((1) * c 31 + (-1) * c 40 + (-1) * c 81 + (-1) * c 86 = 0) ∧
    ((1) * c 31 + (-1) * c 42 + (-1) * c 81 + (-1) * c 87 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((1) * c 31 + (-1) * c 43 + (-1) * c 81 + (-1) * c 89 = 0) ∧
    ((1) * c 31 + (-1) * c 44 + (-1) * c 81 + (-1) * c 90 = 0) ∧
    ((1) * c 45 + (-1) * c 52 + (-1) * c 91 + (-1) * c 94 = 0) ∧
    ((1) * c 45 + (-1) * c 53 + (-1) * c 91 + (-1) * c 95 = 0) ∧
    ((1) * c 45 + (-1) * c 56 + (-1) * c 91 + (-1) * c 98 = 0) ∧
    ((1) * c 45 + (-1) * c 57 + (-1) * c 91 + (-1) * c 99 = 0) ∧
    ((-2) * c 101 = 0) ∧
    ((-1) * c 100 + (-1) * c 102 + (1) * c 108 + (1) * c 109 = 0) ∧
    ((-1) * c 100 + (-1) * c 103 + (1) * c 108 + (1) * c 110 = 0) ∧
    ((-1) * c 100 + (-1) * c 104 + (1) * c 108 + (1) * c 112 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 100 + (-1) * c 105 + (1) * c 108 + (-1) * c 111 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 100 + (-1) * c 106 + (1) * c 108 + (1) * c 113 = 0) ∧
    ((-1) * c 100 + (-1) * c 107 + (1) * c 108 + (1) * c 114 = 0) ∧
    ((1) * c 100 + (-1) * c 102 + (-1) * c 108 + (-1) * c 109 = 0) ∧
    ((1) * c 100 + (-1) * c 103 + (-1) * c 108 + (-1) * c 110 = 0) ∧
    ((1) * c 100 + (-1) * c 106 + (-1) * c 108 + (-1) * c 113 = 0) ∧
    ((1) * c 100 + (-1) * c 107 + (-1) * c 108 + (-1) * c 114 = 0) ∧
    ((-1) * c 117 + (1) * c 118 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 117 + (-1) * c 118 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-1) * c 122 + (1) * c 123 + (-1) * c 126 + (1) * c 130 = 0) ∧
    ((-1) * c 122 + (-1) * c 123 + (1) * c 126 + (-1) * c 130 = 0) ∧
    ((-2) * c 127 = 0) ∧
    ((-1) * c 126 + (-1) * c 128 + (1) * c 130 + (1) * c 131 = 0) ∧
    ((-1) * c 126 + (-1) * c 129 + (1) * c 130 + (1) * c 132 = 0) ∧
    ((1) * c 126 + (-1) * c 128 + (-1) * c 130 + (-1) * c 131 = 0) ∧
    ((1) * c 126 + (-1) * c 129 + (-1) * c 130 + (-1) * c 132 = 0) ∧
    ((-1) * c 1 + (1) * c 2 + (-1) * c 16 + (1) * c 31 = 0) ∧
    ((-1) * c 1 + (-1) * c 2 + (1) * c 16 + (-1) * c 31 = 0) ∧
    ((-1) * c 4 + (1) * c 8 + (-1) * c 58 + (1) * c 100 = 0) := by
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) (decode ![0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) (decode ![1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationYZ, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D4Invariants
