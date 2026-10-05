import S11D5Invariants.ConstraintPolynomial

/-! A bounded block of necessary equations from explicit orthogonal rotations.
The rational certificate is split only to limit compilation resources. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem necessaryBlock16 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ((-1) * c 72 + (-1) * c 84 + (1) * c 94 + (1) * c 110 + (-1) * c 270 + (1) * c 310 = 0) ∧
    (((-544/625)) * c 0 + ((108/625)) * c 1 + ((108/625)) * c 5 + ((144/625)) * c 6 + ((144/625)) * c 25 + ((144/625)) * c 29 + ((192/625)) * c 30 + ((144/625)) * c 115 + ((192/625)) * c 116 + ((256/625)) * c 135 = 0) := by
  refine ⟨?_,?_⟩
  · have raw := hQ (rotationWV 0 1) (rotationWV_orthogonal (by norm_num)) (decode ![0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationWV, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  · have raw := hQ (rotationXY (3/5) (4/5)) (rotationXY_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    rw [hc, hc, coordinates_rotationXY, coordinates_decode] at raw
    simp only [decode, Matrix.cons_val] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]

end
end S11D5Invariants
