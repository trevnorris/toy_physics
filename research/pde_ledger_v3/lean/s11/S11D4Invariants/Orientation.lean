import S11D4Invariants.Rotation

/-! Four explicit quadratic forms and full-group sufficiency.
The orientation form is the Pfaffian polynomial of G-Gᵀ, with fixed
orientation (0,1,2,3); its determinant transformation is proved for all R. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

def orientation (G : Mat) : ℝ :=
  (G 0 1-G 1 0)*(G 2 3-G 3 2) - (G 0 2-G 2 0)*(G 1 3-G 3 1) +
  (G 0 3-G 3 0)*(G 1 2-G 2 1)

theorem orientation_conjugate (R G : Mat) :
    orientation (conjugate R G) = R.det * orientation G := by
  simp only [orientation, conjugate, Matrix.mul_apply, Matrix.transpose_apply,
    Fin.sum_univ_four, det_four]
  ring

end
end S11D4Invariants
