import S11D4Invariants.Quadratic

/-! Full group predicates, elementary rotations and trace identities.
Generated coordinate identities are checked in Lean. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 1600000

def conjugate (R G : Mat) : Mat := R * G * R.transpose
def Orthogonal (R : Mat) : Prop := R.transpose * R = 1
def Proper (R : Mat) : Prop := Orthogonal R ∧ R.det = 1

def SOInvariant (Q : Quad) : Prop := ∀ R, Proper R → ∀ G, Q (conjugate R G) = Q G
def OInvariant (Q : Quad) : Prop := ∀ R, Orthogonal R → ∀ G, Q (conjugate R G) = Q G

def reflection : Mat := ![![-1,0,0,0],![0,1,0,0],![0,0,1,0],![0,0,0,1]]
def ReflectionOdd (Q : Quad) : Prop := ∀ G, Q (conjugate reflection G) = -Q G

theorem det_four (R : Mat) : R.det = -(R 0 0)*(R 2 2)*(R 3 1)*(R 1 3) + (R 0 0)*(R 2 2)*(R 3 3)*(R 1 1) + (R 0 0)*(R 2 3)*(R 3 1)*(R 1 2) - (R 0 0)*(R 2 3)*(R 3 2)*(R 1 1) + (R 0 0)*(R 3 2)*(R 1 3)*(R 2 1) - (R 0 0)*(R 3 3)*(R 1 2)*(R 2 1) + (R 0 1)*(R 2 2)*(R 3 0)*(R 1 3) - (R 0 1)*(R 2 2)*(R 3 3)*(R 1 0) - (R 0 1)*(R 2 3)*(R 3 0)*(R 1 2) + (R 0 1)*(R 2 3)*(R 3 2)*(R 1 0) - (R 0 1)*(R 3 2)*(R 1 3)*(R 2 0) + (R 0 1)*(R 3 3)*(R 1 2)*(R 2 0) - (R 2 2)*(R 3 0)*(R 0 3)*(R 1 1) + (R 2 2)*(R 3 1)*(R 0 3)*(R 1 0) + (R 2 3)*(R 3 0)*(R 0 2)*(R 1 1) - (R 2 3)*(R 3 1)*(R 0 2)*(R 1 0) - (R 3 0)*(R 0 2)*(R 1 3)*(R 2 1) + (R 3 0)*(R 0 3)*(R 1 2)*(R 2 1) + (R 3 1)*(R 0 2)*(R 1 3)*(R 2 0) - (R 3 1)*(R 0 3)*(R 1 2)*(R 2 0) - (R 3 2)*(R 0 3)*(R 1 0)*(R 2 1) + (R 3 2)*(R 0 3)*(R 1 1)*(R 2 0) + (R 3 3)*(R 0 2)*(R 1 0)*(R 2 1) - (R 3 3)*(R 0 2)*(R 1 1)*(R 2 0) := by
  rw [Matrix.det_succ_row_zero]
  simp [Fin.sum_univ_four, Matrix.det_fin_three, Matrix.submatrix_apply, Fin.succAbove,
    Fin.succ, Fin.castSucc]
  ring

def rotationXY (a b : ℝ) : Mat := ![![a,-b,0,0],![b,a,0,0],![0,0,1,0],![0,0,0,1]]

theorem rotationXY_proper {a b : ℝ} (h : a^2+b^2=1) : Proper (rotationXY a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
    fin_cases i <;> fin_cases j <;> simp [rotationXY] <;> nlinarith
  · rw [det_four]
    simp [rotationXY]
    nlinarith

theorem coordinates_rotationXY (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationXY a b) G) = ![a^2*(G 0 0) - a*b*(G 0 1) - a*b*(G 1 0) + b^2*(G 1 1),a^2*(G 0 1) + a*b*(G 0 0) - a*b*(G 1 1) - b^2*(G 1 0),a*(G 0 2) - b*(G 1 2),a*(G 0 3) - b*(G 1 3),a^2*(G 1 0) + a*b*(G 0 0) - a*b*(G 1 1) - b^2*(G 0 1),a^2*(G 1 1) + a*b*(G 0 1) + a*b*(G 1 0) + b^2*(G 0 0),a*(G 1 2) + b*(G 0 2),a*(G 1 3) + b*(G 0 3),a*(G 2 0) - b*(G 2 1),a*(G 2 1) + b*(G 2 0),(G 2 2),(G 2 3),a*(G 3 0) - b*(G 3 1),a*(G 3 1) + b*(G 3 0),(G 3 2),(G 3 3)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] <;>
    norm_num [rotationXY, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

def rotationYZ (a b : ℝ) : Mat := ![![1,0,0,0],![0,a,-b,0],![0,b,a,0],![0,0,0,1]]

theorem rotationYZ_proper {a b : ℝ} (h : a^2+b^2=1) : Proper (rotationYZ a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
    fin_cases i <;> fin_cases j <;> simp [rotationYZ] <;> nlinarith
  · rw [det_four]
    simp [rotationYZ]
    nlinarith

theorem coordinates_rotationYZ (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationYZ a b) G) = ![(G 0 0),a*(G 0 1) - b*(G 0 2),a*(G 0 2) + b*(G 0 1),(G 0 3),a*(G 1 0) - b*(G 2 0),a^2*(G 1 1) - a*b*(G 1 2) - a*b*(G 2 1) + b^2*(G 2 2),a^2*(G 1 2) - a*b*(G 2 2) + a*b*(G 1 1) - b^2*(G 2 1),a*(G 1 3) - b*(G 2 3),a*(G 2 0) + b*(G 1 0),a^2*(G 2 1) - a*b*(G 2 2) + a*b*(G 1 1) - b^2*(G 1 2),a^2*(G 2 2) + a*b*(G 1 2) + a*b*(G 2 1) + b^2*(G 1 1),a*(G 2 3) + b*(G 1 3),(G 3 0),a*(G 3 1) - b*(G 3 2),a*(G 3 2) + b*(G 3 1),(G 3 3)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] <;>
    norm_num [rotationYZ, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

def rotationZW (a b : ℝ) : Mat := ![![1,0,0,0],![0,1,0,0],![0,0,a,-b],![0,0,b,a]]

theorem rotationZW_proper {a b : ℝ} (h : a^2+b^2=1) : Proper (rotationZW a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
    fin_cases i <;> fin_cases j <;> simp [rotationZW] <;> nlinarith
  · rw [det_four]
    simp [rotationZW]
    nlinarith

theorem coordinates_rotationZW (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationZW a b) G) = ![(G 0 0),(G 0 1),a*(G 0 2) - b*(G 0 3),a*(G 0 3) + b*(G 0 2),(G 1 0),(G 1 1),a*(G 1 2) - b*(G 1 3),a*(G 1 3) + b*(G 1 2),a*(G 2 0) - b*(G 3 0),a*(G 2 1) - b*(G 3 1),a^2*(G 2 2) - a*b*(G 2 3) - a*b*(G 3 2) + b^2*(G 3 3),a^2*(G 2 3) + a*b*(G 2 2) - a*b*(G 3 3) - b^2*(G 3 2),a*(G 3 0) + b*(G 2 0),a*(G 3 1) + b*(G 2 1),a^2*(G 3 2) + a*b*(G 2 2) - a*b*(G 3 3) - b^2*(G 2 3),a^2*(G 3 3) + a*b*(G 2 3) + a*b*(G 3 2) + b^2*(G 2 2)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] <;>
    norm_num [rotationZW, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

theorem reflection_orthogonal : Orthogonal reflection := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
  fin_cases i <;> fin_cases j <;> norm_num [reflection, Matrix.cons_val_two, Matrix.cons_val_three]

theorem reflection_det : reflection.det = -1 := by
  rw [det_four]
  norm_num [reflection, Matrix.cons_val_two, Matrix.cons_val_three]

theorem trace_conjugate {R : Mat} (hR : Orthogonal R) (G : Mat) :
    (conjugate R G).trace = G.trace := by
  unfold conjugate
  rw [Matrix.trace_mul_cycle, hR, Matrix.one_mul]

theorem conjugate_mul {R : Mat} (hR : Orthogonal R) (G H : Mat) :
    conjugate R G * conjugate R H = conjugate R (G*H) := by
  unfold conjugate
  calc R*G*R.transpose*(R*H*R.transpose) = R*G*(R.transpose*R)*H*R.transpose := by
         simp only [Matrix.mul_assoc]
       _ = R*(G*H)*R.transpose := by
         rw [hR, Matrix.mul_one]
         simp only [Matrix.mul_assoc]

theorem conjugate_transpose (R G : Mat) :
    (conjugate R G).transpose = conjugate R G.transpose := by
  simp only [conjugate, Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.mul_assoc]

end
end S11D4Invariants
