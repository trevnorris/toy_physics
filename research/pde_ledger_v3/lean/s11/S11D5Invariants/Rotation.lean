import S11D5Invariants.Quadratic

/-! Full-group predicates and the odd-dimensional conjugation identity.
The finite rotation constraints below are necessary conditions only. -/
namespace S11D5Invariants
noncomputable section

theorem sum_five {A : Type*} [AddCommMonoid A] (f : Fin 5 → A) :
    (∑ i, f i) = f 0 + f 1 + f 2 + f 3 + f 4 := by
  rw [Fin.sum_univ_succ]
  simp [Fin.sum_univ_four, add_assoc]

def conjugate (R G : Mat) : Mat := R * G * R.transpose
def Orthogonal (R : Mat) : Prop := R.transpose * R = 1
def Proper (R : Mat) : Prop := Orthogonal R ∧ R.det = 1

def SOInvariant (Q : Quad) : Prop := ∀ R, Proper R → ∀ G, Q (conjugate R G) = Q G
def OInvariant (Q : Quad) : Prop := ∀ R, Orthogonal R → ∀ G, Q (conjugate R G) = Q G

def reflection : Mat := Matrix.diagonal ![-1,1,1,1,1]
def ReflectionOdd (Q : Quad) : Prop := ∀ G, Q (conjugate reflection G) = -Q G

theorem conjugate_neg (R G : Mat) : conjugate (-R) G = conjugate R G := by
  simp only [conjugate, Matrix.transpose_neg, neg_mul, mul_neg, neg_neg]

theorem orthogonal_det {R : Mat} (hR : Orthogonal R) : R.det = 1 ∨ R.det = -1 := by
  have hd : R.det * R.det = 1 := by
    simpa only [Matrix.det_mul, Matrix.det_transpose, Matrix.det_one] using congrArg Matrix.det hR
  exact mul_self_eq_one_iff.mp hd

theorem SO_iff_O (Q : Quad) : SOInvariant Q ↔ OInvariant Q := by
  constructor
  · intro hQ R hR G
    rcases orthogonal_det hR with hd | hd
    · exact hQ R ⟨hR,hd⟩ G
    · have hn : Proper (-R) := by
        constructor
        · simpa only [Orthogonal, Matrix.transpose_neg, neg_mul_neg] using hR
        · rw [Matrix.det_neg]
          norm_num [hd]
      simpa only [conjugate_neg] using hQ (-R) hn G
  · intro hQ R hR G
    exact hQ R hR.1 G

def rotationXY (a b : ℝ) : Mat := ![![a,-b,0,0,0],![b,a,0,0,0],![0,0,1,0,0],![0,0,0,1,0],![0,0,0,0,1]]

theorem rotationXY_orthogonal {a b : ℝ} (h : a^2+b^2=1) : Orthogonal (rotationXY a b) := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  fin_cases i <;> fin_cases j <;>
    norm_num [rotationXY, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> nlinarith

theorem coordinates_rotationXY (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationXY a b) G) = ![a^2*(G 0 0) - a*b*(G 0 1) - a*b*(G 1 0) + b^2*(G 1 1),a^2*(G 0 1) + a*b*(G 0 0) - a*b*(G 1 1) - b^2*(G 1 0),a*(G 0 2) - b*(G 1 2),a*(G 0 3) - b*(G 1 3),a*(G 0 4) - b*(G 1 4),a^2*(G 1 0) + a*b*(G 0 0) - a*b*(G 1 1) - b^2*(G 0 1),a^2*(G 1 1) + a*b*(G 0 1) + a*b*(G 1 0) + b^2*(G 0 0),a*(G 1 2) + b*(G 0 2),a*(G 1 3) + b*(G 0 3),a*(G 1 4) + b*(G 0 4),a*(G 2 0) - b*(G 2 1),a*(G 2 1) + b*(G 2 0),(G 2 2),(G 2 3),(G 2 4),a*(G 3 0) - b*(G 3 1),a*(G 3 1) + b*(G 3 0),(G 3 2),(G 3 3),(G 3 4),a*(G 4 0) - b*(G 4 1),a*(G 4 1) + b*(G 4 0),(G 4 2),(G 4 3),(G 4 4)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, sum_five] <;>
    norm_num [rotationXY, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> ring

def rotationYZ (a b : ℝ) : Mat := ![![1,0,0,0,0],![0,a,-b,0,0],![0,b,a,0,0],![0,0,0,1,0],![0,0,0,0,1]]

theorem rotationYZ_orthogonal {a b : ℝ} (h : a^2+b^2=1) : Orthogonal (rotationYZ a b) := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  fin_cases i <;> fin_cases j <;>
    norm_num [rotationYZ, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> nlinarith

theorem coordinates_rotationYZ (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationYZ a b) G) = ![(G 0 0),a*(G 0 1) - b*(G 0 2),a*(G 0 2) + b*(G 0 1),(G 0 3),(G 0 4),a*(G 1 0) - b*(G 2 0),a^2*(G 1 1) - a*b*(G 2 1) - a*b*(G 1 2) + b^2*(G 2 2),a^2*(G 1 2) - a*b*(G 2 2) + a*b*(G 1 1) - b^2*(G 2 1),a*(G 1 3) - b*(G 2 3),a*(G 1 4) - b*(G 2 4),a*(G 2 0) + b*(G 1 0),a^2*(G 2 1) - a*b*(G 2 2) + a*b*(G 1 1) - b^2*(G 1 2),a^2*(G 2 2) + a*b*(G 2 1) + a*b*(G 1 2) + b^2*(G 1 1),a*(G 2 3) + b*(G 1 3),a*(G 2 4) + b*(G 1 4),(G 3 0),a*(G 3 1) - b*(G 3 2),a*(G 3 2) + b*(G 3 1),(G 3 3),(G 3 4),(G 4 0),a*(G 4 1) - b*(G 4 2),a*(G 4 2) + b*(G 4 1),(G 4 3),(G 4 4)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, sum_five] <;>
    norm_num [rotationYZ, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> ring

def rotationZW (a b : ℝ) : Mat := ![![1,0,0,0,0],![0,1,0,0,0],![0,0,a,-b,0],![0,0,b,a,0],![0,0,0,0,1]]

theorem rotationZW_orthogonal {a b : ℝ} (h : a^2+b^2=1) : Orthogonal (rotationZW a b) := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  fin_cases i <;> fin_cases j <;>
    norm_num [rotationZW, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> nlinarith

theorem coordinates_rotationZW (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationZW a b) G) = ![(G 0 0),(G 0 1),a*(G 0 2) - b*(G 0 3),a*(G 0 3) + b*(G 0 2),(G 0 4),(G 1 0),(G 1 1),a*(G 1 2) - b*(G 1 3),a*(G 1 3) + b*(G 1 2),(G 1 4),a*(G 2 0) - b*(G 3 0),a*(G 2 1) - b*(G 3 1),a^2*(G 2 2) - a*b*(G 2 3) - a*b*(G 3 2) + b^2*(G 3 3),a^2*(G 2 3) + a*b*(G 2 2) - a*b*(G 3 3) - b^2*(G 3 2),a*(G 2 4) - b*(G 3 4),a*(G 3 0) + b*(G 2 0),a*(G 3 1) + b*(G 2 1),a^2*(G 3 2) + a*b*(G 2 2) - a*b*(G 3 3) - b^2*(G 2 3),a^2*(G 3 3) + a*b*(G 2 3) + a*b*(G 3 2) + b^2*(G 2 2),a*(G 3 4) + b*(G 2 4),(G 4 0),(G 4 1),a*(G 4 2) - b*(G 4 3),a*(G 4 3) + b*(G 4 2),(G 4 4)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, sum_five] <;>
    norm_num [rotationZW, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> ring

def rotationWV (a b : ℝ) : Mat := ![![1,0,0,0,0],![0,1,0,0,0],![0,0,1,0,0],![0,0,0,a,-b],![0,0,0,b,a]]

theorem rotationWV_orthogonal {a b : ℝ} (h : a^2+b^2=1) : Orthogonal (rotationWV a b) := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  fin_cases i <;> fin_cases j <;>
    norm_num [rotationWV, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> nlinarith

theorem coordinates_rotationWV (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotationWV a b) G) = ![(G 0 0),(G 0 1),(G 0 2),a*(G 0 3) - b*(G 0 4),a*(G 0 4) + b*(G 0 3),(G 1 0),(G 1 1),(G 1 2),a*(G 1 3) - b*(G 1 4),a*(G 1 4) + b*(G 1 3),(G 2 0),(G 2 1),(G 2 2),a*(G 2 3) - b*(G 2 4),a*(G 2 4) + b*(G 2 3),a*(G 3 0) - b*(G 4 0),a*(G 3 1) - b*(G 4 1),a*(G 3 2) - b*(G 4 2),a^2*(G 3 3) - a*b*(G 3 4) - a*b*(G 4 3) + b^2*(G 4 4),a^2*(G 3 4) + a*b*(G 3 3) - a*b*(G 4 4) - b^2*(G 4 3),a*(G 4 0) + b*(G 3 0),a*(G 4 1) + b*(G 3 1),a*(G 4 2) + b*(G 3 2),a^2*(G 4 3) + a*b*(G 3 3) - a*b*(G 4 4) - b^2*(G 3 4),a^2*(G 4 4) + a*b*(G 3 4) + a*b*(G 4 3) + b^2*(G 3 3)] := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, sum_five] <;>
    norm_num [rotationWV, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] <;> ring

theorem reflection_orthogonal : Orthogonal reflection := by
  unfold Orthogonal reflection
  rw [Matrix.diagonal_transpose, Matrix.diagonal_mul_diagonal]
  change Matrix.diagonal (fun i : Fin 5 =>
    (![-1,1,1,1,1] : Fin 5 → ℝ) i * (![-1,1,1,1,1] : Fin 5 → ℝ) i) =
    Matrix.diagonal (fun _ : Fin 5 => (1 : ℝ))
  apply congrArg Matrix.diagonal
  funext i
  fin_cases i <;> norm_num [Matrix.cons_val]

theorem reflection_det : reflection.det = -1 := by
  norm_num [reflection, Matrix.det_diagonal, Fin.prod_univ_succ]

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
end S11D5Invariants
