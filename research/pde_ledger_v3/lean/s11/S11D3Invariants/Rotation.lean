import S11D3Invariants.Quadratic

/-! Full group predicates and trace invariants on real 3x3 matrices.
Finite test rotations do not replace the quantifier over the full group. -/
namespace S11D3Invariants
noncomputable section

def conjugate (R G : Mat) : Mat := R * G * R.transpose
def Orthogonal (R : Mat) : Prop := R.transpose * R = 1
def Proper (R : Mat) : Prop := Orthogonal R ∧ R.det = 1

def SOInvariant (Q : Quad) : Prop := ∀ R, Proper R → ∀ G, Q (conjugate R G) = Q G
def OInvariant (Q : Quad) : Prop := ∀ R, Orthogonal R → ∀ G, Q (conjugate R G) = Q G

def reflection : Mat := ![![-1,0,0],![0,1,0],![0,0,1]]
def ReflectionOdd (Q : Quad) : Prop := ∀ G, Q (conjugate reflection G) = -Q G

def rotationXY (a b : ℝ) : Mat := ![![a,-b,0],![b,a,0],![0,0,1]]
def rotationYZ (a b : ℝ) : Mat := ![![1,0,0],![0,a,-b],![0,b,a]]

theorem rotationXY_proper {a b : ℝ} (h : a^2+b^2=1) : Proper (rotationXY a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]
    fin_cases i <;> fin_cases j <;> simp [rotationXY] <;> nlinarith
  · rw [Matrix.det_fin_three]
    simp [rotationXY]
    nlinarith

theorem rotationYZ_proper {a b : ℝ} (h : a^2+b^2=1) : Proper (rotationYZ a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]
    fin_cases i <;> fin_cases j <;> simp [rotationYZ] <;> nlinarith
  · rw [Matrix.det_fin_three]
    simp [rotationYZ]
    nlinarith

theorem reflection_orthogonal : Orthogonal reflection := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]
  fin_cases i <;> fin_cases j <;> norm_num [reflection, Matrix.cons_val_two]

theorem reflection_det : reflection.det = -1 := by
  rw [Matrix.det_fin_three]
  norm_num [reflection, Matrix.cons_val_two]

/-- Entry formula avoids reliance on simplification through concrete matrix literals. -/
theorem conjugate_apply (R G : Mat) (i j : Fin 3) :
    conjugate R G i j = ∑ k : Fin 3, (∑ l : Fin 3, R i l * G l k) * R j k := by
  rfl

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

/-- Unnormalized basis: trace square, trace of square, Frobenius square. -/
def traceSquare : Quad :=
  monomial 0 0 + monomial 4 4 + monomial 8 8 +
  2 • monomial 0 4 + 2 • monomial 0 8 + 2 • monomial 4 8

def traceOfSquare : Quad :=
  monomial 0 0 + monomial 4 4 + monomial 8 8 +
  2 • monomial 1 3 + 2 • monomial 2 6 + 2 • monomial 5 7

def frobeniusSquare : Quad :=
  monomial 0 0 + monomial 1 1 + monomial 2 2 + monomial 3 3 +
  monomial 4 4 + monomial 5 5 + monomial 6 6 + monomial 7 7 + monomial 8 8

theorem traceSquare_apply (G : Mat) : traceSquare G = G.trace ^ 2 := by
  simp [traceSquare, monomial, coordinateMap, coordinates, Matrix.trace, Fin.sum_univ_three]
  ring

theorem traceOfSquare_apply (G : Mat) : traceOfSquare G = (G*G).trace := by
  simp [traceOfSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_three]
  ring

theorem frobeniusSquare_apply (G : Mat) : frobeniusSquare G = (G*G.transpose).trace := by
  simp [frobeniusSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_three]
  ring

def invariantForm (v : Fin 3 → ℝ) : Quad :=
  v 0 • traceSquare + v 1 • traceOfSquare + v 2 • frobeniusSquare

theorem invariantForm_apply (v : Fin 3 → ℝ) (G : Mat) :
    invariantForm v G = v 0 * G.trace ^ 2 + v 1 * (G*G).trace +
      v 2 * (G*G.transpose).trace := by
  simp [invariantForm, traceSquare_apply, traceOfSquare_apply, frobeniusSquare_apply]

theorem invariantForm_O (v : Fin 3 → ℝ) : OInvariant (invariantForm v) := by
  intro R hR G
  simp only [invariantForm_apply, trace_conjugate hR, conjugate_transpose,
    conjugate_mul hR]

theorem invariantForm_SO (v : Fin 3 → ℝ) : SOInvariant (invariantForm v) :=
  fun R hR => invariantForm_O v R hR.1

end
end S11D3Invariants
