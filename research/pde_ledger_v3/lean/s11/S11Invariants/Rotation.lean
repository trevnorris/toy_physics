import S11Invariants.Quadratic

/-! Full orthogonal and proper-orthogonal conjugations, with no sampled
replacement for the group quantifier. -/

namespace S11Invariants
noncomputable section

def conjugate (R G : Mat) : Mat := R * G * R.transpose
def Orthogonal (R : Mat) : Prop := R.transpose * R = 1
def Proper (R : Mat) : Prop := Orthogonal R ∧ R.det = 1
def rotation (a b : ℝ) : Mat := ![![a,-b],![b,a]]
def reflection : Mat := ![![-1,0],![0,1]]

def SOInvariant (Q : Quad) : Prop := ∀ R, Proper R → ∀ G, Q (conjugate R G) = Q G
def OInvariant (Q : Quad) : Prop := ∀ R, Orthogonal R → ∀ G, Q (conjugate R G) = Q G
def ReflectionOdd (Q : Quad) : Prop := ∀ G, Q (conjugate reflection G) = -Q G

theorem rotation_proper {a b : ℝ} (h : a ^ 2 + b ^ 2 = 1) : Proper (rotation a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two]
    fin_cases i <;> fin_cases j <;> simp [rotation] <;> nlinarith
  · rw [Matrix.det_fin_two]
    simp [rotation]
    nlinarith

theorem proper_rotation {R : Mat} (h : Proper R) :
    ∃ a b : ℝ, a ^ 2 + b ^ 2 = 1 ∧ R = rotation a b := by
  have h00 := congrArg (fun M : Mat => M 0 0) h.1
  have h01 := congrArg (fun M : Mat => M 0 1) h.1
  have hd := h.2
  simp [Matrix.mul_apply, Fin.sum_univ_two] at h00 h01
  simp [Matrix.det_fin_two] at hd
  have hda : R 1 1 = R 0 0 := by
    nlinarith [congrArg (fun z => R 1 1 * z) h00,
      congrArg (fun z => R 1 0 * z) h01, congrArg (fun z => R 0 0 * z) hd]
  have hbc : R 0 1 = -R 1 0 := by
    nlinarith [congrArg (fun z => R 0 1 * z) h00,
      congrArg (fun z => R 0 0 * z) h01, congrArg (fun z => R 1 0 * z) hd]
  refine ⟨R 0 0, R 1 0, by nlinarith, ?_⟩
  ext i j
  fin_cases i <;> fin_cases j <;> simp [rotation, hda, hbc]

theorem coordinates_rotation (a b : ℝ) (G : Mat) :
    coordinates (conjugate (rotation a b) G) =
    ![(a^2+b^2)*coordinates G 0, (a^2+b^2)*coordinates G 1,
      (a^2-b^2)*coordinates G 2 - 2*a*b*coordinates G 3,
      2*a*b*coordinates G 2 + (a^2-b^2)*coordinates G 3] := by
  ext i
  fin_cases i <;> simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Matrix.cons_val] <;> simp [rotation, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

theorem spin_two_norm (a b x y : ℝ) :
    ((a^2-b^2)*x - 2*a*b*y)^2 + (2*a*b*x + (a^2-b^2)*y)^2 =
      (a^2+b^2)^2 * (x^2+y^2) := by ring

theorem reflection_orthogonal : Orthogonal reflection := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two]
  fin_cases i <;> fin_cases j <;> norm_num [reflection]

theorem reflection_det : reflection.det = -1 := by rw [Matrix.det_fin_two]; norm_num [reflection]

theorem reflection_square : reflection * reflection = 1 := by
  ext i j
  simp only [Matrix.mul_apply, Fin.sum_univ_two]
  fin_cases i <;> fin_cases j <;> norm_num [reflection]

theorem coordinates_reflection (G : Mat) :
    coordinates (conjugate reflection G) =
      ![coordinates G 0, -coordinates G 1, coordinates G 2, -coordinates G 3] := by
  ext i
  fin_cases i <;> simp only [coordinates, conjugate, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Matrix.cons_val] <;> simp [reflection, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

theorem orthogonal_det {R : Mat} (h : Orthogonal R) : R.det = 1 ∨ R.det = -1 := by
  have hd := congrArg Matrix.det h
  simp only [Matrix.det_mul, Matrix.det_transpose, Matrix.det_one] at hd
  have hf : (R.det - 1) * (R.det + 1) = 0 := by nlinarith
  rcases mul_eq_zero.mp hf with hp | hn
  · left; linarith
  · right; linarith

theorem orthogonal_mul {R S : Mat} (hR : Orthogonal R) (hS : Orthogonal S) :
    Orthogonal (R*S) := by
  unfold Orthogonal at *
  rw [Matrix.transpose_mul]
  calc S.transpose * R.transpose * (R*S) = S.transpose * (R.transpose*R) * S := by
        simp only [Matrix.mul_assoc]
       _ = 1 := by rw [hR, Matrix.mul_one, hS]

theorem conjugate_mul (R S G : Mat) :
    conjugate (R*S) G = conjugate R (conjugate S G) := by
  simp only [conjugate, Matrix.transpose_mul, Matrix.mul_assoc]

theorem reflection_twice (G : Mat) : conjugate reflection (conjugate reflection G) = G := by
  rw [← conjugate_mul, reflection_square]
  simp [conjugate]

theorem orthogonal_invariance_iff (Q : Quad) :
    OInvariant Q ↔ SOInvariant Q ∧ ∀ G, Q (conjugate reflection G) = Q G := by
  constructor
  · intro h
    exact ⟨fun R hR => h R hR.1, h reflection reflection_orthogonal⟩
  · rintro ⟨hSO, hF⟩ R hR G
    rcases orthogonal_det hR with hp | hn
    · exact hSO R ⟨hR,hp⟩ G
    · have hRF : Proper (R*reflection) :=
        ⟨orthogonal_mul hR reflection_orthogonal, by rw [Matrix.det_mul, hn, reflection_det]; norm_num⟩
      have he := hSO (R*reflection) hRF (conjugate reflection G)
      rw [conjugate_mul, reflection_twice, hF] at he
      exact he

end
end S11Invariants
