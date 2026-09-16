import S11Invariants.Census

/-! Admissible examples and the concrete native-orientation counterexample. -/
namespace S11Invariants
noncomputable section

def oddPairing : Quad := invariantForm ![0,0,0,1]

theorem oddPairing_SO : SOInvariant oddPairing := invariantForm_SO _

theorem oddPairing_odd : ReflectionOdd oddPairing :=
  (invariantForm_odd _).mpr (by norm_num [Matrix.cons_val_two])

theorem oddPairing_not_O : ¬ OInvariant oddPairing := by
  rw [oddPairing, invariantForm_O]
  norm_num [Matrix.cons_val_three]

/-- Omitting the odd pairing loses an actual SO invariant. -/
theorem oddPairing_not_even_span :
    ¬ ∃ v : Fin 3 → ℝ, oddPairing = invariantForm ![v 0,v 1,v 2,0] := by
  rw [← O_classification]
  exact oddPairing_not_O

def wrongNativeForm : Quad :=
  polynomialForm ![1/2,0,0,0,-1/4,0,0,1/2,0,1/4]

theorem wrongNativeForm_apply (G : Mat) :
    wrongNativeForm G = G 0 0 ^ 2 + G 0 1 * G 1 0 + G 1 1 ^ 2 := by
  simp only [wrongNativeForm, polynomialForm_apply, polynomial, coordinates, Matrix.cons_val]
  ring

theorem wrongNativeForm_not_SO : ¬ SOInvariant wrongNativeForm := by
  intro h
  have hc := invariant_coefficients h ![1/2,0,0,0,-1/4,0,0,1/2,0,1/4]
    (fun G => polynomialForm_apply _ G)
  have he := hc.2.2.2.2.1
  change (1/2 : ℝ) = 1/4 at he
  norm_num at he

/-- The same rational proper rotation and matrix used in the native probe. -/
theorem wrongNativeForm_witness :
    Proper (rotation (3/5) (4/5)) ∧
    wrongNativeForm !![1,0;0,0] = 1 ∧
    wrongNativeForm (conjugate (rotation (3/5) (4/5)) !![1,0;0,0]) = 481/625 := by
  refine ⟨rotation_proper (by norm_num), ?_, ?_⟩
  · rw [wrongNativeForm_apply]
    norm_num
  · change polynomialForm ![1/2,0,0,0,-1/4,0,0,1/2,0,1/4]
      (conjugate (rotation (3/5) (4/5)) !![1,0;0,0]) = 481/625
    rw [polynomialForm_apply, coordinates_rotation]
    have hg : coordinates (!![1,0;0,0] : Mat) = ![1,0,1,0] := by
      ext i
      fin_cases i <;> norm_num [coordinates]
    rw [hg]
    simp only [polynomial, Matrix.cons_val]
    norm_num

end
end S11Invariants
