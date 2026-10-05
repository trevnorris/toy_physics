import S11D5Invariants.Census

/-! Explicit nonemptiness and omission controls for D5.1–D5.2. -/
namespace S11D5Invariants
noncomputable section

theorem trace_not_omittable :
    ¬ ∃ b c : ℝ, invariantForm ![1,0,0] = invariantForm ![0,b,c] := by
  rintro ⟨b,c,h⟩
  have he := congrFun (invariantForm_injective h) 0
  norm_num at he

theorem traceOfSquare_not_omittable :
    ¬ ∃ a c : ℝ, invariantForm ![0,1,0] = invariantForm ![a,0,c] := by
  rintro ⟨a,c,h⟩
  have he := congrFun (invariantForm_injective h) 1
  norm_num at he

theorem frobenius_not_omittable :
    ¬ ∃ a b : ℝ, invariantForm ![0,0,1] = invariantForm ![a,b,0] := by
  rintro ⟨a,b,h⟩
  have he := congrFun (invariantForm_injective h) 2
  norm_num [Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] at he

theorem nonzero_invariant_exists : ∃ Q : Quad, SOInvariant Q ∧ Q ≠ 0 := by
  refine ⟨invariantForm ![1,0,0], invariantForm_SO _, ?_⟩
  intro h
  have he := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0,0,0],![0,0,0,0,0],![0,0,0,0,0],![0,0,0,0,0],![0,0,0,0,0]])) h
  simp only [invariantForm_apply, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
    Matrix.transpose_apply, sum_five] at he
  norm_num [Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four] at he

theorem zero_invariant : SOInvariant (0 : Quad) ∧ ReflectionOdd (0 : Quad) :=
  (odd_classification 0).mpr rfl

theorem nonzero_odd_impossible : ¬ ∃ Q : Quad, SOInvariant Q ∧ ReflectionOdd Q ∧ Q ≠ 0 := by
  rintro ⟨Q,hSO,hOdd,hne⟩
  exact hne ((odd_classification Q).mp ⟨hSO,hOdd⟩)

theorem single_entry_not_invariant : ¬ SOInvariant (monomial 0 0) := by
  intro hQ
  have he := ((SO_iff_O _).mp hQ) (rotationXY 0 1)
    (rotationXY_orthogonal (by norm_num)) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  change coordinates (conjugate (rotationXY 0 1) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])) 0 *
    coordinates (conjugate (rotationXY 0 1) (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])) 0 = coordinates (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 0 * coordinates (decode ![1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 0 at he
  rw [coordinates_rotationXY] at he
  simp only [coordinates, decode, Matrix.cons_val] at he
  norm_num at he

def offDiagonal : Mat := decode ![0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
def twoDiagonal : Mat := decode ![1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
def offSymmetric : Mat := decode ![0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

theorem trace_normalization : traceSquare twoDiagonal = 4 := by
  rw [traceSquare_apply]
  simp only [Matrix.trace, Matrix.diag_apply, sum_five]
  simp only [twoDiagonal, decode, Matrix.cons_val]
  norm_num

theorem pair_normalization : traceOfSquare offSymmetric = 2 := by
  rw [traceOfSquare_apply]
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, sum_five]
  simp only [offSymmetric, decode, Matrix.cons_val]
  norm_num

theorem frobenius_normalization : frobeniusSquare offDiagonal = 1 := by
  rw [frobeniusSquare_apply]
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  simp only [offDiagonal, decode, Matrix.cons_val]
  norm_num

end
end S11D5Invariants
