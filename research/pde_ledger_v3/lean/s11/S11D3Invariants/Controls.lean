import S11D3Invariants.Census

/-! Explicit nonemptiness and omission controls for J1-J2. -/
namespace S11D3Invariants
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
  norm_num [Matrix.cons_val_two] at he

theorem nonzero_invariant_exists : ∃ Q : Quad, SOInvariant Q ∧ Q ≠ 0 := by
  refine ⟨invariantForm ![1,0,0], invariantForm_SO _, ?_⟩
  intro h
  have he := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0],![0,0,0],![0,0,0]])) h
  simp only [invariantForm_apply, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
    Matrix.transpose_apply, Fin.sum_univ_three] at he
  norm_num [Matrix.cons_val_two] at he

theorem zero_invariant : SOInvariant (0 : Quad) ∧ ReflectionOdd (0 : Quad) :=
  (odd_classification 0).mpr rfl

theorem nonzero_odd_impossible : ¬ ∃ Q : Quad, SOInvariant Q ∧ ReflectionOdd Q ∧ Q ≠ 0 := by
  rintro ⟨Q,hSO,hOdd,hne⟩
  exact hne ((odd_classification Q).mp ⟨hSO,hOdd⟩)

theorem single_entry_not_invariant : ¬ SOInvariant (monomial 0 0) := by
  intro hQ
  have image : conjugate (rotationXY 0 1) (Matrix.of ![![1,0,0],![0,0,0],![0,0,0]]) =
      (Matrix.of ![![0,0,0],![0,1,0],![0,0,0]]) := by
    ext i j
    rw [conjugate_apply]
    fin_cases i <;> fin_cases j <;>
      norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
  have he := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num))
    (Matrix.of ![![1,0,0],![0,0,0],![0,0,0]])
  rw [image] at he
  norm_num [monomial, coordinateMap, coordinates, Matrix.cons_val] at he

end
end S11D3Invariants
