import S11D4Invariants.Census

/-! Nonvacuity, completeness and full-group sensitivity controls for D4. -/
namespace S11D4Invariants
noncomputable section

theorem trace_not_omittable :
    ¬ ∃ b c d : ℝ, invariantForm ![1,0,0,0] = invariantForm ![0,b,c,d] := by
  rintro ⟨b,c,d,h⟩
  have he := congrFun (invariantForm_injective h) 0
  norm_num at he

theorem traceOfSquare_not_omittable :
    ¬ ∃ a c d : ℝ, invariantForm ![0,1,0,0] = invariantForm ![a,0,c,d] := by
  rintro ⟨a,c,d,h⟩
  have he := congrFun (invariantForm_injective h) 1
  norm_num at he

theorem frobenius_not_omittable :
    ¬ ∃ a b d : ℝ, invariantForm ![0,0,1,0] = invariantForm ![a,b,0,d] := by
  rintro ⟨a,b,d,h⟩
  have he := congrFun (invariantForm_injective h) 2
  norm_num [Matrix.cons_val_two] at he

theorem orientation_not_omittable :
    ¬ ∃ a b c : ℝ, invariantForm ![0,0,0,1] = invariantForm ![a,b,c,0] := by
  rintro ⟨a,b,c,h⟩
  have he := congrFun (invariantForm_injective h) 3
  norm_num [Matrix.cons_val_three] at he

theorem orientation_nonzero : invariantForm ![0,0,0,1] ≠ 0 := by
  intro h
  have he := congrArg (fun Q : Quad => Q (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,1],![0,0,0,0]])) h
  norm_num [invariantForm_apply, orientation, Matrix.cons_val_two, Matrix.cons_val_three] at he

theorem nonzero_odd_exists : ∃ Q : Quad, SOInvariant Q ∧ ReflectionOdd Q ∧ Q ≠ 0 := by
  refine ⟨invariantForm ![0,0,0,1], invariantForm_SO _, ?_, orientation_nonzero⟩
  exact (invariantForm_odd _).mpr (by norm_num [Matrix.cons_val_two])

theorem zero_invariant : SOInvariant (0 : Quad) ∧ ReflectionOdd (0 : Quad) := by
  constructor
  · intro R _ G; rfl
  · intro G; simp

theorem single_entry_not_invariant : ¬ SOInvariant (monomial 0 0) := by
  intro hQ
  let G : Mat := ![![1,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]]
  have he := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) G
  change coordinates (conjugate (rotationXY 0 1) G) 0 *
    coordinates (conjugate (rotationXY 0 1) G) 0 = coordinates G 0 * coordinates G 0 at he
  rw [coordinates_rotationXY] at he
  norm_num [coordinates, G] at he

end
end S11D4Invariants
