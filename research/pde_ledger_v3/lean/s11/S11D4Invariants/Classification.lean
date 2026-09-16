import S11D4Invariants.Constraints

/-! Exact full-group density classification and its reflection-even/odd split.
No Euler–Lagrange or integration-by-parts quotient is taken. -/
namespace S11D4Invariants
noncomputable section

theorem SO_classification (Q : Quad) :
    SOInvariant Q ↔ ∃ v : Fin 4 → ℝ, Q = invariantForm v := by
  constructor
  · intro hQ
    obtain ⟨c,hc⟩ := quadratic_representation Q
    refine ⟨![c 5/2,c 19/2,c 16,c 26], ?_⟩
    ext G
    simpa only [invariantForm_apply, Matrix.cons_val] using invariant_polynomial hQ c hc G
  · rintro ⟨v,rfl⟩
    exact invariantForm_SO v

theorem invariantForm_injective : Function.Injective invariantForm := by
  intro v w h
  have hF := congrArg (fun Q : Quad => Q (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]])) h
  have hD := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]])) h
  have hT := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0,0],![0,1,0,0],![0,0,0,0],![0,0,0,0]])) h
  have hP := congrArg (fun Q : Quad => Q (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,1],![0,0,0,0]])) h
  simp only [invariantForm_apply, orientation, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
    Matrix.transpose_apply, Fin.sum_univ_four] at hF hD hT hP
  norm_num [Matrix.cons_val_two, Matrix.cons_val_three] at hF hD hT hP
  ext i
  fin_cases i <;> dsimp <;> linarith

theorem SO_unique (Q : Quad) (hQ : SOInvariant Q) :
    ∃! v : Fin 4 → ℝ, Q = invariantForm v := by
  obtain ⟨v,hv⟩ := (SO_classification Q).mp hQ
  exact ⟨v,hv,fun w hw => invariantForm_injective (hw.symm.trans hv)⟩

theorem invariantForm_O (v : Fin 4 → ℝ) : OInvariant (invariantForm v) ↔ v 3 = 0 := by
  constructor
  · intro hQ
    have h := hQ reflection reflection_orthogonal
      (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,1],![0,0,0,0]])
    rw [invariantForm_reflection, invariantForm_apply] at h
    simp only [orientation, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] at h
    norm_num [Matrix.cons_val_two, Matrix.cons_val_three] at h
    linarith
  · intro h R hR G
    simp only [invariantForm_apply, h, zero_mul, add_zero, trace_conjugate hR,
      conjugate_transpose, conjugate_mul hR]

theorem invariantForm_odd (v : Fin 4 → ℝ) :
    ReflectionOdd (invariantForm v) ↔ v 0 = 0 ∧ v 1 = 0 ∧ v 2 = 0 := by
  constructor
  · intro hQ
    have hF := hQ (Matrix.of ![![0,1,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]])
    have hD := hQ (Matrix.of ![![1,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]])
    have hT := hQ (Matrix.of ![![1,0,0,0],![0,1,0,0],![0,0,0,0],![0,0,0,0]])
    rw [invariantForm_reflection, invariantForm_apply] at hF hD hT
    simp only [orientation, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] at hF hD hT
    norm_num [Matrix.cons_val_two, Matrix.cons_val_three] at hF hD hT
    exact ⟨by linarith, by linarith, by linarith⟩
  · rintro ⟨h0,h1,h2⟩ G
    rw [invariantForm_reflection, invariantForm_apply]
    simp only [h0,h1,h2,zero_mul,zero_add,zero_sub]

theorem O_classification (Q : Quad) :
    OInvariant Q ↔ ∃ v : Fin 3 → ℝ, Q = invariantForm ![v 0,v 1,v 2,0] := by
  constructor
  · intro hQ
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp (fun R hR => hQ R hR.1)
    have hv := (invariantForm_O v).mp hQ
    refine ⟨![v 0,v 1,v 2], ?_⟩
    congr 1
    ext i
    fin_cases i <;> simp [Matrix.cons_val_two, hv]
  · rintro ⟨v,rfl⟩
    exact (invariantForm_O _).mpr (by simp [Matrix.cons_val_three])

theorem odd_classification (Q : Quad) :
    SOInvariant Q ∧ ReflectionOdd Q ↔ ∃ d : ℝ, Q = invariantForm ![0,0,0,d] := by
  constructor
  · rintro ⟨hSO,hodd⟩
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp hSO
    obtain ⟨h0,h1,h2⟩ := (invariantForm_odd v).mp hodd
    refine ⟨v 3, ?_⟩
    congr 1
    ext i
    fin_cases i <;> simp [h0,h1,h2]
  · rintro ⟨d,rfl⟩
    exact ⟨invariantForm_SO _, (invariantForm_odd _).mpr (by norm_num [Matrix.cons_val_two])⟩

end
end S11D4Invariants
