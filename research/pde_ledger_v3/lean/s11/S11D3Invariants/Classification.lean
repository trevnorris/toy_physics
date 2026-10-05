import S11D3Invariants.Constraints

/-! Exact full-group classification and absence of a reflection-odd extra.
These are density identities, with no quotient by total divergences. -/
namespace S11D3Invariants
noncomputable section

theorem SO_classification (Q : Quad) :
    SOInvariant Q ↔ ∃ v : Fin 3 → ℝ, Q = invariantForm v := by
  constructor
  · intro hQ
    obtain ⟨c,hc⟩ := quadratic_representation Q
    refine ⟨![c 4/2,c 11/2,c 9], ?_⟩
    ext G
    simpa only [invariantForm_apply, Matrix.cons_val] using invariant_polynomial hQ c hc G
  · rintro ⟨v,rfl⟩
    exact invariantForm_SO v

theorem SO_iff_O (Q : Quad) : SOInvariant Q ↔ OInvariant Q := by
  constructor
  · intro hQ
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp hQ
    exact invariantForm_O v
  · intro hQ R hR
    exact hQ R hR.1

theorem O_classification (Q : Quad) :
    OInvariant Q ↔ ∃ v : Fin 3 → ℝ, Q = invariantForm v := by
  rw [← SO_iff_O, SO_classification]

theorem invariantForm_injective : Function.Injective invariantForm := by
  intro v w h
  have hF := congrArg (fun Q : Quad => Q (Matrix.of ![![0,1,0],![0,0,0],![0,0,0]])) h
  have hD := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0],![0,0,0],![0,0,0]])) h
  have hT := congrArg (fun Q : Quad => Q (Matrix.of ![![1,0,0],![0,1,0],![0,0,0]])) h
  simp only [invariantForm_apply, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply,
    Matrix.transpose_apply, Fin.sum_univ_three] at hF hD hT
  norm_num [Matrix.cons_val_two] at hF hD hT
  ext i
  fin_cases i <;> dsimp <;> linarith

theorem SO_unique (Q : Quad) (hQ : SOInvariant Q) :
    ∃! v : Fin 3 → ℝ, Q = invariantForm v := by
  obtain ⟨v,hv⟩ := (SO_classification Q).mp hQ
  exact ⟨v,hv,fun w hw => invariantForm_injective (hw.symm.trans hv)⟩

theorem odd_classification (Q : Quad) :
    SOInvariant Q ∧ ReflectionOdd Q ↔ Q = 0 := by
  constructor
  · rintro ⟨hSO,hodd⟩
    ext G
    have he := (SO_iff_O Q).mp hSO reflection reflection_orthogonal G
    have ho := hodd G
    change Q G = 0
    linarith
  · rintro rfl
    constructor
    · intro R _ G; rfl
    · intro G; simp

end
end S11D3Invariants
