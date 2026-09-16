import S11D3Invariants.Classification

/-! Dimensions of the actual invariant submodules. -/
namespace S11D3Invariants
noncomputable section

def soSpace : Submodule ℝ Quad where
  carrier := {Q | SOInvariant Q}
  zero_mem' := by intro R _ G; rfl
  add_mem' := by
    intro Q P hQ hP R hR G
    change Q (conjugate R G) + P (conjugate R G) = Q G + P G
    rw [hQ R hR G, hP R hR G]
  smul_mem' := by
    intro s Q hQ R hR G
    change s * Q (conjugate R G) = s * Q G
    rw [hQ R hR G]

def oSpace : Submodule ℝ Quad where
  carrier := {Q | OInvariant Q}
  zero_mem' := by intro R _ G; rfl
  add_mem' := by
    intro Q P hQ hP R hR G
    change Q (conjugate R G) + P (conjugate R G) = Q G + P G
    rw [hQ R hR G, hP R hR G]
  smul_mem' := by
    intro s Q hQ R hR G
    change s * Q (conjugate R G) = s * Q G
    rw [hQ R hR G]

def oddSpace : Submodule ℝ Quad where
  carrier := {Q | SOInvariant Q ∧ ReflectionOdd Q}
  zero_mem' := (odd_classification 0).mpr rfl
  add_mem' := by
    intro Q P hQ hP
    rw [(odd_classification Q).mp hQ, (odd_classification P).mp hP, add_zero]
    exact (odd_classification 0).mpr rfl
  smul_mem' := by
    intro s Q hQ
    rw [(odd_classification Q).mp hQ, smul_zero]
    exact (odd_classification 0).mpr rfl

theorem so_eq_o : soSpace = oSpace := by
  ext Q
  exact SO_iff_O Q

theorem odd_eq_bot : oddSpace = ⊥ := by
  ext Q
  exact odd_classification Q

def invariantMap : (Fin 3 → ℝ) →ₗ[ℝ] Quad where
  toFun := invariantForm
  map_add' v w := by simp [invariantForm, add_smul]; abel
  map_smul' s v := by simp [invariantForm, smul_add, smul_smul]

def soMap : (Fin 3 → ℝ) →ₗ[ℝ] soSpace :=
  invariantMap.codRestrict soSpace invariantForm_SO

theorem soMap_bijective : Function.Bijective soMap := by
  constructor
  · intro v w h
    exact invariantForm_injective (congrArg Subtype.val h)
  · intro Q
    obtain ⟨v,hv⟩ := (SO_classification Q.val).mp Q.property
    exact ⟨v, Subtype.ext hv.symm⟩

def soEquiv : (Fin 3 → ℝ) ≃ₗ[ℝ] soSpace := LinearEquiv.ofBijective soMap soMap_bijective

theorem so_dimension : Module.finrank ℝ soSpace = 3 := by
  rw [← soEquiv.finrank_eq]
  simp

theorem o_dimension : Module.finrank ℝ oSpace = 3 := by
  rw [← so_eq_o, so_dimension]

theorem odd_dimension : Module.finrank ℝ oddSpace = 0 := by
  rw [odd_eq_bot]
  simp

theorem census : (Module.finrank ℝ soSpace, Module.finrank ℝ oSpace,
    Module.finrank ℝ oddSpace) = (3,3,0) := by
  rw [so_dimension, o_dimension, odd_dimension]

end
end S11D3Invariants
