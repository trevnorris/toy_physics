import S11Invariants.Classification

/-! Dimensions of the actual invariant submodules, obtained from explicit
linear equivalences, plus the unique even/odd decomposition. -/
namespace S11Invariants
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
  zero_mem' := ⟨soSpace.zero_mem, by intro G; simp⟩
  add_mem' := by
    intro Q P hQ hP
    refine ⟨soSpace.add_mem hQ.1 hP.1, fun G => ?_⟩
    change Q (conjugate reflection G) + P (conjugate reflection G) = -(Q G + P G)
    rw [hQ.2 G, hP.2 G, neg_add]
  smul_mem' := by
    intro s Q hQ
    refine ⟨soSpace.smul_mem s hQ.1, fun G => ?_⟩
    change s * Q (conjugate reflection G) = -(s * Q G)
    rw [hQ.2 G, mul_neg]

def invariantMap : (Fin 4 → ℝ) →ₗ[ℝ] Quad where
  toFun := invariantForm
  map_add' v w := by simp [invariantForm, add_smul]; abel
  map_smul' s v := by simp [invariantForm, smul_add, smul_smul]

def soMap : (Fin 4 → ℝ) →ₗ[ℝ] soSpace :=
  invariantMap.codRestrict soSpace invariantForm_SO

theorem soMap_bijective : Function.Bijective soMap := by
  constructor
  · intro v w h
    exact invariantForm_injective (congrArg Subtype.val h)
  · intro Q
    obtain ⟨v,hv⟩ := (SO_classification Q.val).mp Q.property
    exact ⟨v, Subtype.ext hv.symm⟩

def soEquiv : (Fin 4 → ℝ) ≃ₗ[ℝ] soSpace := LinearEquiv.ofBijective soMap soMap_bijective

theorem so_dimension : Module.finrank ℝ soSpace = 4 := by
  rw [← soEquiv.finrank_eq]
  simp

def evenEmbedding : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 4 → ℝ) where
  toFun v := ![v 0,v 1,v 2,0]
  map_add' v w := by ext i; fin_cases i <;> simp
  map_smul' s v := by ext i; fin_cases i <;> simp

def oMap : (Fin 3 → ℝ) →ₗ[ℝ] oSpace :=
  (invariantMap.comp evenEmbedding).codRestrict oSpace (fun v =>
    (invariantForm_O _).mpr (by simp [evenEmbedding, Matrix.cons_val_three]))

theorem oMap_bijective : Function.Bijective oMap := by
  constructor
  · intro v w h
    have hv := invariantForm_injective (congrArg Subtype.val h)
    ext i
    fin_cases i
    · exact congrFun hv 0
    · exact congrFun hv 1
    · exact congrFun hv 2
  · intro Q
    obtain ⟨v,hv⟩ := (O_classification Q.val).mp Q.property
    exact ⟨v, Subtype.ext hv.symm⟩

def oEquiv : (Fin 3 → ℝ) ≃ₗ[ℝ] oSpace := LinearEquiv.ofBijective oMap oMap_bijective

theorem o_dimension : Module.finrank ℝ oSpace = 3 := by
  rw [← oEquiv.finrank_eq]
  simp

def oddEmbedding : ℝ →ₗ[ℝ] (Fin 4 → ℝ) where
  toFun d := ![0,0,0,d]
  map_add' v w := by ext i; fin_cases i <;> simp
  map_smul' s v := by ext i; fin_cases i <;> simp

def oddMap : ℝ →ₗ[ℝ] oddSpace :=
  (invariantMap.comp oddEmbedding).codRestrict oddSpace (fun d =>
    (odd_classification _).mpr ⟨d,rfl⟩)

theorem oddMap_bijective : Function.Bijective oddMap := by
  constructor
  · intro v w h
    have hv := invariantForm_injective (congrArg Subtype.val h)
    exact congrFun hv 3
  · intro Q
    obtain ⟨v,hv⟩ := (odd_classification Q.val).mp Q.property
    exact ⟨v, Subtype.ext hv.symm⟩

def oddEquiv : ℝ ≃ₗ[ℝ] oddSpace := LinearEquiv.ofBijective oddMap oddMap_bijective

theorem odd_dimension : Module.finrank ℝ oddSpace = 1 := by
  rw [← oddEquiv.finrank_eq]
  simp

/-- The even and odd invariant spaces have zero intersection. -/
theorem even_odd_disjoint : Disjoint oSpace oddSpace := by
  rw [Submodule.disjoint_def]
  intro Q he ho
  ext G
  have hp := he reflection reflection_orthogonal G
  have hn := ho.2 G
  change Q G = 0
  linarith

/-- Every SO-invariant form splits into its O-invariant and reflection-odd parts.
Together with disjointness, this gives a unique decomposition. -/
theorem even_odd_span : oSpace ⊔ oddSpace = soSpace := by
  apply le_antisymm
  · apply sup_le
    · intro Q hQ
      exact ((orthogonal_invariance_iff Q).mp hQ).1
    · intro Q hQ
      exact hQ.1
  · intro Q hQ
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp hQ
    apply Submodule.mem_sup.mpr
    refine ⟨invariantForm ![v 0,v 1,v 2,0], ?_, invariantForm ![0,0,0,v 3], ?_, ?_⟩
    · exact (invariantForm_O _).mpr (by simp [Matrix.cons_val_three])
    · exact (odd_classification _).mpr ⟨v 3,rfl⟩
    · change invariantMap ![v 0,v 1,v 2,0] + invariantMap ![0,0,0,v 3] = invariantMap v
      rw [← invariantMap.map_add]
      congr 1
      ext i
      fin_cases i <;> simp

end
end S11Invariants
