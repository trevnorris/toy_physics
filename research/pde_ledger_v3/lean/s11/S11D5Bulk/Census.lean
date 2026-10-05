import S11D5Bulk.Bulk

/-! D5B.2: exhaustive equality and nullness of the already classified D5 family. -/
namespace S11D5Bulk
noncomputable section
open S10Pilot

def BulkEquivalent (v w : Coeff) : Prop :=
  ∀ u : Point 5 → Vec 5, SmoothField u → ∀ x, eulerLagrange v u x = eulerLagrange w u x

theorem bulkEquivalent_iff (v w : Coeff) :
    BulkEquivalent v w ↔ v 2 = w 2 ∧ v 0 + v 1 = w 0 + w 1 := by
  constructor
  · intro h
    have ht := congrFun (h (planeWave 0 ![1,0,0,0,0] ![0,1,0,0,0])
      (smooth_planeWave _ _ _) 0) 1
    have hl := congrFun (h (planeWave 0 ![1,0,0,0,0] ![1,0,0,0,0])
      (smooth_planeWave _ _ _) 0) 0
    simp [eulerLagrange_planeWave, phase, modalOperator, normSq, dot,
      S11D5Invariants.sum_five, Matrix.cons_val] at ht hl
    constructor <;> linarith
  · rintro ⟨hc,hs⟩ u hu x
    rw [eulerLagrange_eq hu, eulerLagrange_eq hu, hc, hs]

def responseMap : Coeff →ₗ[ℝ] Vec 2 where
  toFun v := ![v 2, v 0 + v 1]
  map_add' v w := by
    ext i
    fin_cases i <;> simp
    ring
  map_smul' s v := by
    ext i
    fin_cases i <;> simp
    ring

theorem bulkEquivalent_response (v w : Coeff) :
    BulkEquivalent v w ↔ responseMap v = responseMap w := by
  rw [bulkEquivalent_iff]
  simp [responseMap, funext_iff, Fin.forall_fin_two]

theorem response_surjective : Function.Surjective responseMap := by
  intro s
  refine ⟨![s 1,0,s 0], ?_⟩
  ext i
  fin_cases i <;> simp [responseMap]

def VariationallyNull (v : Coeff) : Prop :=
  ∀ u : Point 5 → Vec 5, SmoothField u → ActionStationary v u

theorem variationallyNull_iff (v : Coeff) :
    VariationallyNull v ↔ v 2 = 0 ∧ v 0 + v 1 = 0 := by
  have he : VariationallyNull v ↔ BulkEquivalent v 0 := by
    constructor
    · intro h u hu x
      rw [(actionStationary_iff_eulerLagrange hu v).mp (h u hu) x, eulerLagrange_eq hu]
      simp
    · intro h u hu
      apply (actionStationary_iff_eulerLagrange hu v).mpr
      intro x
      rw [h u hu x, eulerLagrange_eq hu]
      simp
  rw [he, bulkEquivalent_iff]
  simp

theorem null_parameterization (v : Coeff) :
    VariationallyNull v ↔ ∃ t : ℝ, v = ![t,-t,0] := by
  rw [variationallyNull_iff]
  constructor
  · rintro ⟨hc,hs⟩
    refine ⟨v 0, ?_⟩
    ext i
    fin_cases i
    · rfl
    · change v 1 = -v 0
      linarith
    · exact hc
  · rintro ⟨t,rfl⟩
    simp

theorem nullSpace_identification (v : Coeff) :
    v ∈ LinearMap.ker responseMap ↔ VariationallyNull v := by
  rw [LinearMap.mem_ker, variationallyNull_iff]
  simp [responseMap, funext_iff, Fin.forall_fin_two]

theorem bulk_response_dimension : Module.finrank ℝ (LinearMap.range responseMap) = 2 := by
  rw [LinearMap.range_eq_top.mpr response_surjective]
  simp [Vec]

theorem null_dimension : Module.finrank ℝ (LinearMap.ker responseMap) = 1 := by
  have h := responseMap.finrank_range_add_finrank_ker
  rw [bulk_response_dimension] at h
  norm_num [Coeff, Vec] at h ⊢
  omega

theorem null_nonzero_example : VariationallyNull ![1,-1,0] ∧ (![1,-1,0] : Coeff) ≠ 0 := by
  constructor
  · rw [variationallyNull_iff]
    norm_num [Matrix.cons_val_two]
  · intro h
    have := congrFun h 0
    norm_num at this

theorem exists_nonzero_firstVariation (v : Coeff) (hv : ¬(v 2 = 0 ∧ v 0 + v 1 = 0)) :
    ∃ u h : Point 5 → Vec 5, SmoothField u ∧ TestField h ∧
      deriv (relativeAction v u h) 0 ≠ 0 := by
  have hn : ¬VariationallyNull v := fun h => hv ((variationallyNull_iff v).mp h)
  unfold VariationallyNull ActionStationary at hn
  push Not at hn
  obtain ⟨u,hu,h,hh,hn⟩ := hn
  exact ⟨u,h,hu,hh,hn⟩

end
end S11D5Bulk
