import S11OddDynamics.Variation

/-! E2/E3: non-null first variation and exhaustive longitudinal/transverse
mixing of the D2 odd contribution. No full XFORM_EXTRA spectrum is asserted. -/
namespace S11OddDynamics
noncomputable section
open S10Pilot

theorem dot_turn (k : Vec 2) : dot k (turn k) = 0 := by
  simp [dot, turn, Fin.sum_univ_two]
  ring

theorem turn_dot (k : Vec 2) : dot (turn k) k = 0 := by
  rw [dot_comm, dot_turn]

theorem normSq_turn (k : Vec 2) : normSq (turn k) = normSq k := by
  simp [normSq, dot, turn, Fin.sum_univ_two]
  ring

theorem normSq_eq_zero_iff (k : Vec 2) : normSq k = 0 ↔ k = 0 :=
  dotProduct_self_eq_zero

/-- These two directions cover every amplitude; at k=0 this identity degenerates. -/
theorem polarization_reconstruction (k a : Vec 2) :
    normSq k • a = dot k a • k + dot (turn k) a • turn k := by
  ext i
  fin_cases i <;> simp [normSq, dot, turn, Fin.sum_univ_two] <;> ring

theorem polarization_decomposition {k : Vec 2} (hk : k ≠ 0) (a : Vec 2) :
    a = (dot k a / normSq k) • k + (dot (turn k) a / normSq k) • turn k := by
  have hn := ne_of_gt (dot_self_pos hk)
  have h := congrArg (fun v : Vec 2 => (normSq k)⁻¹ • v) (polarization_reconstruction k a)
  simpa [smul_add, smul_smul, hn, div_eq_mul_inv, mul_comm] using h

theorem modalOperator_longitudinal (beta : ℝ) (k : Vec 2) :
    modalOperator beta k k = (-beta / 2 * normSq k) • turn k := by
  simp [modalOperator, turn_dot, normSq, smul_smul]

theorem modalOperator_transverse (beta : ℝ) (k : Vec 2) :
    modalOperator beta k (turn k) = (-beta / 2 * normSq k) • k := by
  have hn : dot (turn k) (turn k) = normSq k := normSq_turn k
  simp [modalOperator, dot_turn, hn, smul_smul]

theorem cross_pairing (beta : ℝ) (k : Vec 2) :
    dot (turn k) (modalOperator beta k k) = -beta / 2 * normSq k ^ 2 ∧
    dot k (modalOperator beta k (turn k)) = -beta / 2 * normSq k ^ 2 := by
  rw [modalOperator_longitudinal, modalOperator_transverse, dot_smul_right, dot_smul_right]
  change _ * normSq (turn k) = _ ∧ _ * normSq k = _
  rw [normSq_turn]
  constructor <;> ring

/-- Both cross-sector matrix elements are nonzero, with no generic-k assumption. -/
def Mixes (beta : ℝ) (k : Vec 2) : Prop :=
  dot (turn k) (modalOperator beta k k) ≠ 0 ∧
  dot k (modalOperator beta k (turn k)) ≠ 0

theorem mixing_iff (beta : ℝ) (k : Vec 2) :
    Mixes beta k ↔ beta ≠ 0 ∧ k ≠ 0 := by
  unfold Mixes
  rw [(cross_pairing beta k).1, (cross_pairing beta k).2]
  simp [normSq_eq_zero_iff]

theorem modalOperator_zero_beta (k a : Vec 2) : modalOperator 0 k a = 0 := by
  simp [modalOperator]

theorem modalOperator_zero_wavevector (beta : ℝ) (a : Vec 2) :
    modalOperator beta 0 a = 0 := by
  simp [modalOperator, turn]

/-- Exhaustive, disjoint parameter/wavevector cases, including the zero intersection. -/
theorem mixing_cases (beta : ℝ) (k : Vec 2) :
    (beta = 0 ∧ ∀ a, modalOperator beta k a = 0) ∨
    (beta ≠ 0 ∧ k = 0 ∧ ∀ a, modalOperator beta k a = 0) ∨
    (beta ≠ 0 ∧ k ≠ 0 ∧ Mixes beta k) := by
  by_cases hb : beta = 0
  · exact Or.inl ⟨hb, by subst beta; exact modalOperator_zero_beta k⟩
  · by_cases hk : k = 0
    · exact Or.inr (Or.inl ⟨hb, hk, by subst k; exact modalOperator_zero_wavevector beta⟩)
    · exact Or.inr (Or.inr ⟨hb, hk, (mixing_iff beta k).2 ⟨hb,hk⟩⟩)

/-- A smooth static wave is enough to witness a nonzero bulk variational derivative. -/
def witnessField : Point 2 → Vec 2 := planeWave 0 ![1,0] ![1,0]

theorem witness_smooth : SmoothField witnessField := smooth_planeWave _ _ _

theorem witness_eulerLagrange (beta : ℝ) :
    eulerLagrange beta witnessField 0 = ![0,-beta / 2] := by
  rw [witnessField, eulerLagrange_planeWave]
  ext i
  fin_cases i <;> simp [phase, dot, modalOperator, turn, Fin.sum_univ_two]

theorem witness_not_stationary {beta : ℝ} (hb : beta ≠ 0) :
    ¬ ActionStationary beta witnessField := by
  rw [actionStationary_iff_eulerLagrange witness_smooth]
  intro h
  have he := congrFun (h 0) 1
  rw [witness_eulerLagrange] at he
  simp only [Matrix.cons_val, Pi.zero_apply] at he
  exact hb (by linarith)

/-- Nonzero first variation against an allowed compact test field, not merely
nonzero evaluation of a supplied residual. Existence follows from the proved
integral/EL equivalence and the smooth witness. -/
theorem exists_nonzero_firstVariation {beta : ℝ} (hb : beta ≠ 0) :
    ∃ h : Point 2 → Vec 2, TestField h ∧ deriv (relativeAction beta witnessField h) 0 ≠ 0 := by
  by_contra hn
  apply witness_not_stationary hb
  intro h hh
  by_contra hf
  exact hn ⟨h, hh, hf⟩

def VariationallyNull (beta : ℝ) : Prop :=
  ∀ u : Point 2 → Vec 2, SmoothField u → ActionStationary beta u

theorem zero_beta_stationary (u : Point 2 → Vec 2) : ActionStationary 0 u := by
  intro h _
  have he : relativeAction 0 u h = fun _ : ℝ => 0 := by
    funext s
    simp [relativeAction, lagrangian]
  rw [he]
  simp

/-- The all-fields nullness claim has exactly one coefficient locus. -/
theorem variationallyNull_iff (beta : ℝ) : VariationallyNull beta ↔ beta = 0 := by
  constructor
  · intro h
    by_contra hb
    exact witness_not_stationary hb (h witnessField witness_smooth)
  · intro hb
    subst beta
    exact fun u _ => zero_beta_stationary u

end
end S11OddDynamics
