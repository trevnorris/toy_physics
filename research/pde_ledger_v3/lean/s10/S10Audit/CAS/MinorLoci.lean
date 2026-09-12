import S10Audit.CAS.MinorSupport

namespace S10Audit.CAS
open S10Pilot
noncomputable section

theorem wavevector_zero_iff (k : Vec 3) : k = 0 ↔ k 0 = 0 ∧ k 1 = 0 ∧ k 2 = 0 := by
  constructor
  · intro h; simp [h]
  · rintro ⟨h0, h1, h2⟩
    funext i
    fin_cases i
    · exact h0
    · exact h1
    · exact h2

theorem normSq_zero_iff (k : Vec 3) : normSq k = 0 ↔ k 0 = 0 ∧ k 1 = 0 ∧ k 2 = 0 := by
  rw [← wavevector_zero_iff]
  exact dotProduct_self_eq_zero

theorem perp_pair_zero_iff (k : Vec 3) : k 1 ^ 2 + k 2 ^ 2 = 0 ↔ k 1 = 0 ∧ k 2 = 0 := by
  constructor
  · intro h
    constructor <;> nlinarith [sq_nonneg (k 1), sq_nonneg (k 2)]
  · rintro ⟨h1, h2⟩; simp [h1, h2]

open Lean Elab Tactic in
elab "minor_zero" : tactic => do
  evalTactic (← `(tactic| norm_num [staticRankMinors, staticTransverseMinors,
    ordinaryRankMinors, ordinaryTransverseMinors, extraRankMinors, extraTransverseMinors,
    Fin.forall_fin_succ, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three,
    Matrix.cons_val_four, Matrix.head_cons, Matrix.tail_cons, mul_eq_zero, div_eq_zero_iff, sub_eq_zero,
    normSq_zero_iff, perp_pair_zero_iff, *]))
  if !(← getGoals).isEmpty then evalTactic (← `(tactic| tauto))

theorem staticRank_locus (mu : ℝ) (k : Vec 3) (hm : mu ≠ 0) :
    (∀ i, staticRankMinors mu k i = 0) ↔ k = 0 := by
  rw [wavevector_zero_iff]
  minor_zero

theorem staticTransverse_locus (mu : ℝ) (k : Vec 3) (hm : mu ≠ 0) :
    (∀ i, staticTransverseMinors mu k i = 0) ↔ k = 0 := by
  rw [wavevector_zero_iff]
  minor_zero

theorem ordinaryRank_locus (mu sigma : ℝ) (k : Vec 3) (hm : mu ≠ 0) (hs1 : sigma ≠ 1) :
    (∀ i, ordinaryRankMinors mu sigma k i = 0) ↔ k 1 = 0 ∧ k 2 = 0 := by
  minor_zero

theorem ordinaryTransverse_locus (mu sigma : ℝ) (k : Vec 3) (hm : mu ≠ 0) (hs1 : sigma ≠ 1) :
    (∀ i, ordinaryTransverseMinors mu sigma k i = 0) ↔ k 1 = 0 ∧ k 2 = 0 := by
  minor_zero

theorem extraRank_locus (mu sigma : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    (∀ i, extraRankMinors mu sigma k i = 0) ↔ k 1 = 0 ∧ k 2 = 0 := by
  minor_zero

theorem extraTransverse_locus (mu sigma : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    (∀ i, extraTransverseMinors mu sigma k i = 0) ↔ k 0 = 0 ∨ (k 1 = 0 ∧ k 2 = 0) := by
  minor_zero

theorem exceptional_strata_disjoint (k : Vec 3) (hk : k ≠ 0) :
    ¬ (k 0 = 0 ∧ k 1 = 0 ∧ k 2 = 0) := by
  rwa [← wavevector_zero_iff]

theorem extraTransverse_nonzero_partition (mu sigma : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) (hk : k ≠ 0) :
    (∀ i, extraTransverseMinors mu sigma k i = 0) ↔
      (k 0 = 0 ∧ (k 1 ≠ 0 ∨ k 2 ≠ 0)) ∨ (k 0 ≠ 0 ∧ k 1 = 0 ∧ k 2 = 0) := by
  rw [extraTransverse_locus mu sigma k hm hs hs1]
  have hn := exceptional_strata_disjoint k hk
  tauto

theorem perpendicular_drops_only_transverse (mu sigma : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) (hk : k ≠ 0) (h0 : k 0 = 0) :
    (¬ ∀ i, extraRankMinors mu sigma k i = 0) ∧
      (∀ i, extraTransverseMinors mu sigma k i = 0) := by
  rw [extraRank_locus mu sigma k hm hs hs1, extraTransverse_locus mu sigma k hm hs hs1]
  have hn := exceptional_strata_disjoint k hk
  tauto

theorem parallel_drops_both (mu sigma : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) (h1 : k 1 = 0) (h2 : k 2 = 0) :
    (∀ i, extraRankMinors mu sigma k i = 0) ∧
      (∀ i, extraTransverseMinors mu sigma k i = 0) := by
  rw [extraRank_locus mu sigma k hm hs hs1, extraTransverse_locus mu sigma k hm hs hs1]
  exact ⟨⟨h1, h2⟩, Or.inr ⟨h1, h2⟩⟩

theorem static_no_nonzero_rank_drop (mu : ℝ) (k : Vec 3) (hm : mu ≠ 0) (hk : k ≠ 0) :
    (¬ ∀ i, staticRankMinors mu k i = 0) ∧ (¬ ∀ i, staticTransverseMinors mu k i = 0) := by
  rw [staticRank_locus mu k hm, staticTransverse_locus mu k hm]
  exact ⟨hk, hk⟩

theorem positive_or_negative (x : ℝ) : (0 < x ∨ x < 0) ↔ x ≠ 0 := by
  constructor
  · rintro (h | h)
    · exact ne_of_gt h
    · exact ne_of_lt h
  · intro h
    exact (lt_or_gt_of_ne h).symm

open Lean Elab Tactic in
elab "locus_eval" h:Lean.Parser.Tactic.simpLemma : tactic => do
  evalTactic (← `(tactic| simp only [$h, wavevector_zero_iff, positive_or_negative]))
  if !(← getGoals).isEmpty then evalTactic (← `(tactic| tauto))

end
end S10Audit.CAS
