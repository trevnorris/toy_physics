import S10Audit.CAS.RootBindings

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section
set_option backward.isDefEq.respectTransparency false

def coincidenceDomain (rho mu sigma lambda : ℝ) : Prop :=
  0 < rho ∧ 0 < mu ∧ 0 < sigma ∧ sigma ≠ 1 ∧ 0 < lambda

def pairLeft : Fin 3 → Fin 3 := ![0, 0, 1]
def pairRight : Fin 3 → Fin 3 := ![1, 2, 2]

def coincidenceGeometry (i : Fin 3) (k : Vec 3) : Prop :=
  if i = 2 then k 1 = 0 ∧ k 2 = 0 else k = 0

def allowedCoincidence (i : Fin 3) (k : Vec 3) : Prop :=
  coincidenceGeometry i k ∧ k ≠ 0

theorem normSq_positive_iff (k : Vec 3) : 0 < normSq k ↔ k ≠ 0 := by
  constructor
  · intro h he
    subst k
    simp [normSq, dot] at h
  · exact dot_self_pos

theorem coincidence_pair_geometry (rho mu sigma : ℝ) (k : Vec 3) (i : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    referenceRoot rho mu sigma k (pairLeft i) = referenceRoot rho mu sigma k (pairRight i) ↔
      coincidenceGeometry i k := by
  fin_cases i
  · change 0 = referenceRoot rho mu sigma k 1 ↔ k = 0
    constructor
    · intro h
      by_contra hk
      exact (root_nonzero rho mu sigma k hr hm hs hk).1 h.symm
    · intro h
      subst k
      simp [referenceRoot, coneValue, normSq, dot]
  · change 0 = referenceRoot rho mu sigma k 2 ↔ k = 0
    constructor
    · intro h
      by_contra hk
      exact (root_nonzero rho mu sigma k hr hm hs hk).2 h.symm
    · intro h
      subst k
      simp [referenceRoot, extraConeValue, extraValue, extraNumerator, perpSq, normSq, dot]
  · change referenceRoot rho mu sigma k 1 = referenceRoot rho mu sigma k 2 ↔ k 1 = 0 ∧ k 2 = 0
    rw [eq_comm, root_coincidence rho mu sigma k hr hm (ne_of_gt hs) hs1]
    have hq : perpSq 0 k = k 1 ^ 2 + k 2 ^ 2 := by
      simp only [perpSq, normSq, dot, Fin.sum_univ_three]
      ring
    rw [hq, perp_pair_zero_iff]

theorem allowedCoincidence_geometry (i : Fin 3) (k : Vec 3) :
    allowedCoincidence i k ↔ (if i = 2 then k 0 ≠ 0 ∧ k 1 = 0 ∧ k 2 = 0 else False) := by
  unfold allowedCoincidence coincidenceGeometry
  by_cases hi : i = 2
  · simp only [hi, if_true, ne_eq, wavevector_zero_iff]
    tauto
  · simp [hi]

theorem allowedCoincidence_exists (i : Fin 3) :
    (∃ k : Vec 3, allowedCoincidence i k) ↔ i = 2 := by
  constructor
  · rintro ⟨k, hk⟩
    by_contra h
    simp [allowedCoincidence_geometry, h] at hk
  · intro h
    subst i
    refine ⟨![1, 0, 0], ?_⟩
    rw [allowedCoincidence_geometry]
    norm_num [Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

theorem positive_sigma_split (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (0 < sigma ∧ sigma < 1) ∨ 1 < sigma := by
  rcases lt_or_gt_of_ne hs1 with h | h
  · exact Or.inl ⟨hs, h⟩
  · exact Or.inr h

theorem sumSquares_positive_iff (k : Vec 3) :
    0 < k 0 ^ 2 + k 1 ^ 2 + k 2 ^ 2 ↔ k ≠ 0 := by
  have he : k 0 ^ 2 + k 1 ^ 2 + k 2 ^ 2 = normSq k := by
    simp only [normSq, dot, Fin.sum_univ_three]
    ring
  rw [he, normSq_positive_iff]

open Lean Elab Tactic in
elab "coincidence_eval" h:Lean.Parser.Tactic.simpLemma : tactic => do
  evalTactic (← `(tactic| simp only [$h, sumSquares_positive_iff, allowedCoincidence, coincidenceGeometry,
    PYLoci.parallelPoint, PYLoci.perpendicularPoint, WLLoci.parallelPoint,
    WLLoci.perpendicularPoint, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]))
  if !(← getGoals).isEmpty then evalTactic (← `(tactic| norm_num [Fin.ext_iff, Fin.coe_ofNat_eq_mod]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| simp_all only [ne_eq, wavevector_zero_iff, and_true, true_and, or_false, false_or,
      and_false, false_and, true_or, or_true, not_false_eq_true, not_true_eq_false, gt_iff_lt]))
  if !(← getGoals).isEmpty then evalTactic (← `(tactic| first | tauto | aesop))

end
end S10Audit.CAS
