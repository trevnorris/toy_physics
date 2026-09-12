import S10Audit.CAS.RerunSupport

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section

namespace PYReference

def parallel0 : Fin 1 → Vec 3 := ![![1, 0, 0]]

theorem parallel0_independent : LinearIndependent ℝ parallel0 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 0
    simpa [parallel0, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem parallel0_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    parallel0 j ∈ rootMode sigma PYLoci.parallelPoint 0 := by
  change normalizedOperator 0 sigma (0) PYLoci.parallelPoint (parallel0 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [parallel0, PYLoci.parallelPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem parallel0_dimensions (sigma : ℝ) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.parallelPoint 0) = 1 ∧
      Module.finrank ℝ (rootMode sigma PYLoci.parallelPoint 0 ⊓ transverseSpace PYLoci.parallelPoint : Submodule ℝ (Vec 3)) = 0 := by
  exact zero_counts (sigma := sigma) PYLoci.parallelPoint_nonzero

theorem parallel0_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range parallel0) = rootMode sigma PYLoci.parallelPoint 0 := by
  exact complete_of_dimension parallel0 _ parallel0_independent
    (parallel0_members sigma (ne_of_gt hs)) (parallel0_dimensions sigma hs hs1).1

def parallel1 : Fin 2 → Vec 3 := ![![0, 1, 0], ![0, 0, 1]]

theorem parallel1_independent : LinearIndependent ℝ parallel1 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 1
    simpa [parallel1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h
  · have h := congrFun hc 2
    simpa [parallel1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem parallel1_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 2) :
    parallel1 j ∈ rootMode sigma PYLoci.parallelPoint 1 := by
  change normalizedOperator 0 sigma (normSq PYLoci.parallelPoint) PYLoci.parallelPoint (parallel1 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [parallel1, PYLoci.parallelPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem parallel1_dimensions (sigma : ℝ) (hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.parallelPoint 1) = 2 ∧
      Module.finrank ℝ (rootMode sigma PYLoci.parallelPoint 1 ⊓ transverseSpace PYLoci.parallelPoint : Submodule ℝ (Vec 3)) = 2 := by
  apply parallel_counts (ne_of_gt hs) PYLoci.parallelPoint_nonzero
  norm_num [perpSq, normSq, dot, Fin.sum_univ_three, PYLoci.parallelPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

theorem parallel1_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range parallel1) = rootMode sigma PYLoci.parallelPoint 1 := by
  exact complete_of_dimension parallel1 _ parallel1_independent
    (parallel1_members sigma (ne_of_gt hs)) (parallel1_dimensions sigma hs hs1).1

def perpendicular0 : Fin 1 → Vec 3 := ![![0, 1, 0]]

theorem perpendicular0_independent : LinearIndependent ℝ perpendicular0 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 1
    simpa [perpendicular0, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular0_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular0 j ∈ rootMode sigma PYLoci.perpendicularPoint 0 := by
  change normalizedOperator 0 sigma (0) PYLoci.perpendicularPoint (perpendicular0 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular0, PYLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem perpendicular0_dimensions (sigma : ℝ) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 0) = 1 ∧
      Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 0 ⊓ transverseSpace PYLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 0 := by
  exact zero_counts (sigma := sigma) PYLoci.perpendicularPoint_nonzero

theorem perpendicular0_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular0) = rootMode sigma PYLoci.perpendicularPoint 0 := by
  exact complete_of_dimension perpendicular0 _ perpendicular0_independent
    (perpendicular0_members sigma (ne_of_gt hs)) (perpendicular0_dimensions sigma hs hs1).1

def perpendicular1 : Fin 1 → Vec 3 := ![![0, 0, 1]]

theorem perpendicular1_independent : LinearIndependent ℝ perpendicular1 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 2
    simpa [perpendicular1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular1_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular1 j ∈ rootMode sigma PYLoci.perpendicularPoint 1 := by
  change normalizedOperator 0 sigma (normSq PYLoci.perpendicularPoint) PYLoci.perpendicularPoint (perpendicular1 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular1, PYLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem perpendicular1_dimensions (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 1) = 1 ∧
      Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 1 ⊓ transverseSpace PYLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 1 := by
  have h := perpendicular_counts hs hs1 PYLoci.perpendicularPoint_nonzero (by norm_num [PYLoci.perpendicularPoint] : PYLoci.perpendicularPoint 0 = 0)
  exact ⟨h.1, h.2.1⟩

theorem perpendicular1_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular1) = rootMode sigma PYLoci.perpendicularPoint 1 := by
  exact complete_of_dimension perpendicular1 _ perpendicular1_independent
    (perpendicular1_members sigma (ne_of_gt hs)) (perpendicular1_dimensions sigma hs hs1).1

def perpendicular2 : Fin 1 → Vec 3 := ![![1, 0, 0]]

theorem perpendicular2_independent : LinearIndependent ℝ perpendicular2 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 0
    simpa [perpendicular2, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular2_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular2 j ∈ rootMode sigma PYLoci.perpendicularPoint 2 := by
  change normalizedOperator 0 sigma (extraValue 0 sigma PYLoci.perpendicularPoint) PYLoci.perpendicularPoint (perpendicular2 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular2, PYLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]
  all_goals cas_equal

theorem perpendicular2_dimensions (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 2) = 1 ∧
      Module.finrank ℝ (rootMode sigma PYLoci.perpendicularPoint 2 ⊓ transverseSpace PYLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 1 := by
  have h := perpendicular_counts hs hs1 PYLoci.perpendicularPoint_nonzero (by norm_num [PYLoci.perpendicularPoint] : PYLoci.perpendicularPoint 0 = 0)
  exact ⟨h.2.2.1, h.2.2.2⟩

theorem perpendicular2_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular2) = rootMode sigma PYLoci.perpendicularPoint 2 := by
  exact complete_of_dimension perpendicular2 _ perpendicular2_independent
    (perpendicular2_members sigma (ne_of_gt hs)) (perpendicular2_dimensions sigma hs hs1).1

end PYReference

namespace WLReference

def parallel0 : Fin 1 → Vec 3 := ![![1, 0, 0]]

theorem parallel0_independent : LinearIndependent ℝ parallel0 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 0
    simpa [parallel0, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem parallel0_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    parallel0 j ∈ rootMode sigma WLLoci.parallelPoint 0 := by
  change normalizedOperator 0 sigma (0) WLLoci.parallelPoint (parallel0 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [parallel0, WLLoci.parallelPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem parallel0_dimensions (sigma : ℝ) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma WLLoci.parallelPoint 0) = 1 ∧
      Module.finrank ℝ (rootMode sigma WLLoci.parallelPoint 0 ⊓ transverseSpace WLLoci.parallelPoint : Submodule ℝ (Vec 3)) = 0 := by
  exact zero_counts (sigma := sigma) WLLoci.parallelPoint_nonzero

theorem parallel0_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range parallel0) = rootMode sigma WLLoci.parallelPoint 0 := by
  exact complete_of_dimension parallel0 _ parallel0_independent
    (parallel0_members sigma (ne_of_gt hs)) (parallel0_dimensions sigma hs hs1).1

def parallel1 : Fin 2 → Vec 3 := ![![0, 0, 1], ![0, 1, 0]]

theorem parallel1_independent : LinearIndependent ℝ parallel1 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 2
    simpa [parallel1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h
  · have h := congrFun hc 1
    simpa [parallel1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem parallel1_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 2) :
    parallel1 j ∈ rootMode sigma WLLoci.parallelPoint 1 := by
  change normalizedOperator 0 sigma (normSq WLLoci.parallelPoint) WLLoci.parallelPoint (parallel1 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [parallel1, WLLoci.parallelPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem parallel1_dimensions (sigma : ℝ) (hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma WLLoci.parallelPoint 1) = 2 ∧
      Module.finrank ℝ (rootMode sigma WLLoci.parallelPoint 1 ⊓ transverseSpace WLLoci.parallelPoint : Submodule ℝ (Vec 3)) = 2 := by
  apply parallel_counts (ne_of_gt hs) WLLoci.parallelPoint_nonzero
  norm_num [perpSq, normSq, dot, Fin.sum_univ_three, WLLoci.parallelPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

theorem parallel1_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range parallel1) = rootMode sigma WLLoci.parallelPoint 1 := by
  exact complete_of_dimension parallel1 _ parallel1_independent
    (parallel1_members sigma (ne_of_gt hs)) (parallel1_dimensions sigma hs hs1).1

def perpendicular0 : Fin 1 → Vec 3 := ![![0, -54, 1]]

theorem perpendicular0_independent : LinearIndependent ℝ perpendicular0 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 2
    simpa [perpendicular0, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular0_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular0 j ∈ rootMode sigma WLLoci.perpendicularPoint 0 := by
  change normalizedOperator 0 sigma (0) WLLoci.perpendicularPoint (perpendicular0 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular0, WLLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem perpendicular0_dimensions (sigma : ℝ) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 0) = 1 ∧
      Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 0 ⊓ transverseSpace WLLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 0 := by
  exact zero_counts (sigma := sigma) WLLoci.perpendicularPoint_nonzero

theorem perpendicular0_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular0) = rootMode sigma WLLoci.perpendicularPoint 0 := by
  exact complete_of_dimension perpendicular0 _ perpendicular0_independent
    (perpendicular0_members sigma (ne_of_gt hs)) (perpendicular0_dimensions sigma hs hs1).1

def perpendicular1 : Fin 1 → Vec 3 := ![![0, 1/54, 1]]

theorem perpendicular1_independent : LinearIndependent ℝ perpendicular1 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 2
    simpa [perpendicular1, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular1_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular1 j ∈ rootMode sigma WLLoci.perpendicularPoint 1 := by
  change normalizedOperator 0 sigma (normSq WLLoci.perpendicularPoint) WLLoci.perpendicularPoint (perpendicular1 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular1, WLLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]

theorem perpendicular1_dimensions (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 1) = 1 ∧
      Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 1 ⊓ transverseSpace WLLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 1 := by
  have h := perpendicular_counts hs hs1 WLLoci.perpendicularPoint_nonzero (by norm_num [WLLoci.perpendicularPoint] : WLLoci.perpendicularPoint 0 = 0)
  exact ⟨h.1, h.2.1⟩

theorem perpendicular1_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular1) = rootMode sigma WLLoci.perpendicularPoint 1 := by
  exact complete_of_dimension perpendicular1 _ perpendicular1_independent
    (perpendicular1_members sigma (ne_of_gt hs)) (perpendicular1_dimensions sigma hs hs1).1

def perpendicular2 : Fin 1 → Vec 3 := ![![1, 0, 0]]

theorem perpendicular2_independent : LinearIndependent ℝ perpendicular2 := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  fin_cases j
  · have h := congrFun hc 0
    simpa [perpendicular2, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h

theorem perpendicular2_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin 1) :
    perpendicular2 j ∈ rootMode sigma WLLoci.perpendicularPoint 2 := by
  change normalizedOperator 0 sigma (extraValue 0 sigma WLLoci.perpendicularPoint) WLLoci.perpendicularPoint (perpendicular2 j) = 0
  ext i
  fin_cases j
  all_goals fin_cases i
  all_goals
    norm_num [perpendicular2, WLLoci.perpendicularPoint, normalizedOperator, extraValue, extraNumerator, perpSq,
      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]
  all_goals cas_equal

theorem perpendicular2_dimensions (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 2) = 1 ∧
      Module.finrank ℝ (rootMode sigma WLLoci.perpendicularPoint 2 ⊓ transverseSpace WLLoci.perpendicularPoint : Submodule ℝ (Vec 3)) = 1 := by
  have h := perpendicular_counts hs hs1 WLLoci.perpendicularPoint_nonzero (by norm_num [WLLoci.perpendicularPoint] : WLLoci.perpendicularPoint 0 = 0)
  exact ⟨h.2.2.1, h.2.2.2⟩

theorem perpendicular2_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Submodule.span ℝ (Set.range perpendicular2) = rootMode sigma WLLoci.perpendicularPoint 2 := by
  exact complete_of_dimension perpendicular2 _ perpendicular2_independent
    (perpendicular2_members sigma (ne_of_gt hs)) (perpendicular2_dimensions sigma hs hs1).1

end WLReference

open Lean Elab Tactic in
elab "rerun_equal" : tactic => do
  evalTactic (← `(tactic| norm_num [PYReference.parallel0, PYReference.parallel1, PYReference.perpendicular0, PYReference.perpendicular1, PYReference.perpendicular2, WLReference.parallel0, WLReference.parallel1, WLReference.perpendicular0, WLReference.perpendicular1, WLReference.perpendicular2, PYLoci.parallelPoint, PYLoci.perpendicularPoint, WLLoci.parallelPoint, WLLoci.perpendicularPoint,
    referenceMatrix, referenceRoot, Matrix.det_fin_three, Matrix.mulVec, dotProduct,
    extraConeValue, extraValue, extraNumerator, perpSq, coneValue, normSq, dot, Fin.sum_univ_three,
    Fin.ext_iff, Fin.coe_ofNat_eq_mod, -Fin.val_eq_zero_iff,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| field_simp))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| ring_nf))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| norm_num))

open Lean Elab Tactic in
elab "rerun_eval" h:Lean.Parser.Tactic.rwRule : tactic => do
  evalTactic (← `(tactic| rw [$h]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| rerun_equal))

end
end S10Audit.CAS
