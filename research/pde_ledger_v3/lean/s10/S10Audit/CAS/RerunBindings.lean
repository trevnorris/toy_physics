import S10Audit.CAS.PYRerun
import S10Audit.CAS.WLRerun

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYRerun

def parallel0Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n35, n35, n35], ![n35, n39, n35], ![n35, n35, n39]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem parallel0Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root1_n1_matrix_cell0 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell1 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell2 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell3 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell4 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell5 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell6 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell7 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n1_matrix_cell8 rho mu sigma 0 PYLoci.parallelPoint

def parallel0Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n35, n35, n35], ![n35, n39, n35], ![n35, n35, n39], ![n40, n24, n24]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem parallel0Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) 2, PYLoci.parallelPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root1_n3_stacked_matrix_cell0 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell1 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell2 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell3 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell4 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell5 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell6 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell7 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell8 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell9 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell10 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n3_stacked_matrix_cell11 rho mu sigma 0 PYLoci.parallelPoint

def parallel0Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n40, n24, n24]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem parallel0Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Basis rho mu sigma = PYReference.parallel0 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root1_n6_nullspace_basis_cell0 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n6_nullspace_basis_cell1 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root1_n6_nullspace_basis_cell2 rho mu sigma 0 PYLoci.parallelPoint

theorem parallel0_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (parallel0Basis rho mu sigma) := by
  rw [parallel0Basis_reference rho mu sigma hr hs]
  exact PYReference.parallel0_independent

theorem parallel0_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (parallel0Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel0Basis rho mu sigma)) := by
  rw [parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs),
    parallel0Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.parallel0_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma PYLoci.parallelPoint a 0 hr hm]

theorem parallel0_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 0) (c • PYLoci.parallelPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel0Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c PYLoci.parallelPoint a hc,
    ← parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact parallel0_basis_complete rho mu sigma a hr hm hs hs1

def parallel1Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n42, n35, n35], ![n35, n35, n35], ![n35, n35, n35]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem parallel1Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root2_n1_matrix_cell0 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell1 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell2 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell3 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell4 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell5 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell6 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell7 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n1_matrix_cell8 rho mu sigma 0 PYLoci.parallelPoint _hr

def parallel1Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n42, n35, n35], ![n35, n35, n35], ![n35, n35, n35], ![n40, n24, n24]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem parallel1Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) 2, PYLoci.parallelPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root2_n3_stacked_matrix_cell0 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell1 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell2 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell3 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell4 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell5 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell6 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell7 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell8 rho mu sigma 0 PYLoci.parallelPoint _hr
  · exact q8_stratum1_root2_n3_stacked_matrix_cell9 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n3_stacked_matrix_cell10 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n3_stacked_matrix_cell11 rho mu sigma 0 PYLoci.parallelPoint

def parallel1Basis (rho mu sigma : ℝ) : Fin 2 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.parallelPoint)
    ((![![n24, n40, n24], ![n24, n24, n40]] : Fin 2 → Fin 3 → Expr Symbol) i j)

theorem parallel1Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Basis rho mu sigma = PYReference.parallel1 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum1_root2_n6_nullspace_basis_cell0 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n6_nullspace_basis_cell1 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n6_nullspace_basis_cell2 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n6_nullspace_basis_cell3 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n6_nullspace_basis_cell4 rho mu sigma 0 PYLoci.parallelPoint
  · exact q8_stratum1_root2_n6_nullspace_basis_cell5 rho mu sigma 0 PYLoci.parallelPoint

theorem parallel1_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (parallel1Basis rho mu sigma) := by
  rw [parallel1Basis_reference rho mu sigma hr hs]
  exact PYReference.parallel1_independent

theorem parallel1_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (parallel1Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel1Basis rho mu sigma)) := by
  rw [parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs),
    parallel1Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.parallel1_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma PYLoci.parallelPoint a 1 hr hm]

theorem parallel1_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 1) (c • PYLoci.parallelPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel1Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c PYLoci.parallelPoint a hc,
    ← parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact parallel1_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular0Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n39, n35, n35], ![n35, n35, n35], ![n35, n35, n39]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root1_n1_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n1_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint

def perpendicular0Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n39, n35, n35], ![n35, n35, n35], ![n35, n35, n39], ![n24, n40, n24]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) 2, PYLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root1_n3_stacked_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell9 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell10 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n3_stacked_matrix_cell11 rho mu sigma 0 PYLoci.perpendicularPoint

def perpendicular0Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n24, n40, n24]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Basis rho mu sigma = PYReference.perpendicular0 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root1_n6_nullspace_basis_cell0 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n6_nullspace_basis_cell1 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root1_n6_nullspace_basis_cell2 rho mu sigma 0 PYLoci.perpendicularPoint

theorem perpendicular0_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular0Basis rho mu sigma) := by
  rw [perpendicular0Basis_reference rho mu sigma hr hs]
  exact PYReference.perpendicular0_independent

theorem perpendicular0_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular0Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular0Basis rho mu sigma)) := by
  rw [perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular0Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular0_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma PYLoci.perpendicularPoint a 0 hr hm]

theorem perpendicular0_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 0) (c • PYLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular0Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c PYLoci.perpendicularPoint a hc,
    ← perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular0_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular1Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n45, n35, n35], ![n35, n46, n35], ![n35, n35, n35]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root2_n1_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n1_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint _hr

def perpendicular1Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n45, n35, n35], ![n35, n46, n35], ![n35, n35, n35], ![n24, n40, n24]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) 2, PYLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root2_n3_stacked_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint _hr
  · exact q8_stratum2_root2_n3_stacked_matrix_cell9 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root2_n3_stacked_matrix_cell10 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root2_n3_stacked_matrix_cell11 rho mu sigma 0 PYLoci.perpendicularPoint

def perpendicular1Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n24, n24, n40]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Basis rho mu sigma = PYReference.perpendicular1 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root2_n6_nullspace_basis_cell0 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root2_n6_nullspace_basis_cell1 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root2_n6_nullspace_basis_cell2 rho mu sigma 0 PYLoci.perpendicularPoint

theorem perpendicular1_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular1Basis rho mu sigma) := by
  rw [perpendicular1Basis_reference rho mu sigma hr hs]
  exact PYReference.perpendicular1_independent

theorem perpendicular1_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular1Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular1Basis rho mu sigma)) := by
  rw [perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular1Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular1_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma PYLoci.perpendicularPoint a 1 hr hm]

theorem perpendicular1_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 1) (c • PYLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular1Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c PYLoci.perpendicularPoint a hc,
    ← perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular1_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular2Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n35, n35, n35], ![n35, n48, n35], ![n35, n35, n50]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root3_n1_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint _hs _hr
  · exact q8_stratum2_root3_n1_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n1_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint _hs _hr

def perpendicular2Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n35, n35, n35], ![n35, n48, n35], ![n35, n35, n50], ![n24, n40, n24]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) 2, PYLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root3_n3_stacked_matrix_cell0 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell1 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell2 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell3 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell4 rho mu sigma 0 PYLoci.perpendicularPoint _hs _hr
  · exact q8_stratum2_root3_n3_stacked_matrix_cell5 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell6 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell7 rho mu sigma 0 PYLoci.perpendicularPoint _hr _hs
  · exact q8_stratum2_root3_n3_stacked_matrix_cell8 rho mu sigma 0 PYLoci.perpendicularPoint _hs _hr
  · exact q8_stratum2_root3_n3_stacked_matrix_cell9 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root3_n3_stacked_matrix_cell10 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root3_n3_stacked_matrix_cell11 rho mu sigma 0 PYLoci.perpendicularPoint

def perpendicular2Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 PYLoci.perpendicularPoint)
    ((![![n40, n24, n24]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Basis rho mu sigma = PYReference.perpendicular2 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact q8_stratum2_root3_n6_nullspace_basis_cell0 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root3_n6_nullspace_basis_cell1 rho mu sigma 0 PYLoci.perpendicularPoint
  · exact q8_stratum2_root3_n6_nullspace_basis_cell2 rho mu sigma 0 PYLoci.perpendicularPoint

theorem perpendicular2_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular2Basis rho mu sigma) := by
  rw [perpendicular2Basis_reference rho mu sigma hr hs]
  exact PYReference.perpendicular2_independent

theorem perpendicular2_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular2Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular2Basis rho mu sigma)) := by
  rw [perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular2Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular2_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma PYLoci.perpendicularPoint a 2 hr hm]

theorem perpendicular2_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 2) (c • PYLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular2Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c PYLoci.perpendicularPoint a hc,
    ← perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular2_basis_complete rho mu sigma a hr hm hs hs1

end PYRerun

namespace WLRerun

def parallel0Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n33, n33, n33], ![n33, n38, n33], ![n33, n33, n38]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem parallel0Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root1_n1_matrix_cell0 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell1 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell2 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell3 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell4 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell5 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell6 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell7 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n1_matrix_cell8 rho mu sigma 0 WLLoci.parallelPoint

def parallel0Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n33, n33, n33], ![n33, n38, n33], ![n33, n33, n38], ![n39, n13, n13]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem parallel0Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) 2, WLLoci.parallelPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root1_n3_stacked_matrix_cell0 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell1 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell2 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell3 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell4 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell5 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell6 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell7 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell8 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell9 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell10 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n3_stacked_matrix_cell11 rho mu sigma 0 WLLoci.parallelPoint

def parallel0Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n40, n13, n13]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem parallel0Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel0Basis rho mu sigma = WLReference.parallel0 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root1_n6_nullspace_basis_cell0 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n6_nullspace_basis_cell1 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root1_n6_nullspace_basis_cell2 rho mu sigma 0 WLLoci.parallelPoint

theorem parallel0_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (parallel0Basis rho mu sigma) := by
  rw [parallel0Basis_reference rho mu sigma hr hs]
  exact WLReference.parallel0_independent

theorem parallel0_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (parallel0Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel0Basis rho mu sigma)) := by
  rw [parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs),
    parallel0Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.parallel0_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma WLLoci.parallelPoint a 0 hr hm]

theorem parallel0_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 0) (c • WLLoci.parallelPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel0Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c WLLoci.parallelPoint a hc,
    ← parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact parallel0_basis_complete rho mu sigma a hr hm hs hs1

def parallel1Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n42, n33, n33], ![n33, n33, n33], ![n33, n33, n33]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem parallel1Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root2_n1_matrix_cell0 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell1 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell2 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell3 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell4 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell5 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell6 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell7 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n1_matrix_cell8 rho mu sigma 0 WLLoci.parallelPoint _hr

def parallel1Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n42, n33, n33], ![n33, n33, n33], ![n33, n33, n33], ![n39, n13, n13]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem parallel1Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) 2, WLLoci.parallelPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root2_n3_stacked_matrix_cell0 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell1 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell2 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell3 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell4 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell5 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell6 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell7 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell8 rho mu sigma 0 WLLoci.parallelPoint _hr
  · exact stratum1_root2_n3_stacked_matrix_cell9 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n3_stacked_matrix_cell10 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n3_stacked_matrix_cell11 rho mu sigma 0 WLLoci.parallelPoint

def parallel1Basis (rho mu sigma : ℝ) : Fin 2 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.parallelPoint)
    ((![![n13, n13, n40], ![n13, n40, n13]] : Fin 2 → Fin 3 → Expr Symbol) i j)

theorem parallel1Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    parallel1Basis rho mu sigma = WLReference.parallel1 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum1_root2_n6_nullspace_basis_cell0 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n6_nullspace_basis_cell1 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n6_nullspace_basis_cell2 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n6_nullspace_basis_cell3 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n6_nullspace_basis_cell4 rho mu sigma 0 WLLoci.parallelPoint
  · exact stratum1_root2_n6_nullspace_basis_cell5 rho mu sigma 0 WLLoci.parallelPoint

theorem parallel1_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (parallel1Basis rho mu sigma) := by
  rw [parallel1Basis_reference rho mu sigma hr hs]
  exact WLReference.parallel1_independent

theorem parallel1_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (parallel1Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel1Basis rho mu sigma)) := by
  rw [parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs),
    parallel1Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.parallel1_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma WLLoci.parallelPoint a 1 hr hm]

theorem parallel1_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 1) (c • WLLoci.parallelPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (parallel1Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c WLLoci.parallelPoint a hc,
    ← parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact parallel1_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular0Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n49, n33, n33], ![n33, n52, n55], ![n33, n55, n38]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root1_n1_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n1_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint

def perpendicular0Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n49, n33, n33], ![n33, n52, n55], ![n33, n55, n38], ![n13, n39, n56]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) 2, WLLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root1_n3_stacked_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell9 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell10 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n3_stacked_matrix_cell11 rho mu sigma 0 WLLoci.perpendicularPoint

def perpendicular0Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n13, n58, n40]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular0Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular0Basis rho mu sigma = WLReference.perpendicular0 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root1_n6_nullspace_basis_cell0 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n6_nullspace_basis_cell1 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root1_n6_nullspace_basis_cell2 rho mu sigma 0 WLLoci.perpendicularPoint

theorem perpendicular0_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular0Basis rho mu sigma) := by
  rw [perpendicular0Basis_reference rho mu sigma hr hs]
  exact WLReference.perpendicular0_independent

theorem perpendicular0_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular0Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular0Basis rho mu sigma)) := by
  rw [perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular0Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular0_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma WLLoci.perpendicularPoint a 0 hr hm]

theorem perpendicular0_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 0) (c • WLLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular0Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c WLLoci.perpendicularPoint a hc,
    ← perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular0_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular1Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n62, n33, n33], ![n33, n63, n55], ![n33, n55, n64]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root2_n1_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n1_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint _hr

def perpendicular1Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n62, n33, n33], ![n33, n63, n55], ![n33, n55, n64], ![n13, n39, n56]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) 2, WLLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root2_n3_stacked_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint _hr
  · exact stratum2_root2_n3_stacked_matrix_cell9 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root2_n3_stacked_matrix_cell10 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root2_n3_stacked_matrix_cell11 rho mu sigma 0 WLLoci.perpendicularPoint

def perpendicular1Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n13, n70, n40]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular1Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular1Basis rho mu sigma = WLReference.perpendicular1 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root2_n6_nullspace_basis_cell0 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root2_n6_nullspace_basis_cell1 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root2_n6_nullspace_basis_cell2 rho mu sigma 0 WLLoci.perpendicularPoint

theorem perpendicular1_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular1Basis rho mu sigma) := by
  rw [perpendicular1Basis_reference rho mu sigma hr hs]
  exact WLReference.perpendicular1_independent

theorem perpendicular1_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular1Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular1Basis rho mu sigma)) := by
  rw [perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular1Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular1_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma WLLoci.perpendicularPoint a 1 hr hm]

theorem perpendicular1_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 1) (c • WLLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular1Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c WLLoci.perpendicularPoint a hc,
    ← perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular1_basis_complete rho mu sigma a hr hm hs hs1

def perpendicular2Matrix (rho mu sigma : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n33, n33, n33], ![n33, n77, n55], ![n33, n55, n83]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Matrix_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Matrix rho mu sigma = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root3_n1_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint _hs _hr
  · exact stratum2_root3_n1_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n1_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint _hs _hr

def perpendicular2Stack (rho mu sigma : ℝ) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n33, n33, n33], ![n33, n77, n55], ![n33, n55, n83], ![n13, n39, n56]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Stack_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Stack rho mu sigma = ![((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) 0, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) 1, ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) 2, WLLoci.perpendicularPoint] := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root3_n3_stacked_matrix_cell0 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell1 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell2 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell3 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell4 rho mu sigma 0 WLLoci.perpendicularPoint _hs _hr
  · exact stratum2_root3_n3_stacked_matrix_cell5 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell6 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell7 rho mu sigma 0 WLLoci.perpendicularPoint _hr _hs
  · exact stratum2_root3_n3_stacked_matrix_cell8 rho mu sigma 0 WLLoci.perpendicularPoint _hs _hr
  · exact stratum2_root3_n3_stacked_matrix_cell9 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root3_n3_stacked_matrix_cell10 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root3_n3_stacked_matrix_cell11 rho mu sigma 0 WLLoci.perpendicularPoint

def perpendicular2Basis (rho mu sigma : ℝ) : Fin 1 → Vec 3 :=
  fun i j => Expr.eval (values rho mu sigma 0 WLLoci.perpendicularPoint)
    ((![![n40, n13, n13]] : Fin 1 → Fin 3 → Expr Symbol) i j)

theorem perpendicular2Basis_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    perpendicular2Basis rho mu sigma = WLReference.perpendicular2 := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact stratum2_root3_n6_nullspace_basis_cell0 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root3_n6_nullspace_basis_cell1 rho mu sigma 0 WLLoci.perpendicularPoint
  · exact stratum2_root3_n6_nullspace_basis_cell2 rho mu sigma 0 WLLoci.perpendicularPoint

theorem perpendicular2_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    LinearIndependent ℝ (perpendicular2Basis rho mu sigma) := by
  rw [perpendicular2Basis_reference rho mu sigma hr hs]
  exact WLReference.perpendicular2_independent

theorem perpendicular2_basis_complete (rho mu sigma : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    (perpendicular2Matrix rho mu sigma).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular2Basis rho mu sigma)) := by
  rw [perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs),
    perpendicular2Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular2_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma WLLoci.perpendicularPoint a 2 hr hm]

theorem perpendicular2_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma
      (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 2) (c • WLLoci.perpendicularPoint)).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ (Set.range (perpendicular2Basis rho mu sigma)) := by
  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c WLLoci.perpendicularPoint a hc,
    ← perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)]
  exact perpendicular2_basis_complete rho mu sigma a hr hm hs hs1

end WLRerun

end
end S10Audit.CAS
