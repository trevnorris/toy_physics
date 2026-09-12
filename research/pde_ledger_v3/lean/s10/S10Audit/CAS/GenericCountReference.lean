import S10Audit.CAS.GenericCountSupport

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYGenericCountReference

def root0Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n180, PY.n40, PY.n43], ![PY.n40, PY.n183, PY.n52], ![PY.n43, PY.n52, PY.n185]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root0Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root0Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 0) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root1_n1_matrix_cell0 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell1 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell2 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell3 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell4 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell5 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell6 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell7 rho mu sigma z k
  · exact PY.root1_n1_matrix_cell8 rho mu sigma z k

def root0Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n180, PY.n40, PY.n43], ![PY.n40, PY.n183, PY.n52], ![PY.n43, PY.n52, PY.n185], ![PY.n15, PY.n0, PY.n4]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root0Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root0Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 0) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root1_n3_stacked_matrix_cell0 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell1 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell2 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell3 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell4 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell5 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell6 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell7 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell8 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact PY.root1_n3_stacked_matrix_cell11 rho mu sigma z k

def root0Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => PY.basis rho mu sigma z k 0

theorem root0_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root0Basis rho mu sigma z k) := by
  exact PY.basis_independent rho mu sigma z k 0 h

theorem root0_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root0Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root0Basis rho mu sigma z k)) := by
  rw [root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 0 hr hm h]
  have hb : Set.range (root0Basis rho mu sigma z k) = {referenceBasis sigma k 0} := by
    change Set.range (fun _ : Fin 1 => PY.basis rho mu sigma z k 0) = _
    rw [Set.range_const, PY.basis_reference rho mu sigma z k 0 h]
  rw [hb]

theorem root0_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root0Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 0 := by
  rw [root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 0 hr hm

theorem root0_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root0Stack rho mu sigma z k).mulVecLin = rootMode sigma k 0 ⊓ transverseSpace k := by
  rw [root0Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root0_kernel rho mu sigma z k hr hm h]

theorem root0_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root0Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root0_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 0 h).1

theorem root0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root0Stack rho mu sigma z k) = 0 := by
  rw [matrixNullity, root0_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 0 h).2

theorem root0_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root0Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root0Matrix rho mu sigma z k) (root0_nullity rho mu sigma z k hr hm h)

theorem root0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root0Stack rho mu sigma z k).rank = 3 := by
  exact matrix_rank_of_nullity (root0Stack rho mu sigma z k) (root0_transverse_nullity rho mu sigma z k hr hm h)

theorem root0_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (PYGenericCountReference.root0Matrix rho mu sigma z k) (PYGenericCountReference.root0Stack rho mu sigma z k) = (1 : ℤ) := by
  rw [nullityDifference, root0_nullity rho mu sigma z k hr hm h, root0_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root0_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (PYGenericCountReference.root0Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (PYGenericCountReference.root0Basis rho mu sigma z k) (PYGenericCountReference.root0Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root0_basis_count, root0_nullity rho mu sigma z k hr hm h]
  norm_num

def root1Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n204, PY.n40, PY.n43], ![PY.n40, PY.n53, PY.n52], ![PY.n43, PY.n52, PY.n34]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root1Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root1Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 1) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root2_n1_matrix_cell0 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell1 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell2 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell3 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell4 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell5 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell6 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell7 rho mu sigma z k _hr
  · exact PY.root2_n1_matrix_cell8 rho mu sigma z k _hr

def root1Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n204, PY.n40, PY.n43], ![PY.n40, PY.n53, PY.n52], ![PY.n43, PY.n52, PY.n34], ![PY.n15, PY.n0, PY.n4]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root1Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root1Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 1) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root2_n3_stacked_matrix_cell0 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell1 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell2 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell3 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell4 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell5 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell6 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell7 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell8 rho mu sigma z k _hr
  · exact PY.root2_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact PY.root2_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact PY.root2_n3_stacked_matrix_cell11 rho mu sigma z k

def root1Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => PY.basis rho mu sigma z k 1

theorem root1_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root1Basis rho mu sigma z k) := by
  exact PY.basis_independent rho mu sigma z k 1 h

theorem root1_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root1Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root1Basis rho mu sigma z k)) := by
  rw [root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 1 hr hm h]
  have hb : Set.range (root1Basis rho mu sigma z k) = {referenceBasis sigma k 1} := by
    change Set.range (fun _ : Fin 1 => PY.basis rho mu sigma z k 1) = _
    rw [Set.range_const, PY.basis_reference rho mu sigma z k 1 h]
  rw [hb]

theorem root1_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root1Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 1 := by
  rw [root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 1 hr hm

theorem root1_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root1Stack rho mu sigma z k).mulVecLin = rootMode sigma k 1 ⊓ transverseSpace k := by
  rw [root1Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root1_kernel rho mu sigma z k hr hm h]

theorem root1_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root1Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root1_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 1 h).1

theorem root1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root1Stack rho mu sigma z k) = 1 := by
  rw [matrixNullity, root1_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 1 h).2

theorem root1_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root1Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root1Matrix rho mu sigma z k) (root1_nullity rho mu sigma z k hr hm h)

theorem root1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root1Stack rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root1Stack rho mu sigma z k) (root1_transverse_nullity rho mu sigma z k hr hm h)

theorem root1_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (PYGenericCountReference.root1Matrix rho mu sigma z k) (PYGenericCountReference.root1Stack rho mu sigma z k) = (0 : ℤ) := by
  rw [nullityDifference, root1_nullity rho mu sigma z k hr hm h, root1_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root1_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (PYGenericCountReference.root1Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (PYGenericCountReference.root1Basis rho mu sigma z k) (PYGenericCountReference.root1Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root1_basis_count, root1_nullity rho mu sigma z k hr hm h]
  norm_num

def root2Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n222, PY.n40, PY.n43], ![PY.n40, PY.n227, PY.n52], ![PY.n43, PY.n52, PY.n231]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root2Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root2Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 2) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root3_n1_matrix_cell0 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell1 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell2 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell3 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell4 rho mu sigma z k _hs _hr
  · exact PY.root3_n1_matrix_cell5 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell6 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell7 rho mu sigma z k _hr _hs
  · exact PY.root3_n1_matrix_cell8 rho mu sigma z k _hs _hr

def root2Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![PY.n222, PY.n40, PY.n43], ![PY.n40, PY.n227, PY.n52], ![PY.n43, PY.n52, PY.n231], ![PY.n15, PY.n0, PY.n4]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root2Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root2Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 2) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact PY.root3_n3_stacked_matrix_cell0 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell1 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell2 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell3 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell4 rho mu sigma z k _hs _hr
  · exact PY.root3_n3_stacked_matrix_cell5 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell6 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell7 rho mu sigma z k _hr _hs
  · exact PY.root3_n3_stacked_matrix_cell8 rho mu sigma z k _hs _hr
  · exact PY.root3_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact PY.root3_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact PY.root3_n3_stacked_matrix_cell11 rho mu sigma z k

def root2Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => PY.basis rho mu sigma z k 2

theorem root2_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root2Basis rho mu sigma z k) := by
  exact PY.basis_independent rho mu sigma z k 2 h

theorem root2_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root2Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root2Basis rho mu sigma z k)) := by
  rw [root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 2 hr hm h]
  have hb : Set.range (root2Basis rho mu sigma z k) = {referenceBasis sigma k 2} := by
    change Set.range (fun _ : Fin 1 => PY.basis rho mu sigma z k 2) = _
    rw [Set.range_const, PY.basis_reference rho mu sigma z k 2 h]
  rw [hb]

theorem root2_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root2Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 2 := by
  rw [root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 2 hr hm

theorem root2_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root2Stack rho mu sigma z k).mulVecLin = rootMode sigma k 2 ⊓ transverseSpace k := by
  rw [root2Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root2_kernel rho mu sigma z k hr hm h]

theorem root2_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root2Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root2_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 2 h).1

theorem root2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root2Stack rho mu sigma z k) = 0 := by
  rw [matrixNullity, root2_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 2 h).2

theorem root2_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root2Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root2Matrix rho mu sigma z k) (root2_nullity rho mu sigma z k hr hm h)

theorem root2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root2Stack rho mu sigma z k).rank = 3 := by
  exact matrix_rank_of_nullity (root2Stack rho mu sigma z k) (root2_transverse_nullity rho mu sigma z k hr hm h)

theorem root2_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (PYGenericCountReference.root2Matrix rho mu sigma z k) (PYGenericCountReference.root2Stack rho mu sigma z k) = (1 : ℤ) := by
  rw [nullityDifference, root2_nullity rho mu sigma z k hr hm h, root2_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root2_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (PYGenericCountReference.root2Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (PYGenericCountReference.root2Basis rho mu sigma z k) (PYGenericCountReference.root2Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root2_basis_count, root2_nullity rho mu sigma z k hr hm h]
  norm_num

end PYGenericCountReference

namespace WLGenericCountReference

def root0Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n93, WL.n35, WL.n36], ![WL.n35, WL.n94, WL.n40], ![WL.n36, WL.n40, WL.n95]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root0Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root0Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 0) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root1_n1_matrix_cell0 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell1 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell2 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell3 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell4 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell5 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell6 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell7 rho mu sigma z k
  · exact WL.root1_n1_matrix_cell8 rho mu sigma z k

def root0Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n93, WL.n35, WL.n36], ![WL.n35, WL.n94, WL.n40], ![WL.n36, WL.n40, WL.n95], ![WL.n14, WL.n0, WL.n2]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root0Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root0Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 0) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root1_n3_stacked_matrix_cell0 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell1 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell2 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell3 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell4 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell5 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell6 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell7 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell8 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact WL.root1_n3_stacked_matrix_cell11 rho mu sigma z k

def root0Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => WL.basis rho mu sigma z k 0

theorem root0_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root0Basis rho mu sigma z k) := by
  exact WL.basis_independent rho mu sigma z k 0 h

theorem root0_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root0Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root0Basis rho mu sigma z k)) := by
  rw [root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 0 hr hm h]
  have hb : Set.range (root0Basis rho mu sigma z k) = {referenceBasis sigma k 0} := by
    change Set.range (fun _ : Fin 1 => WL.basis rho mu sigma z k 0) = _
    rw [Set.range_const, WL.basis_reference rho mu sigma z k 0 h]
  rw [hb]

theorem root0_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root0Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 0 := by
  rw [root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 0 hr hm

theorem root0_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root0Stack rho mu sigma z k).mulVecLin = rootMode sigma k 0 ⊓ transverseSpace k := by
  rw [root0Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root0Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root0_kernel rho mu sigma z k hr hm h]

theorem root0_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root0Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root0_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 0 h).1

theorem root0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root0Stack rho mu sigma z k) = 0 := by
  rw [matrixNullity, root0_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 0 h).2

theorem root0_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root0Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root0Matrix rho mu sigma z k) (root0_nullity rho mu sigma z k hr hm h)

theorem root0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root0Stack rho mu sigma z k).rank = 3 := by
  exact matrix_rank_of_nullity (root0Stack rho mu sigma z k) (root0_transverse_nullity rho mu sigma z k hr hm h)

theorem root0_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (WLGenericCountReference.root0Matrix rho mu sigma z k) (WLGenericCountReference.root0Stack rho mu sigma z k) = (1 : ℤ) := by
  rw [nullityDifference, root0_nullity rho mu sigma z k hr hm h, root0_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root0_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (WLGenericCountReference.root0Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (WLGenericCountReference.root0Basis rho mu sigma z k) (WLGenericCountReference.root0Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root0_basis_count, root0_nullity rho mu sigma z k hr hm h]
  norm_num

def root1Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n114, WL.n35, WL.n36], ![WL.n35, WL.n115, WL.n40], ![WL.n36, WL.n40, WL.n116]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root1Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root1Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 1) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root2_n1_matrix_cell0 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell1 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell2 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell3 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell4 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell5 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell6 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell7 rho mu sigma z k _hr
  · exact WL.root2_n1_matrix_cell8 rho mu sigma z k _hr

def root1Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n114, WL.n35, WL.n36], ![WL.n35, WL.n115, WL.n40], ![WL.n36, WL.n40, WL.n116], ![WL.n14, WL.n0, WL.n2]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root1Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root1Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 1) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root2_n3_stacked_matrix_cell0 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell1 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell2 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell3 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell4 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell5 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell6 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell7 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell8 rho mu sigma z k _hr
  · exact WL.root2_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact WL.root2_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact WL.root2_n3_stacked_matrix_cell11 rho mu sigma z k

def root1Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => WL.basis rho mu sigma z k 1

theorem root1_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root1Basis rho mu sigma z k) := by
  exact WL.basis_independent rho mu sigma z k 1 h

theorem root1_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root1Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root1Basis rho mu sigma z k)) := by
  rw [root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 1 hr hm h]
  have hb : Set.range (root1Basis rho mu sigma z k) = {referenceBasis sigma k 1} := by
    change Set.range (fun _ : Fin 1 => WL.basis rho mu sigma z k 1) = _
    rw [Set.range_const, WL.basis_reference rho mu sigma z k 1 h]
  rw [hb]

theorem root1_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root1Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 1 := by
  rw [root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 1 hr hm

theorem root1_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root1Stack rho mu sigma z k).mulVecLin = rootMode sigma k 1 ⊓ transverseSpace k := by
  rw [root1Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root1Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root1_kernel rho mu sigma z k hr hm h]

theorem root1_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root1Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root1_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 1 h).1

theorem root1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root1Stack rho mu sigma z k) = 1 := by
  rw [matrixNullity, root1_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 1 h).2

theorem root1_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root1Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root1Matrix rho mu sigma z k) (root1_nullity rho mu sigma z k hr hm h)

theorem root1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root1Stack rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root1Stack rho mu sigma z k) (root1_transverse_nullity rho mu sigma z k hr hm h)

theorem root1_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (WLGenericCountReference.root1Matrix rho mu sigma z k) (WLGenericCountReference.root1Stack rho mu sigma z k) = (0 : ℤ) := by
  rw [nullityDifference, root1_nullity rho mu sigma z k hr hm h, root1_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root1_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (WLGenericCountReference.root1Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (WLGenericCountReference.root1Basis rho mu sigma z k) (WLGenericCountReference.root1Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root1_basis_count, root1_nullity rho mu sigma z k hr hm h]
  norm_num

def root2Matrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n133, WL.n35, WL.n36], ![WL.n35, WL.n137, WL.n40], ![WL.n36, WL.n40, WL.n140]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem root2Matrix_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root2Matrix rho mu sigma z k = ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 2) k) := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root3_n1_matrix_cell0 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell1 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell2 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell3 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell4 rho mu sigma z k _hs _hr
  · exact WL.root3_n1_matrix_cell5 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell6 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell7 rho mu sigma z k _hr _hs
  · exact WL.root3_n1_matrix_cell8 rho mu sigma z k _hs _hr

def root2Stack (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 4) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![WL.n133, WL.n35, WL.n36], ![WL.n35, WL.n137, WL.n40], ![WL.n36, WL.n40, WL.n140], ![WL.n14, WL.n0, WL.n2]] : Fin 4 → Fin 3 → Expr Symbol) i j)

theorem root2Stack_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :
    root2Stack rho mu sigma z k = appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k 2) k) k := by
  ext i j
  fin_cases i
  all_goals fin_cases j
  · exact WL.root3_n3_stacked_matrix_cell0 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell1 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell2 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell3 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell4 rho mu sigma z k _hs _hr
  · exact WL.root3_n3_stacked_matrix_cell5 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell6 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell7 rho mu sigma z k _hr _hs
  · exact WL.root3_n3_stacked_matrix_cell8 rho mu sigma z k _hs _hr
  · exact WL.root3_n3_stacked_matrix_cell9 rho mu sigma z k
  · exact WL.root3_n3_stacked_matrix_cell10 rho mu sigma z k
  · exact WL.root3_n3_stacked_matrix_cell11 rho mu sigma z k

def root2Basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 1 → Vec 3 :=
  fun _ => WL.basis rho mu sigma z k 2

theorem root2_basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    LinearIndependent ℝ (root2Basis rho mu sigma z k) := by
  exact WL.basis_independent rho mu sigma z k 2 h

theorem root2_basis_complete (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) (a : Vec 3) :
    (root2Matrix rho mu sigma z k).mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range (root2Basis rho mu sigma z k)) := by
  rw [root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_basis_complete rho mu sigma k a 2 hr hm h]
  have hb : Set.range (root2Basis rho mu sigma z k) = {referenceBasis sigma k 2} := by
    change Set.range (fun _ : Fin 1 => WL.basis rho mu sigma z k 2) = _
    rw [Set.range_const, WL.basis_reference rho mu sigma z k 2 h]
  rw [hb]

theorem root2_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root2Matrix rho mu sigma z k).mulVecLin = rootMode sigma k 2 := by
  rw [root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1)]
  ext a
  exact coordinate_kernel_mode rho mu sigma k a 2 hr hm

theorem root2_stack_kernel (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    LinearMap.ker (root2Stack rho mu sigma z k).mulVecLin = rootMode sigma k 2 ⊓ transverseSpace k := by
  rw [root2Stack_reference rho mu sigma z k hr (ne_of_gt h.1),
    ← root2Matrix_reference rho mu sigma z k hr (ne_of_gt h.1),
    appendConstraint_kernel, root2_kernel rho mu sigma z k hr hm h]

theorem root2_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root2Matrix rho mu sigma z k) = 1 := by
  rw [matrixNullity, root2_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 2 h).1

theorem root2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : matrixNullity (root2Stack rho mu sigma z k) = 0 := by
  rw [matrixNullity, root2_stack_kernel rho mu sigma z k hr hm h]
  exact (genericChart_mode_dimensions sigma k 2 h).2

theorem root2_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root2Matrix rho mu sigma z k).rank = 2 := by
  exact matrix_rank_of_nullity (root2Matrix rho mu sigma z k) (root2_nullity rho mu sigma z k hr hm h)

theorem root2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : (root2Stack rho mu sigma z k).rank = 3 := by
  exact matrix_rank_of_nullity (root2Stack rho mu sigma z k) (root2_transverse_nullity rho mu sigma z k hr hm h)

theorem root2_difference (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : nullityDifference (WLGenericCountReference.root2Matrix rho mu sigma z k) (WLGenericCountReference.root2Stack rho mu sigma z k) = (1 : ℤ) := by
  rw [nullityDifference, root2_nullity rho mu sigma z k hr hm h, root2_transverse_nullity rho mu sigma z k hr hm h]
  norm_num

theorem root2_basis_count (rho mu sigma z : ℝ) (k : Vec 3) : basisCount (WLGenericCountReference.root2Basis rho mu sigma z k) = (1 : ℕ) := by
  rfl

theorem root2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) : basisCountResidual (WLGenericCountReference.root2Basis rho mu sigma z k) (WLGenericCountReference.root2Matrix rho mu sigma z k) = (0 : ℤ) := by
  rw [basisCountResidual, root2_basis_count, root2_nullity rho mu sigma z k hr hm h]
  norm_num

end WLGenericCountReference

end
end S10Audit.CAS
