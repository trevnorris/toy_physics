import S10Audit.CAS.PY
import S10Audit.CAS.WL
import S10Audit.CAS.BasisCompletion

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PY

def matrixA (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n13, n18, n20], ![n18, n24, n27], ![n20, n27, n29]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixA_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixA rho mu sigma z k = ((-1) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_a_cell0 rho mu sigma z k
  · exact q2_matrix_a_cell1 rho mu sigma z k
  · exact q2_matrix_a_cell2 rho mu sigma z k
  · exact q2_matrix_a_cell3 rho mu sigma z k
  · exact q2_matrix_a_cell4 rho mu sigma z k
  · exact q2_matrix_a_cell5 rho mu sigma z k
  · exact q2_matrix_a_cell6 rho mu sigma z k
  · exact q2_matrix_a_cell7 rho mu sigma z k
  · exact q2_matrix_a_cell8 rho mu sigma z k

def matrixB (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n37, n40, n43], ![n40, n49, n52], ![n43, n52, n55]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixB_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixB rho mu sigma z k = ((1/2) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_b_cell0 rho mu sigma z k
  · exact q2_matrix_b_cell1 rho mu sigma z k
  · exact q2_matrix_b_cell2 rho mu sigma z k
  · exact q2_matrix_b_cell3 rho mu sigma z k
  · exact q2_matrix_b_cell4 rho mu sigma z k
  · exact q2_matrix_b_cell5 rho mu sigma z k
  · exact q2_matrix_b_cell6 rho mu sigma z k
  · exact q2_matrix_b_cell7 rho mu sigma z k
  · exact q2_matrix_b_cell8 rho mu sigma z k

def matrixResidual (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n68, n73, n76], ![n73, n82, n86], ![n76, n86, n88]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixResidual_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixResidual rho mu sigma z k = ((-3/2) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_residual_cell0 rho mu sigma z k
  · exact q2_matrix_residual_cell1 rho mu sigma z k
  · exact q2_matrix_residual_cell2 rho mu sigma z k
  · exact q2_matrix_residual_cell3 rho mu sigma z k
  · exact q2_matrix_residual_cell4 rho mu sigma z k
  · exact q2_matrix_residual_cell5 rho mu sigma z k
  · exact q2_matrix_residual_cell6 rho mu sigma z k
  · exact q2_matrix_residual_cell7 rho mu sigma z k
  · exact q2_matrix_residual_cell8 rho mu sigma z k

theorem routes (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixA rho mu sigma z k = (-2 : ℝ) • matrixB rho mu sigma z k := by
  rw [matrixA_reference, matrixB_reference, smul_smul]
  norm_num

theorem matrixB_action (rho mu sigma omega : ℝ) (k : Vec 3) :
    matrixB rho mu sigma (omega ^ 2) k =
      (1/2 : ℝ) • actionMatrix .anisotropic 0 rho mu sigma 1 omega k := by
  rw [matrixB_reference, referenceMatrix_action]

theorem matrixB_kernel (rho mu sigma z : ℝ) (k a : Vec 3) :
    (matrixB rho mu sigma z k).mulVec a = 0 ↔
      (referenceMatrix rho mu sigma z k).mulVec a = 0 := by
  rw [matrixB_reference, Matrix.smul_mulVec, smul_eq_zero]
  norm_num

end PY

namespace PY

def basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 3 → Vec 3 :=
  fun r j => Expr.eval (values rho mu sigma z k)
    ((![![n192, n193, n163], ![n162, n216, n163], ![n240, n193, n163]] : Fin 3 → Fin 3 → Expr Symbol) r j)

theorem basis_reference (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : basis rho mu sigma z k r = referenceBasis sigma k r := by
  rcases h with ⟨hs, _hs1, h0, h1, h2⟩
  have hsn : sigma ≠ 0 := ne_of_gt hs
  ext j
  fin_cases r <;> fin_cases j

  · exact root1_n6_nullspace_basis_cell0 rho mu sigma z k h2
  · exact root1_n6_nullspace_basis_cell1 rho mu sigma z k h2
  · exact root1_n6_nullspace_basis_cell2 rho mu sigma z k h2
  · exact root2_n6_nullspace_basis_cell0 rho mu sigma z k h1
  · exact root2_n6_nullspace_basis_cell1 rho mu sigma z k h1
  · exact root2_n6_nullspace_basis_cell2 rho mu sigma z k h1
  · exact root3_n6_nullspace_basis_cell0 rho mu sigma z k h0 h2 hsn
  · exact root3_n6_nullspace_basis_cell1 rho mu sigma z k h2 hsn h0
  · exact root3_n6_nullspace_basis_cell2 rho mu sigma z k hsn h0 h2

theorem basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) :
    LinearIndependent ℝ (fun _ : Fin 1 => basis rho mu sigma z k r) := by
  rw [basis_reference rho mu sigma z k r h]
  exact referenceBasis_independent sigma k r (genericChart_basisDomain sigma k r h)

theorem basis_complete (rho mu sigma : ℝ) (k a : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    (matrixB rho mu sigma (referenceRoot rho mu sigma k r) k).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ {basis rho mu sigma (referenceRoot rho mu sigma k r) k r} := by
  rw [matrixB_kernel, basis_reference rho mu sigma _ k r h]
  exact referenceMatrix_basis_complete rho mu sigma k a r hr hm h

end PY

namespace WL

def matrixA (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n12, n17, n20], ![n17, n24, n27], ![n20, n27, n30]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixA_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixA rho mu sigma z k = ((-1) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_a_cell0 rho mu sigma z k
  · exact q2_matrix_a_cell1 rho mu sigma z k
  · exact q2_matrix_a_cell2 rho mu sigma z k
  · exact q2_matrix_a_cell3 rho mu sigma z k
  · exact q2_matrix_a_cell4 rho mu sigma z k
  · exact q2_matrix_a_cell5 rho mu sigma z k
  · exact q2_matrix_a_cell6 rho mu sigma z k
  · exact q2_matrix_a_cell7 rho mu sigma z k
  · exact q2_matrix_a_cell8 rho mu sigma z k

def matrixB (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n34, n35, n36], ![n35, n39, n40], ![n36, n40, n43]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixB_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixB rho mu sigma z k = ((1/2) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_b_cell0 rho mu sigma z k
  · exact q2_matrix_b_cell1 rho mu sigma z k
  · exact q2_matrix_b_cell2 rho mu sigma z k
  · exact q2_matrix_b_cell3 rho mu sigma z k
  · exact q2_matrix_b_cell4 rho mu sigma z k
  · exact q2_matrix_b_cell5 rho mu sigma z k
  · exact q2_matrix_b_cell6 rho mu sigma z k
  · exact q2_matrix_b_cell7 rho mu sigma z k
  · exact q2_matrix_b_cell8 rho mu sigma z k

def matrixResidual (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    ((![![n46, n51, n54], ![n51, n56, n60], ![n54, n60, n62]] : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem matrixResidual_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixResidual rho mu sigma z k = ((-3/2) : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
  · exact q2_matrix_residual_cell0 rho mu sigma z k
  · exact q2_matrix_residual_cell1 rho mu sigma z k
  · exact q2_matrix_residual_cell2 rho mu sigma z k
  · exact q2_matrix_residual_cell3 rho mu sigma z k
  · exact q2_matrix_residual_cell4 rho mu sigma z k
  · exact q2_matrix_residual_cell5 rho mu sigma z k
  · exact q2_matrix_residual_cell6 rho mu sigma z k
  · exact q2_matrix_residual_cell7 rho mu sigma z k
  · exact q2_matrix_residual_cell8 rho mu sigma z k

theorem routes (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixA rho mu sigma z k = (-2 : ℝ) • matrixB rho mu sigma z k := by
  rw [matrixA_reference, matrixB_reference, smul_smul]
  norm_num

theorem matrixB_action (rho mu sigma omega : ℝ) (k : Vec 3) :
    matrixB rho mu sigma (omega ^ 2) k =
      (1/2 : ℝ) • actionMatrix .anisotropic 0 rho mu sigma 1 omega k := by
  rw [matrixB_reference, referenceMatrix_action]

theorem matrixB_kernel (rho mu sigma z : ℝ) (k a : Vec 3) :
    (matrixB rho mu sigma z k).mulVec a = 0 ↔
      (referenceMatrix rho mu sigma z k).mulVec a = 0 := by
  rw [matrixB_reference, Matrix.smul_mulVec, smul_eq_zero]
  norm_num

end WL

namespace WL

def basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 3 → Vec 3 :=
  fun r j => Expr.eval (values rho mu sigma z k)
    ((![![n102, n103, n79], ![n78, n128, n79], ![n152, n103, n79]] : Fin 3 → Fin 3 → Expr Symbol) r j)

theorem basis_reference (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : basis rho mu sigma z k r = referenceBasis sigma k r := by
  rcases h with ⟨hs, _hs1, h0, h1, h2⟩
  have hsn : sigma ≠ 0 := ne_of_gt hs
  ext j
  fin_cases r <;> fin_cases j

  · exact root1_n6_nullspace_basis_cell0 rho mu sigma z k h2
  · exact root1_n6_nullspace_basis_cell1 rho mu sigma z k h2
  · exact root1_n6_nullspace_basis_cell2 rho mu sigma z k h2
  · exact root2_n6_nullspace_basis_cell0 rho mu sigma z k h1
  · exact root2_n6_nullspace_basis_cell1 rho mu sigma z k h1
  · exact root2_n6_nullspace_basis_cell2 rho mu sigma z k h1
  · exact root3_n6_nullspace_basis_cell0 rho mu sigma z k h0 h2 hsn
  · exact root3_n6_nullspace_basis_cell1 rho mu sigma z k h2 hsn h0
  · exact root3_n6_nullspace_basis_cell2 rho mu sigma z k hsn h0 h2

theorem basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) :
    LinearIndependent ℝ (fun _ : Fin 1 => basis rho mu sigma z k r) := by
  rw [basis_reference rho mu sigma z k r h]
  exact referenceBasis_independent sigma k r (genericChart_basisDomain sigma k r h)

theorem basis_complete (rho mu sigma : ℝ) (k a : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    (matrixB rho mu sigma (referenceRoot rho mu sigma k r) k).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ {basis rho mu sigma (referenceRoot rho mu sigma k r) k r} := by
  rw [matrixB_kernel, basis_reference rho mu sigma _ k r h]
  exact referenceMatrix_basis_complete rho mu sigma k a r hr hm h

end WL

theorem matrixA_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) :
    PY.matrixA rho mu sigma z k = WL.matrixA rho mu sigma z k := by
  rw [PY.matrixA_reference, WL.matrixA_reference]

theorem matrixB_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) :
    PY.matrixB rho mu sigma z k = WL.matrixB rho mu sigma z k := by
  rw [PY.matrixB_reference, WL.matrixB_reference]

theorem matrixResidual_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) :
    PY.matrixResidual rho mu sigma z k = WL.matrixResidual rho mu sigma z k := by
  rw [PY.matrixResidual_reference, WL.matrixResidual_reference]

theorem basis_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : PY.basis rho mu sigma z k r = WL.basis rho mu sigma z k r := by
  rw [PY.basis_reference rho mu sigma z k r h, WL.basis_reference rho mu sigma z k r h]

end
end S10Audit.CAS
