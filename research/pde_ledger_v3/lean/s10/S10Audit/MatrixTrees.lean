import S10Audit.Checks
import Mathlib.LinearAlgebra.Matrix.Rank

/-! Entrywise trees for the action-derived modal matrices. The squared-frequency
tree is a parameter, so ordinary, negative, anisotropic and static roots share
the same construction. Zero entries retain the units of their matrix slot. -/

namespace S10Audit
open S10Pilot
noncomputable section
variable {n : ℕ}

def stiffnessMatrixTree (p : Package) (i j : Fin (n + 1)) : Tree (n + 1) :=
  let diagonal := Expr.mul (.scalar (if i = j then 1 else 0)) normKTree
  let outer := Expr.mul (.atom (.wavevector i)) (.atom (.wavevector j))
  match p with
  | .fullGradient => diagonal
  | .divergenceOnly => outer
  | _ => .sub diagonal outer

def matrixTree (p : Package) (e : Fin (n + 1)) (z : Tree (n + 1))
    (i j : Fin (n + 1)) : Tree (n + 1) :=
  .sub (.mul (.mul (.atom (.coefficient .rho)) z)
    (.mul (if p = .anisotropic ∧ i = e then .atom (.coefficient .sigma) else .scalar 1)
      (.scalar (if i = j then 1 else 0))))
    (.mul (coefficientTree p) (stiffnessMatrixTree p i j))

def packageOperator (p : Package) {D : ℕ} (e : Fin D) (rho mu sigma scale omega : ℝ)
    (k a : Vec D) : Vec D := match p with
  | .main => S10Pilot.modalOperator rho mu omega k a
  | .fullGradient => S10Controls.modalOperator .fullGradient rho mu omega k a
  | .divergenceOnly => S10Controls.modalOperator .divergenceOnly rho mu omega k a
  | .signFlip => S10ScalarControls.spectralOperator (-1) rho mu (omega ^ 2) k a
  | .anisotropic => S10Anisotropic.modalOperator e sigma rho mu omega k a
  | .coefficientScale => S10ScalarControls.spectralOperator scale rho mu (omega ^ 2) k a

def actionMatrix (p : Package) {D : ℕ} (e : Fin D) (rho mu sigma scale omega : ℝ)
    (k : Vec D) : Matrix (Fin D) (Fin D) ℝ :=
  fun i j => packageOperator p e rho mu sigma scale omega k (S10Anisotropic.unit j) i

theorem matrixTree_eval (p : Package) (e i j : Fin (n + 1))
    (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J)
      (matrixTree p e (.pow (.atom .frequency) 2) i j) =
      actionMatrix p e rho mu sigma scale omega k i j := by
  cases p <;>
    simp only [matrixTree, stiffnessMatrixTree, Expr.eval, atomValues, normKTree_eval,
      coefficientTree, actionMatrix, packageOperator, S10Anisotropic.modalOperator,
      S10Pilot.modalOperator, S10Controls.modalOperator, S10ScalarControls.spectralOperator,
      Pi.add_apply, Pi.smul_apply, smul_eq_mul, S10Anisotropic.dot_unit_right]
  all_goals
    by_cases hij : i = j <;> by_cases hie : i = e <;> by_cases hje : j = e <;>
      simp_all [S10Anisotropic.unit, Expr.eval, atomValues, eq_comm] <;> ring

theorem actionMatrix_mulVec (p : Package) (e : Fin (n + 1))
    (rho mu sigma scale omega : ℝ) (k a : Vec (n + 1)) :
    (actionMatrix p e rho mu sigma scale omega k).mulVec a =
      packageOperator p e rho mu sigma scale omega k a := by
  ext i
  simp only [Matrix.mulVec, dotProduct, actionMatrix]
  cases p <;>
    simp only [packageOperator, S10Anisotropic.modalOperator, S10Pilot.modalOperator,
      S10Controls.modalOperator, S10ScalarControls.spectralOperator,
      Pi.add_apply, Pi.smul_apply, smul_eq_mul, S10Anisotropic.dot_unit_right]
  all_goals
    simp only [S10Anisotropic.unit, Pi.single_apply]
    simp only [sub_mul, add_mul, mul_ite, ite_mul, mul_zero, mul_one,
      Finset.sum_sub_distrib, Finset.sum_add_distrib]
    try simp only [Finset.sum_ite_eq, Finset.sum_ite_eq', Finset.mem_univ, if_true]
    simp only [mul_assoc, ← Finset.mul_sum, S10Pilot.dot]
    simp [sub_mul, ite_mul, Finset.sum_sub_distrib, Finset.sum_ite_irrel,
      mul_assoc, mul_left_comm, ← Finset.mul_sum] <;> ring_nf <;> simp

def matrixDim (D : ℕ) (field : Dim) : Dim := inferredRho D field * timeDim⁻¹ ^ (2 : ℕ)

theorem matrixDim_eq (D : ℕ) (field : Dim) : matrixDim D field = densityDim D / field ^ (2 : ℕ) := by
  ext b
  simp only [matrixDim, inferredRho, Dimension.mul_exponent, Dimension.div_exponent,
    Dimension.inv_exponent, Dimension.npow_exponent, nsmul_eq_mul]
  ring

theorem coefficientTree_hasDim (p : Package) (field : Dim) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (coefficientTree (n := n) p)
      (inferredMu (n + 1) field) := by
  let u := solvedUnits (n + 1) field
  cases p
  case coefficientScale => simpa only [coefficientTree, atomUnits, u, solvedUnits, one_mul] using
    (Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.coefficient .scale))
      (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)))
  case signFlip => simpa only [coefficientTree, atomUnits, u, solvedUnits, one_mul] using
    (Expr.HasDim.mul (Expr.HasDim.scalar (-1))
      (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)))
  all_goals exact .atom _

theorem stiffnessMatrixTree_hasDim (p : Package) (i j : Fin (n + 1)) (u : UnitAssignment) :
    Expr.HasDim (atomUnits u) (stiffnessMatrixTree p i j) (lengthDim⁻¹ ^ (2 : ℕ)) := by
  have hd := Expr.HasDim.mul (.scalar (if i = j then 1 else 0)) (normKTree_hasDim (n := n) u)
  have ho := Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.wavevector i))
    (Expr.HasDim.atom (u := atomUnits u) (.wavevector j))
  have hd' : Expr.HasDim (atomUnits (D := n + 1) u)
      (.mul (.scalar (if i = j then 1 else 0)) normKTree) (lengthDim⁻¹ ^ (2 : ℕ)) := by
    simpa only [one_mul] using hd
  have ho' : Expr.HasDim (atomUnits u) (.mul (.atom (.wavevector i)) (.atom (.wavevector j)))
      (lengthDim⁻¹ ^ (2 : ℕ)) := by simpa only [atomUnits, pow_two] using ho
  cases p
  case fullGradient => exact hd'
  case divergenceOnly => exact ho'
  all_goals exact .sub hd' ho'

theorem matrixTree_hasDim (p : Package) (e i j : Fin (n + 1)) (field : Dim)
    {z : Tree (n + 1)}
    (hz : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) z (timeDim⁻¹ ^ (2 : ℕ))) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (matrixTree p e z i j)
      (matrixDim (n + 1) field) := by
  let u := solvedUnits (n + 1) field
  have hw : Expr.HasDim (atomUnits (D := n + 1) u)
      (if p = .anisotropic ∧ i = e then .atom (.coefficient .sigma) else .scalar 1) 1 := by
    split_ifs <;> constructor
  have ht := Expr.HasDim.mul (Expr.HasDim.mul
    (Expr.HasDim.atom (u := atomUnits u) (.coefficient .rho)) hz)
    (Expr.HasDim.mul hw (.scalar (if i = j then 1 else 0)))
  have hs := Expr.HasDim.mul (coefficientTree_hasDim (n := n) p field)
    (stiffnessMatrixTree_hasDim p i j u)
  have he : inferredMu (n + 1) field * lengthDim⁻¹ ^ (2 : ℕ) = matrixDim (n + 1) field := by
    ext b
    simp only [inferredMu, matrixDim, inferredRho, Dimension.mul_exponent,
      Dimension.div_exponent, Dimension.npow_exponent, Dimension.inv_exponent, nsmul_eq_mul]
    ring
  apply Expr.HasDim.sub
  · simpa only [atomUnits, u, solvedUnits, one_mul, mul_one, matrixDim] using ht
  · exact he ▸ hs

/-- The appended wavevector row has its own unit, independent of the matrix rows. -/
def stackedTree (M : Fin (n + 1) → Fin (n + 1) → Tree (n + 1))
    (r : Fin (n + 2)) (j : Fin (n + 1)) : Tree (n + 1) :=
  Fin.lastCases (.atom (.wavevector j)) (fun i => M i j) r

def stackedRowDim (D : ℕ) (field : Dim) : Fin (D + 1) → Dim :=
  Fin.lastCases lengthDim⁻¹ (fun _ => matrixDim D field)

theorem stackedTree_hasDim (p : Package) (e : Fin (n + 1)) (field : Dim)
    {z : Tree (n + 1)}
    (hz : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) z (timeDim⁻¹ ^ (2 : ℕ)))
    (r : Fin (n + 2)) (j : Fin (n + 1)) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (stackedTree (matrixTree p e z) r j)
      (stackedRowDim (n + 1) field r) := by
  refine Fin.lastCases ?_ (fun i => ?_) r
  · simpa only [stackedTree, stackedRowDim, Fin.lastCases_last, atomUnits] using
      (Expr.HasDim.atom (u := atomUnits (solvedUnits (n + 1) field)) (.wavevector j))
  · simpa only [stackedTree, stackedRowDim, Fin.lastCases_castSucc] using
      matrixTree_hasDim p e i j field hz

def staticRootTree : Tree (n + 1) := .mul (.scalar 0) (.pow (.atom .frequency) 2)

theorem staticRootTree_eval (v : Atom (n + 1) → ℝ) : Expr.eval v staticRootTree = 0 := by
  simp [staticRootTree, Expr.eval]

theorem staticRootTree_hasDim (u : UnitAssignment) :
    Expr.HasDim (atomUnits u) (staticRootTree (n := n)) (timeDim⁻¹ ^ (2 : ℕ)) := by
  simpa only [staticRootTree, atomUnits, one_mul] using
    (Expr.HasDim.mul (.scalar 0) (Expr.HasDim.pow 2 (Expr.HasDim.atom (u := atomUnits u) .frequency)))

end
end S10Audit
