import S10Audit.MatrixCovariance

/-! N5 and N6 residuals accept arbitrary, uniformly dimensioned vector trees.
Their units follow the complete formulas, including both wavevector factors in
the longitudinality residual. -/

namespace S10Audit
open S10Pilot
noncomputable section
variable {n : ℕ}

def mulVecTree (M : Fin (n + 1) → Fin (n + 1) → Tree (n + 1))
    (a : Fin (n + 1) → Tree (n + 1)) (i : Fin (n + 1)) : Tree (n + 1) :=
  .sum n fun j => .mul (M i j) (a j)

def dotKTree (a : Fin (n + 1) → Tree (n + 1)) : Tree (n + 1) :=
  .sum n fun i => .mul (.atom (.wavevector i)) (a i)

def longitudinalResidualTree (a : Fin (n + 1) → Tree (n + 1))
    (i : Fin (n + 1)) : Tree (n + 1) :=
  .sub (.mul normKTree (a i)) (.mul (dotKTree a) (.atom (.wavevector i)))

theorem mulVecTree_eval (v : Atom (n + 1) → ℝ)
    (M : Fin (n + 1) → Fin (n + 1) → Tree (n + 1))
    (a : Fin (n + 1) → Tree (n + 1)) (i : Fin (n + 1)) :
    Expr.eval v (mulVecTree M a i) =
      Matrix.mulVec (fun i j => Expr.eval v (M i j)) (fun j => Expr.eval v (a j)) i := rfl

theorem dotKTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1))
    (J : Jet (n + 1)) (a : Fin (n + 1) → Tree (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (dotKTree a) =
      dot k (fun i => Expr.eval (atomValues rho mu sigma scale omega u k J) (a i)) := rfl

theorem longitudinalResidualTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1))
    (J : Jet (n + 1)) (a : Fin (n + 1) → Tree (n + 1)) (i : Fin (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (longitudinalResidualTree a i) =
      normSq k * Expr.eval (atomValues rho mu sigma scale omega u k J) (a i) -
        dot k (fun j => Expr.eval (atomValues rho mu sigma scale omega u k J) (a j)) * k i := by
  simp only [longitudinalResidualTree, Expr.eval, normKTree_eval, dotKTree_eval, atomValues]

theorem mulVecTree_hasDim (u : UnitAssignment)
    (M : Fin (n + 1) → Fin (n + 1) → Tree (n + 1))
    (a : Fin (n + 1) → Tree (n + 1)) (row : Fin (n + 1) → Dim) (basis : Dim)
    (hM : ∀ i j, Expr.HasDim (atomUnits u) (M i j) (row i))
    (ha : ∀ i, Expr.HasDim (atomUnits u) (a i) basis) (i : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (mulVecTree M a i) (row i * basis) :=
  .sum fun j => .mul (hM i j) (ha j)

theorem n5_hasDim (p : Package) (e i : Fin (n + 1)) (field : Dim) {z : Tree (n + 1)}
    (hz : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) z (timeDim⁻¹ ^ (2 : ℕ))) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field))
      (mulVecTree (matrixTree p e z) (fun j => .atom (.wavevector j)) i)
      (matrixDim (n + 1) field * lengthDim⁻¹) :=
  mulVecTree_hasDim _ _ _ _ lengthDim⁻¹ (fun i j => matrixTree_hasDim p e i j field hz)
    (fun j => Expr.HasDim.atom (u := atomUnits (solvedUnits (n + 1) field)) (Atom.wavevector j)) i

theorem dotKTree_hasDim (u : UnitAssignment) (a : Fin (n + 1) → Tree (n + 1)) (basis : Dim)
    (ha : ∀ i, Expr.HasDim (atomUnits u) (a i) basis) :
    Expr.HasDim (atomUnits u) (dotKTree a) (lengthDim⁻¹ * basis) :=
  .sum fun i => .mul (.atom _) (ha i)

theorem longitudinalResidualTree_hasDim (u : UnitAssignment)
    (a : Fin (n + 1) → Tree (n + 1)) (basis : Dim)
    (ha : ∀ i, Expr.HasDim (atomUnits u) (a i) basis) (i : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (longitudinalResidualTree a i) (lengthDim⁻¹ ^ (2 : ℕ) * basis) := by
  apply Expr.HasDim.sub (.mul (normKTree_hasDim u) (ha i))
  simpa only [atomUnits, pow_two, mul_assoc, mul_left_comm, mul_comm] using
    (Expr.HasDim.mul (dotKTree_hasDim u a basis ha)
      (Expr.HasDim.atom (u := atomUnits u) (.wavevector i)))

end
end S10Audit
