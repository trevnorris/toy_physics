import S10Audit.RootDimensions
import S10Audit.CoefficientInventory
import S10Audit.Curl

/-! Whole-expression unit checks and the Q7 operand connection. -/

namespace S10Audit
open S10Pilot
noncomputable section
variable {n : ℕ}

theorem action_unit_covariance (p : Package) (e : Fin (n + 1)) (field : Dim)
    (χ : Dim →* ℝˣ) (v : Atom (n + 1) → ℝ) :
    Expr.eval (fun a => (χ (atomUnits (solvedUnits (n + 1) field) a) : ℝ) * v a)
      (actionTree p e) = (χ (densityDim (n + 1)) : ℝ) * Expr.eval v (actionTree p e) :=
  (Expr.infer_sound (solved_action_infer p e field)).rescale χ v

theorem dimensionless_factor_is_invisible (p : Package) (e : Fin (n + 1)) (field : Dim) (c : ℝ) :
    Expr.infer (atomUnits (solvedUnits (n + 1) field)) (.mul (.scalar c) (actionTree p e)) =
      some (densityDim (n + 1)) := by
  simpa only [one_mul] using
    (Expr.HasDim.mul (.scalar c) (Expr.infer_sound (solved_action_infer p e field))).infer

theorem bare_field_is_counted (i : Fin (n + 1)) :
    Expr.infer (atomUnits (solvedUnits (n + 1) lengthDim)) (.pow (.atom (.field i)) 2) ≠ some 1 := by
  rw [(bare_field_square_hasDim (solvedUnits (n + 1) lengthDim) i).infer]
  intro h
  have hx := congrArg (fun d : Dim => (d.exponent 0 : ℚ)) (Option.some.inj h)
  norm_num [solvedUnits, lengthDim, dimensions, Dimension.npow_exponent] at hx

def curlSquaredTree : Tree 3 := .sum 2 fun i => .pow (.sum 2 fun j => .sum 2 fun k =>
  .mul (.scalar (leviCivitaSymbol ![i, j, k] : ℝ)) (.atom (.jet j.succ k))) 2

def q7DifferenceTree (p : Package) : Tree 3 := .sub (stiffnessTree p) curlSquaredTree

theorem curlSquaredTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec 3) (J : Jet 3) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) curlSquaredTree = curlSquared J := by
  simp only [curlSquaredTree, Expr.eval, atomValues, curlSquared, normSq, dot, epsilonCurl, pow_two]

theorem q7DifferenceTree_eval (p : Package) (rho mu sigma scale omega : ℝ)
    (u k : Vec 3) (J : Jet 3) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (q7DifferenceTree p) = q7Difference p J := by
  simp only [q7DifferenceTree, Expr.eval, stiffnessTree_eval, curlSquaredTree_eval, q7Difference]

theorem q7DifferenceTree_hasDim (p : Package) (u : UnitAssignment) :
    Expr.HasDim (atomUnits u) (q7DifferenceTree p) ((u.field / lengthDim) ^ (2 : ℕ)) := by
  apply Expr.HasDim.sub (stiffnessTree_hasDim p u)
  apply Expr.HasDim.sum
  intro i
  apply Expr.HasDim.pow
  apply Expr.HasDim.sum
  intro j
  apply Expr.HasDim.sum
  intro k
  simpa only [one_mul, atomUnits, Fin.cases_succ] using
    (Expr.HasDim.mul (.scalar (leviCivitaSymbol ![i, j, k] : ℝ))
      (Expr.HasDim.atom (u := atomUnits u) (.jet j.succ k)))

end
end S10Audit
