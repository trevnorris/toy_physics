import S10Audit.Packages
import S10Audit.Dimensions

/-! Expression trees are evaluated back to the existing action definitions.
Bare displacement and derivative atoms are separate; both carry field units. -/

namespace S10Audit
open S10Pilot
noncomputable section

inductive Coefficient where
  | rho | mu | sigma | scale
  deriving DecidableEq

inductive Atom (D : ℕ) where
  | coefficient (c : Coefficient)
  | field (i : Fin D)
  | jet (j : Fin (D + 1)) (i : Fin D)
  | wavevector (i : Fin D)
  | frequency
  deriving DecidableEq

structure UnitAssignment where
  rho : Dim
  mu : Dim
  sigma : Dim
  scale : Dim
  field : Dim

def atomUnits {D : ℕ} (u : UnitAssignment) : Atom D → Dim
  | .coefficient .rho => u.rho
  | .coefficient .mu => u.mu
  | .coefficient .sigma => u.sigma
  | .coefficient .scale => u.scale
  | .field _ => u.field
  | .jet j _ => u.field / (Fin.cases timeDim (fun _ => lengthDim) j)
  | .wavevector _ => lengthDim⁻¹
  | .frequency => timeDim⁻¹

def atomValues {D : ℕ} (rho mu sigma scale omega : ℝ) (u k : Vec D) (J : Jet D) :
    Atom D → ℝ
  | .coefficient .rho => rho
  | .coefficient .mu => mu
  | .coefficient .sigma => sigma
  | .coefficient .scale => scale
  | .field i => u i
  | .jet j i => J j i
  | .wavevector i => k i
  | .frequency => omega

variable {n : ℕ}
abbrev Tree (D : ℕ) := Expr (Atom D)

def stiffnessTree (p : Package) : Tree (n + 1) := match p with
  | .fullGradient => .sum n fun i => .sum n fun j => .pow (.atom (.jet i.succ j)) 2
  | .divergenceOnly => .pow (.sum n fun i => .atom (.jet i.succ i)) 2
  | _ => .mul (.scalar (1 / 2)) (.sum n fun i => .sum n fun j =>
      .pow (.sub (.atom (.jet i.succ j)) (.atom (.jet j.succ i))) 2)

def kineticTermTree (p : Package) (e i : Fin (n + 1)) : Tree (n + 1) :=
  .mul (if p = .anisotropic ∧ i = e then .atom (.coefficient .sigma) else .scalar 1)
    (.pow (.atom (.jet 0 i)) 2)

def coefficientTree (p : Package) : Tree (n + 1) := match p with
  | .coefficientScale => .mul (.atom (.coefficient .scale)) (.atom (.coefficient .mu))
  | .signFlip => .mul (.scalar (-1)) (.atom (.coefficient .mu))
  | _ => .atom (.coefficient .mu)

def kineticActionTermTree (p : Package) (e i : Fin (n + 1)) : Tree (n + 1) :=
  .mul (.div (.atom (.coefficient .rho)) (.scalar 2)) (kineticTermTree p e i)

def stiffnessActionTermTree (p : Package) : Tree (n + 1) :=
  .mul (.div (coefficientTree p) (.scalar 2)) (stiffnessTree p)

def actionTree (p : Package) (e : Fin (n + 1)) : Tree (n + 1) :=
  .sub (.sum n fun i => kineticActionTermTree p e i) (stiffnessActionTermTree p)

theorem stiffnessTree_eval (p : Package) (rho mu sigma scale omega : ℝ)
    (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (stiffnessTree p) =
      packageStiffness p J := by
  cases p <;> rfl

theorem actionTree_eval (p : Package) (e : Fin (n + 1)) (rho mu sigma scale omega : ℝ)
    (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (actionTree p e) =
      packageAction p e rho mu sigma scale J := by
  rw [packageAction_uses_stiffness]
  have hc : Expr.eval (atomValues rho mu sigma scale omega u k J) (coefficientTree p) =
      stiffnessCoefficient p mu scale := by cases p <;> simp [coefficientTree, Expr.eval,
        atomValues, stiffnessCoefficient]
  simp only [actionTree, kineticActionTermTree, stiffnessActionTermTree, Expr.eval,
    atomValues, stiffnessTree_eval, hc, ← Finset.mul_sum]
  congr 1
  congr 1
  cases p <;> simp [kineticTermTree, Expr.eval, atomValues, packageKinetic,
    S10Anisotropic.kinetic, normSq, dot, pow_two]
  apply Finset.sum_congr rfl
  intro i _
  split_ifs <;> simp [Expr.eval, atomValues]

theorem stiffnessTree_hasDim (p : Package) (u : UnitAssignment) :
    Expr.HasDim (atomUnits u) (stiffnessTree (n := n) p) ((u.field / lengthDim) ^ (2 : ℕ)) := by
  have h (i j : Fin (n + 1)) : Expr.HasDim (atomUnits u)
      (.atom (.jet i.succ j)) (u.field / lengthDim) := Expr.HasDim.atom _
  have hs : Expr.HasDim (atomUnits u) (.sum n fun i => .sum n fun j =>
      .pow (.sub (.atom (.jet i.succ j)) (.atom (.jet j.succ i))) 2)
      ((u.field / lengthDim) ^ (2 : ℕ)) :=
    .sum fun i => .sum fun j => .pow 2 (.sub (h i j) (h j i))
  cases p
  all_goals unfold stiffnessTree
  case fullGradient => exact .sum fun i => .sum fun j => .pow 2 (h i j)
  case divergenceOnly => exact .pow 2 (.sum fun i => h i i)
  all_goals simpa only [one_mul] using (Expr.HasDim.mul (.scalar (1 / 2)) hs)

theorem kineticActionTermTree_infer (p : Package) (e i : Fin (n + 1)) (u : UnitAssignment) :
    Expr.infer (atomUnits u) (kineticActionTermTree p e i) =
      some (u.rho * ((if p = .anisotropic ∧ i = e then u.sigma else 1) *
        (u.field / timeDim) ^ (2 : ℕ))) := by
  unfold kineticActionTermTree kineticTermTree
  split_ifs <;> simp [Expr.infer, atomUnits]

theorem stiffnessActionTermTree_infer (p : Package) (u : UnitAssignment) :
    Expr.infer (atomUnits u) (stiffnessActionTermTree (n := n) p) =
      some ((if p = .coefficientScale then u.scale * u.mu else u.mu) *
        (u.field / lengthDim) ^ (2 : ℕ)) := by
  have hc : Expr.HasDim (atomUnits u) (coefficientTree (n := n) p)
      (if p = .coefficientScale then u.scale * u.mu else u.mu) := by
    cases p
    case coefficientScale => exact .mul (.atom _) (.atom _)
    case signFlip => simpa only [coefficientTree, reduceCtorEq, if_false, one_mul] using
      (Expr.HasDim.mul (Expr.HasDim.scalar (-1)) (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)) :
        Expr.HasDim (atomUnits u) _ (1 * u.mu))
    all_goals exact .atom _
  simpa only [stiffnessActionTermTree, div_one] using
    (Expr.HasDim.mul (.div hc (.scalar 2)) (stiffnessTree_hasDim (n := n) p u)).infer

end
end S10Audit
