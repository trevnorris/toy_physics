import S10Audit.BasisTrees

/-! Dimensionless coordinate and one-dimensional bases. Normalization is only
used on the stated nonzero-denominator domains, preserving the existing spans. -/

namespace S10Audit
open S10Pilot S10Anisotropic
noncomputable section
variable {D n : ℕ}

def unitVectorTree (j r : Fin (n + 1)) : Tree (n + 1) :=
  .scalar (if r = j then 1 else 0)

def longitudinalBasisTree (p r : Fin (n + 1)) : Tree (n + 1) :=
  .div (.atom (.wavevector r)) (.atom (.wavevector p))

def perpKTree (e : Fin (n + 1)) : Tree (n + 1) :=
  .sub normKTree (.pow (.atom (.wavevector e)) 2)

def extraNumeratorTree (e : Fin (n + 1)) : Tree (n + 1) :=
  .add (perpKTree e) (.mul (.atom (.coefficient .sigma)) (.pow (.atom (.wavevector e)) 2))

def extraVectorTree (e r : Fin (n + 1)) : Tree (n + 1) :=
  .sub (.mul (extraNumeratorTree e) (unitVectorTree e r))
    (.mul (.mul (.atom (.coefficient .sigma)) (.atom (.wavevector e))) (.atom (.wavevector r)))

def extraBasisTree (e r : Fin (n + 1)) : Tree (n + 1) :=
  .div (extraVectorTree e r) (perpKTree e)

theorem unit_reconstruct (a : Vec D) : a = ∑ j, a j • unit j := by
  funext r
  simp [unit, Pi.single_apply, Finset.sum_apply, mul_ite]

theorem unit_linearIndependent : LinearIndependent ℝ (unit (D := D)) := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  rw [← unit_reconstruct c] at hc
  exact congrFun hc j

theorem unitVectorTree_eval (v : Atom (n + 1) → ℝ) (j r : Fin (n + 1)) :
    Expr.eval v (unitVectorTree j r) = unit j r := by
  simp [unitVectorTree, Expr.eval, unit, Pi.single_apply, eq_comm]

theorem longitudinalBasisTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1))
    (J : Jet (n + 1)) (p r : Fin (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (longitudinalBasisTree p r) =
      (k p)⁻¹ * k r := by
  simp only [longitudinalBasisTree, Expr.eval, atomValues, div_eq_mul_inv, mul_comm]

theorem extraBasisTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1))
    (J : Jet (n + 1)) (e r : Fin (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (extraBasisTree e r) =
      (perpSq e k)⁻¹ * extraVector e sigma k r := by
  simp only [extraBasisTree, extraVectorTree, extraNumeratorTree, perpKTree,
    Expr.eval, normKTree_eval, unitVectorTree_eval, atomValues, extraVector,
    extraNumerator, perpSq, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, div_eq_mul_inv, mul_comm]

theorem unitVectorTree_hasDim (u : UnitAssignment) (j r : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (unitVectorTree j r) 1 := .scalar _

theorem longitudinalBasisTree_hasDim (u : UnitAssignment) (p r : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (longitudinalBasisTree p r) 1 := by
  simpa only [longitudinalBasisTree, atomUnits, div_self'] using
    (Expr.HasDim.div (Expr.HasDim.atom (u := atomUnits u) (.wavevector r))
      (Expr.HasDim.atom (u := atomUnits u) (.wavevector p)))

theorem extraBasisTree_hasDim (field : Dim) (e r : Fin (n + 1)) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (extraBasisTree e r) 1 := by
  let u := solvedUnits (n + 1) field
  have hp : Expr.HasDim (atomUnits u) (perpKTree e) (lengthDim⁻¹ ^ (2 : ℕ)) :=
    .sub (normKTree_hasDim u) (.pow 2 (.atom _))
  have hn : Expr.HasDim (atomUnits u) (extraNumeratorTree e) (lengthDim⁻¹ ^ (2 : ℕ)) := by
    apply Expr.HasDim.add hp
    simpa only [atomUnits, u, solvedUnits, one_mul] using
      (Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.coefficient .sigma))
        (Expr.HasDim.pow 2 (Expr.HasDim.atom (u := atomUnits u) (.wavevector e))))
  have hv : Expr.HasDim (atomUnits u) (extraVectorTree e r) (lengthDim⁻¹ ^ (2 : ℕ)) := by
    apply Expr.HasDim.sub
    · simpa only [mul_one] using Expr.HasDim.mul hn (unitVectorTree_hasDim u e r)
    · simpa only [atomUnits, u, solvedUnits, one_mul, pow_two] using
        (Expr.HasDim.mul
          (Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.coefficient .sigma))
            (Expr.HasDim.atom (u := atomUnits u) (.wavevector e)))
          (Expr.HasDim.atom (u := atomUnits u) (.wavevector r)))
  simpa only [extraBasisTree, div_self'] using Expr.HasDim.div hv hp

theorem normalized_span (a : Vec D) (c : ℝ) (hc : c ≠ 0) :
    Submodule.span ℝ {c⁻¹ • a} = Submodule.span ℝ {a} :=
  Submodule.span_singleton_smul_eq (isUnit_iff_ne_zero.mpr (inv_ne_zero hc)) a

theorem normalized_linearIndependent (a : Vec D) (ha : a ≠ 0) (c : ℝ) (hc : c ≠ 0) :
    LinearIndependent ℝ (fun _ : Fin 1 => c⁻¹ • a) := by
  rw [linearIndependent_unique_iff]
  exact smul_ne_zero (inv_ne_zero hc) ha

theorem longitudinalBasis_span (k : Vec D) (p : Fin D) (hp : k p ≠ 0) :
    Submodule.span ℝ {(k p)⁻¹ • k} = longitudinalSpace k := normalized_span k (k p) hp

theorem extraBasis_span (e : Fin D) (sigma : ℝ) (k : Vec D) (hq : perpSq e k ≠ 0) :
    Submodule.span ℝ {(perpSq e k)⁻¹ • extraVector e sigma k} =
      Submodule.span ℝ {extraVector e sigma k} := normalized_span _ _ hq

theorem longitudinalBasis_linearIndependent (k : Vec D) (p : Fin D) (hp : k p ≠ 0) :
    LinearIndependent ℝ (fun _ : Fin 1 => (k p)⁻¹ • k) := by
  apply normalized_linearIndependent k _ _ hp
  intro h
  exact hp (congrFun h p)

theorem extraBasis_linearIndependent (e : Fin D) (sigma : ℝ) (k : Vec D)
    (hq : perpSq e k ≠ 0) :
    LinearIndependent ℝ (fun _ : Fin 1 => (perpSq e k)⁻¹ • extraVector e sigma k) :=
  normalized_linearIndependent _ (extraVector_ne_zero hq) _ hq

end
end S10Audit
