import S10Audit.DimensionSolve

/-! Root dimensions walk the full formula, including every wavevector factor.
The zero expression is deliberately not assigned a unique physical dimension. -/

namespace S10Audit
open S10Pilot
noncomputable section
variable {n : ℕ}

def normKTree : Tree (n + 1) := .sum n fun i => .pow (.atom (.wavevector i)) 2

def ordinaryRootTree (p : Package) : Tree (n + 1) :=
  .mul (.div (coefficientTree p) (.atom (.coefficient .rho))) normKTree

def extraRootTree (e : Fin (n + 1)) : Tree (n + 1) :=
  .mul (.div (.atom (.coefficient .mu)) (.atom (.coefficient .rho)))
    (.div (.add (.sub normKTree (.pow (.atom (.wavevector e)) 2))
      (.mul (.atom (.coefficient .sigma)) (.pow (.atom (.wavevector e)) 2)))
      (.atom (.coefficient .sigma)))

theorem normKTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) normKTree = normSq k := by
  simp [normKTree, Expr.eval, atomValues, normSq, dot, pow_two]

theorem ordinaryRootTree_eval (p : Package) (rho mu sigma scale omega : ℝ)
    (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (ordinaryRootTree p) =
      S10Pilot.coneValue rho (stiffnessCoefficient p mu scale) k := by
  have hc : Expr.eval (atomValues rho mu sigma scale omega u k J) (coefficientTree p) =
      stiffnessCoefficient p mu scale := by cases p <;> simp [coefficientTree,
        stiffnessCoefficient, Expr.eval, atomValues]
  simp only [ordinaryRootTree, Expr.eval, hc, atomValues, normKTree_eval, S10Pilot.coneValue]

theorem extraRootTree_eval (e : Fin (n + 1)) (rho mu sigma scale omega : ℝ)
    (u k : Vec (n + 1)) (J : Jet (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (extraRootTree e) =
      S10Anisotropic.extraConeValue e sigma rho mu k := by
  simp only [extraRootTree, Expr.eval, atomValues, normKTree_eval,
    S10Anisotropic.extraConeValue, S10Anisotropic.extraValue,
    S10Anisotropic.extraNumerator, S10Anisotropic.perpSq]

theorem normKTree_hasDim (u : UnitAssignment) :
    Expr.HasDim (atomUnits u) (normKTree (n := n)) (lengthDim⁻¹ ^ (2 : ℕ)) :=
  .sum fun _ => .pow 2 (.atom _)

theorem ordinaryRootTree_hasDim (p : Package) (field : Dim) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (ordinaryRootTree (n := n) p)
      (timeDim⁻¹ ^ (2 : ℕ)) := by
  let u := solvedUnits (n + 1) field
  have hc : Expr.HasDim (atomUnits u) (coefficientTree (n := n) p) u.mu := by
    cases p
    case coefficientScale => simpa only [coefficientTree, u, solvedUnits, one_mul] using
      (Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.coefficient .scale))
        (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)) : Expr.HasDim (atomUnits u) _ (u.scale * u.mu))
    case signFlip => simpa only [coefficientTree, reduceCtorEq, if_false, one_mul] using
      (Expr.HasDim.mul (Expr.HasDim.scalar (-1)) (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)) :
        Expr.HasDim (atomUnits u) _ (1 * u.mu))
    all_goals exact .atom _
  have h := Expr.HasDim.mul (.div hc (Expr.HasDim.atom (u := atomUnits u) (.coefficient .rho))) (normKTree_hasDim (n := n) u)
  have hd : u.mu / u.rho * (lengthDim⁻¹ ^ (2 : ℕ)) = timeDim⁻¹ ^ (2 : ℕ) := by
    dsimp [u, solvedUnits]
    rw [inferred_ratio_field_blind]
    ext b
    simp only [Dimension.mul_exponent, Dimension.div_exponent, Dimension.npow_exponent,
      Dimension.inv_exponent, nsmul_eq_mul]
    ring
  exact hd ▸ h

theorem extraRootTree_hasDim (e : Fin (n + 1)) (field : Dim) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (extraRootTree e)
      (timeDim⁻¹ ^ (2 : ℕ)) := by
  let u := solvedUnits (n + 1) field
  have he : Expr.HasDim (atomUnits u) (.pow (.atom (.wavevector e)) 2)
      (lengthDim⁻¹ ^ (2 : ℕ)) := .pow 2 (.atom _)
  have hs : Expr.HasDim (atomUnits u)
      (.mul (.atom (.coefficient .sigma)) (.pow (.atom (.wavevector e)) 2))
      (lengthDim⁻¹ ^ (2 : ℕ)) := by
    simpa only [u, solvedUnits, one_mul] using
      (Expr.HasDim.mul (Expr.HasDim.atom (u := atomUnits u) (.coefficient .sigma)) he :
        Expr.HasDim (atomUnits u) _ (u.sigma * (lengthDim⁻¹ ^ (2 : ℕ))))
  have hsum := Expr.HasDim.add (.sub (normKTree_hasDim (n := n) u) he) hs
  have h := Expr.HasDim.mul
    (Expr.HasDim.div (Expr.HasDim.atom (u := atomUnits u) (.coefficient .mu)) (Expr.HasDim.atom (u := atomUnits u) (.coefficient .rho)))
    (Expr.HasDim.div hsum (Expr.HasDim.atom (u := atomUnits u) (.coefficient .sigma)))
  have hd : u.mu / u.rho * ((lengthDim⁻¹ ^ (2 : ℕ)) / u.sigma) = timeDim⁻¹ ^ (2 : ℕ) := by
    dsimp [u, solvedUnits]
    rw [inferred_ratio_field_blind]
    ext b
    simp only [Dimension.mul_exponent, Dimension.div_exponent, Dimension.npow_exponent,
      Dimension.inv_exponent, Dimension.one_exponent, nsmul_eq_mul]
    ring
  exact hd ▸ h

theorem root_ratio_hasDim (field : Dim) {root : Tree (n + 1)}
    (h : Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) root (timeDim⁻¹ ^ (2 : ℕ))) :
    Expr.HasDim (atomUnits (solvedUnits (n + 1) field)) (.div root normKTree)
      ((lengthDim / timeDim) ^ (2 : ℕ)) := by
  have hd : (timeDim⁻¹ ^ (2 : ℕ)) / (lengthDim⁻¹ ^ (2 : ℕ)) =
      (lengthDim / timeDim) ^ (2 : ℕ) := by
    ext b
    simp only [Dimension.div_exponent, Dimension.npow_exponent, Dimension.inv_exponent,
      nsmul_eq_mul]
    ring
  exact hd ▸ Expr.HasDim.div h (normKTree_hasDim _)

theorem bare_field_square_hasDim (u : UnitAssignment) (i : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (.pow (.atom (.field i)) 2) (u.field ^ (2 : ℕ)) :=
  .pow 2 (.atom _)

theorem dropped_wavevector_changes_dimension (D : ℕ) (field : Dim) :
    inferredMu D field / inferredRho D field ≠ timeDim⁻¹ ^ (2 : ℕ) := by
  rw [inferred_ratio_field_blind]
  intro h
  have hx := congrArg (fun d : Dim => (d.exponent 0).toRat) h
  norm_num [lengthDim, timeDim, dimensions, Dimension.npow_exponent] at hx

end
end S10Audit
