import S10Audit.DimensionSolve

/-! Unknown coefficient counts are extracted from each action's expression tree.
Only the coefficient-scale declaration removes a present symbol from the solve. -/

namespace S10Audit
noncomputable section
variable {D n : ℕ}

def coefficients : Tree D → Finset Coefficient
  | .atom (.coefficient c) => {c}
  | .atom _ | .scalar _ => ∅
  | .add a b | .sub a b | .mul a b | .div a b => coefficients a ∪ coefficients b
  | .pow a _ => coefficients a
  | .sum _ f => Finset.univ.biUnion fun i => coefficients (f i)

theorem stiffnessTree_coefficients (p : Package) :
    coefficients (stiffnessTree (n := n) p) = ∅ := by
  cases p <;> ext c <;> simp [stiffnessTree, coefficients]

theorem actionTree_coefficients (p : Package) (e : Fin (n + 1)) :
    coefficients (actionTree p e) =
      if p = .anisotropic then {.rho, .mu, .sigma}
      else if p = .coefficientScale then {.rho, .mu, .scale} else {.rho, .mu} := by
  have hs := stiffnessTree_coefficients (n := n) p
  simp only [actionTree, kineticActionTermTree, stiffnessActionTermTree, coefficients, hs]
  cases p <;> ext c <;> cases c <;>
    simp [kineticTermTree, coefficientTree, coefficients]
  · exact ⟨e, by simp [coefficients]⟩
  · intro i; split_ifs <;> simp [coefficients]

def unknownCoefficients (p : Package) (e : Fin (n + 1)) : Finset Coefficient :=
  coefficients (actionTree p e) \ (if p = .coefficientScale then {.scale} else ∅)

theorem unknown_coefficient_count (p : Package) (e : Fin (n + 1)) :
    (unknownCoefficients p e).card = if p = .anisotropic then 3 else 2 := by
  rw [unknownCoefficients, actionTree_coefficients]
  cases p <;> decide +revert

theorem anisotropic_scale_is_unknown (e : Fin (n + 1)) :
    .sigma ∈ unknownCoefficients .anisotropic e := by
  simp [unknownCoefficients, actionTree_coefficients]

theorem declared_scale_is_not_unknown (e : Fin (n + 1)) :
    .scale ∉ unknownCoefficients .coefficientScale e := by
  simp [unknownCoefficients, actionTree_coefficients]

set_option backward.isDefEq.respectTransparency false in
/-- In the rho/mu/sigma order, unweighted kinetic, stiffness, and weighted kinetic
equations have this invertible coefficient block for each L/T/M coordinate. -/
theorem anisotropic_equation_block :
    Matrix.det (![![1, 0, 0], ![0, 1, 0], ![1, 0, 1]] : Matrix (Fin 3) (Fin 3) ℚ) = 1 := by
  rw [Matrix.det_fin_three]
  change (1 * 1 * 1 - 1 * 0 * 0 - 0 * 0 * 1 + 0 * 0 * 1 + 0 * 0 * 0 - 0 * 1 * 1 : ℚ) = 1
  norm_num

end
end S10Audit
