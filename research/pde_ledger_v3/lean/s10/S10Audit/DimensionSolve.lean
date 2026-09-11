import S10Audit.ActionTrees

/-! Coefficient dimensions are solved from the actual action-term trees.
The displacement unit stays a free premise until the final L/T/M specialization. -/

namespace S10Audit
noncomputable section
variable {n : ℕ}

def inferredRho (D : ℕ) (field : Dim) : Dim := densityDim D / (field / timeDim) ^ (2 : ℕ)
def inferredMu (D : ℕ) (field : Dim) : Dim := densityDim D / (field / lengthDim) ^ (2 : ℕ)

def HomogeneousActionTerms (p : Package) (e : Fin (n + 1)) (u : UnitAssignment) : Prop :=
  (∀ i, Expr.infer (atomUnits u) (kineticActionTermTree p e i) = some (densityDim (n + 1))) ∧
    Expr.infer (atomUnits u) (stiffnessActionTermTree (n := n) p) = some (densityDim (n + 1))

theorem action_term_equations (p : Package) (e : Fin (n + 1)) (u : UnitAssignment) :
    HomogeneousActionTerms p e u ↔
      (∀ i, u.rho * ((if p = .anisotropic ∧ i = e then u.sigma else 1) *
        (u.field / timeDim) ^ (2 : ℕ)) = densityDim (n + 1)) ∧
      (if p = .coefficientScale then u.scale * u.mu else u.mu) *
        (u.field / lengthDim) ^ (2 : ℕ) = densityDim (n + 1) := by
  simp only [HomogeneousActionTerms, kineticActionTermTree_infer,
    stiffnessActionTermTree_infer, Option.some.injEq]

theorem dimension_product_solve (a b c : Dim) : a * b = c ↔ a = c / b := by
  constructor
  · intro h; rw [← h]; simp
  · rintro rfl; simp

theorem isotropic_dimension_solve (p : Package) (e : Fin (n + 1)) (u : UnitAssignment)
    (hp : p ≠ .anisotropic) :
    HomogeneousActionTerms p e u ↔
      u.rho = inferredRho (n + 1) u.field ∧
      (if p = .coefficientScale then u.scale * u.mu else u.mu) = inferredMu (n + 1) u.field := by
  rw [action_term_equations]
  simp only [hp, false_and, if_false, one_mul, forall_const, dimension_product_solve,
    inferredRho, inferredMu]

/-- A second inertia component determines sigma; no dimensionless premise for it is used. -/
theorem anisotropic_dimension_solve (e : Fin (n + 2)) (u : UnitAssignment) :
    HomogeneousActionTerms .anisotropic e u ↔
      u.rho = inferredRho (n + 2) u.field ∧
      u.mu = inferredMu (n + 2) u.field ∧ u.sigma = 1 := by
  rw [action_term_equations]
  simp only [true_and]
  obtain ⟨j, hj⟩ : ∃ j : Fin (n + 2), j ≠ e := by
    by_cases he : e = 0
    · exact ⟨⟨1, by omega⟩, by rw [he]; exact Fin.ne_of_val_ne (by simp)⟩
    · exact ⟨0, Ne.symm he⟩
  constructor
  · rintro ⟨hk, hs⟩
    have hr := hk j
    simp only [hj, if_false, one_mul] at hr
    have he := hk e
    simp only [ite_true] at he
    have hu : u.sigma = 1 := by
      apply mul_right_cancel (b := (u.field / timeDim) ^ (2 : ℕ))
      apply mul_left_cancel (a := u.rho)
      simpa only [one_mul] using he.trans hr.symm
    exact ⟨(dimension_product_solve _ _ _).mp hr,
      (dimension_product_solve _ _ _).mp hs, hu⟩
  · rintro ⟨hr, hm, hsig⟩
    simp [hr, hm, hsig, inferredRho, inferredMu]

/-- Without the declared scale unit, only scale*mu is fixed. Every choice is retained. -/
theorem coefficient_scale_free_family (e : Fin (n + 1)) (field scale : Dim) :
    HomogeneousActionTerms .coefficientScale e
      ⟨inferredRho (n + 1) field, inferredMu (n + 1) field / scale, 1, scale, field⟩ := by
  rw [isotropic_dimension_solve _ _ _ (by decide)]
  simp

theorem coefficient_scale_declared_dimensionless (e : Fin (n + 1)) (u : UnitAssignment)
    (hs : u.scale = 1) :
    HomogeneousActionTerms .coefficientScale e u ↔
      u.rho = inferredRho (n + 1) u.field ∧ u.mu = inferredMu (n + 1) u.field := by
  rw [isotropic_dimension_solve _ _ _ (by decide)]
  simp [hs]

def solvedUnits (D : ℕ) (field : Dim) : UnitAssignment :=
  ⟨inferredRho D field, inferredMu D field, 1, 1, field⟩

theorem solved_action_terms (p : Package) (e : Fin (n + 1)) (field : Dim) :
    HomogeneousActionTerms p e (solvedUnits (n + 1) field) := by
  rw [action_term_equations]
  simp [solvedUnits, inferredRho, inferredMu]

/-- This green action check is a consequence of the solve, for every field-unit premise. -/
theorem solved_action_infer (p : Package) (e : Fin (n + 1)) (field : Dim) :
    Expr.infer (atomUnits (solvedUnits (n + 1) field)) (actionTree p e) =
      some (densityDim (n + 1)) := by
  have h := solved_action_terms p e field
  exact (Expr.HasDim.sub (Expr.HasDim.sum fun i => Expr.infer_sound (h.1 i))
    (Expr.infer_sound h.2)).infer

theorem inferred_ratio_field_blind (D : ℕ) (field : Dim) :
    inferredMu D field / inferredRho D field = (lengthDim / timeDim) ^ (2 : ℕ) := by
  unfold inferredMu inferredRho
  ext b
  simp only [Dimension.div_exponent, Dimension.npow_exponent, nsmul_eq_mul]
  ring

theorem inferred_rho_LTM (D : ℕ) : inferredRho D lengthDim = dimensions (-(D : ℚ)) 0 1 := by
  have hn : ((D : Dimension.Exponent) : ℚ) = (D : ℚ) :=
    map_natCast Dimension.Exponent.ringEquivRat D
  ext i
  apply Dimension.Exponent.ringEquivRat.injective
  fin_cases i <;>
    simp [inferredRho, densityDim, energyDim, lengthDim, timeDim, dimensions,
      Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat, nsmul_eq_mul, hn]

theorem inferred_mu_LTM (D : ℕ) : inferredMu D lengthDim = dimensions (2 - (D : ℚ)) (-2) 1 := by
  have hn : ((D : Dimension.Exponent) : ℚ) = (D : ℚ) :=
    map_natCast Dimension.Exponent.ringEquivRat D
  ext i
  apply Dimension.Exponent.ringEquivRat.injective
  fin_cases i <;>
    simp [inferredMu, densityDim, energyDim, lengthDim, dimensions,
      Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat, nsmul_eq_mul, hn]

end
end S10Audit
