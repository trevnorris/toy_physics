import S10Audit.ResidualTrees

/-! Coordinate charts for complete transverse subspaces. Each chart has an
explicit nonzero pivot; no single coordinate is assumed nonzero globally. -/

namespace S10Audit
open S10Pilot S10Anisotropic
noncomputable section
variable {D n : ℕ}

def chartVector (k : Vec D) (p j : Fin D) : Vec D :=
  unit j - (k j / k p) • unit p

theorem chartVector_coordinate (k : Vec D) (p j r : Fin D) (hr : r ≠ p) :
    chartVector k p j r = if r = j then 1 else 0 := by
  simp [chartVector, unit, Pi.single_apply, hr, eq_comm]

theorem chartVector_transverse (k : Vec D) (p j : Fin D) (hp : k p ≠ 0) :
    dot k (chartVector k p j) = 0 := by
  change (dotLinear k) (unit j - (k j / k p) • unit p) = 0
  rw [map_sub, map_smul]
  change dot k (unit j) - (k j / k p) * dot k (unit p) = 0
  rw [dot_unit_right, dot_unit_right, div_mul_cancel₀ _ hp, sub_self]

/-- Equality on the nonpivot coordinates determines the remaining coordinate
inside the transverse subspace. -/
theorem transverse_ext (k a b : Vec D) (p : Fin D) (hp : k p ≠ 0)
    (ha : dot k a = 0) (hb : dot k b = 0) (h : ∀ r, r ≠ p → a r = b r) : a = b := by
  have hd : dot k (a - b) = 0 := by
    change (dotLinear k) (a - b) = 0
    rw [map_sub]
    exact sub_eq_zero.mpr (ha.trans hb.symm)
  have hs : dot k (a - b) = k p * (a p - b p) := by
    apply Finset.sum_eq_single p
    · intro r _ hr
      simp only [Pi.sub_apply, h r hr, sub_self, mul_zero]
    · simp
  have hpab : a p = b p := sub_eq_zero.mp ((mul_eq_zero.mp (hs.symm.trans hd)).resolve_left hp)
  funext r
  by_cases hr : r = p
  · exact hr ▸ hpab
  · exact h r hr

/-- The allowed nonpivot coordinates are arbitrary. Excluding the distinguished
axis gives the ordinary anisotropic subspace; allowing all gives the full
transverse subspace. -/
theorem chart_reconstruct (k a : Vec D) (p : Fin D) (hp : k p ≠ 0)
    (free : Finset (Fin D)) (ha : dot k a = 0)
    (hsupport : ∀ r, r ≠ p → r ∉ free → a r = 0) :
    a = ∑ j ∈ free, a j • chartVector k p j := by
  apply transverse_ext k _ _ p hp ha
  · change (dotLinear k) (∑ j ∈ free, a j • chartVector k p j) = 0
    simp only [map_sum, map_smul]
    change (∑ j ∈ free, a j * dot k (chartVector k p j)) = 0
    simp only [chartVector_transverse k p _ hp, mul_zero, Finset.sum_const_zero]
  · intro r hr
    simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, chartVector_coordinate k p _ r hr,
      mul_ite, mul_one, mul_zero]
    by_cases hrf : r ∈ free
    · simp [hrf]
    · simp [hrf, hsupport r hr hrf]

theorem chart_linearIndependent (k : Vec D) (p : Fin D) (free : Finset (Fin D))
    (hfree : p ∉ free) : LinearIndependent ℝ (fun j : free => chartVector k p j.val) := by
  rw [Fintype.linearIndependent_iff]
  intro c hc j
  have hjp : j.val ≠ p := by intro h; exact hfree (h ▸ j.property)
  have hv := congrFun hc j.val
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply,
    chartVector_coordinate k p _ j.val hjp, mul_ite, mul_one, mul_zero] at hv
  simpa only [Subtype.val_inj, Finset.sum_ite_eq, Finset.mem_univ, if_true] using hv

theorem transverse_chart_reconstruct (k a : Vec D) (p : Fin D) (hp : k p ≠ 0)
    (ha : a ∈ transverseSpace k) :
    a = ∑ j ∈ Finset.univ.erase p, a j • chartVector k p j :=
  chart_reconstruct k a p hp _ ha
    (by intro r hr h; exact False.elim (h (Finset.mem_erase.mpr ⟨hr, Finset.mem_univ _⟩)))

theorem ordinary_chart_reconstruct (e : Fin D) (k a : Vec D) (p : Fin D) (hp : k p ≠ 0)
    (ha : a ∈ ordinarySpace e k) :
    a = ∑ j ∈ (Finset.univ.erase e).erase p, a j • chartVector k p j := by
  have ha' := (mem_ordinarySpace e k a).mp ha
  apply chart_reconstruct k a p hp _ ha'.2
  intro r hr h
  have hre : r = e := by
    by_contra hn
    exact h (Finset.mem_erase.mpr ⟨hr, Finset.mem_erase.mpr ⟨hn, Finset.mem_univ _⟩⟩)
  exact hre ▸ ha'.1

theorem transverse_pivot_exists (k : Vec D) (hk : k ≠ 0) : ∃ p, k p ≠ 0 := by
  by_contra h
  push Not at h
  exact hk (funext h)

theorem chartVector_ordinary (e : Fin D) (k : Vec D) (p j : Fin D)
    (hp : k p ≠ 0) (hpe : p ≠ e) (hje : j ≠ e) :
    chartVector k p j ∈ ordinarySpace e k := by
  rw [mem_ordinarySpace]
  exact ⟨by simp only [chartVector_coordinate k p j e hpe.symm, hje.symm, if_false],
    chartVector_transverse k p j hp⟩

theorem ordinary_pivot_exists (e : Fin D) (k : Vec D) (hq : perpSq e k ≠ 0) :
    ∃ p, p ≠ e ∧ k p ≠ 0 := by
  by_contra h
  push Not at h
  apply hq
  apply (perpSq_zero_iff e k).mpr
  funext r
  by_cases hr : r = e
  · subst r; simp
  · simp [unit, hr, h r hr]

def chartVectorTree (p j r : Fin (n + 1)) : Tree (n + 1) :=
  .sub (.scalar (if r = j then 1 else 0))
    (.mul (.div (.atom (.wavevector j)) (.atom (.wavevector p)))
      (.scalar (if r = p then 1 else 0)))

theorem chartVectorTree_eval (rho mu sigma scale omega : ℝ) (u k : Vec (n + 1))
    (J : Jet (n + 1)) (p j r : Fin (n + 1)) :
    Expr.eval (atomValues rho mu sigma scale omega u k J) (chartVectorTree p j r) =
      chartVector k p j r := by
  simp only [chartVectorTree, Expr.eval, atomValues, chartVector, Pi.sub_apply,
    Pi.smul_apply, smul_eq_mul, unit, Pi.single_apply]

theorem chartVectorTree_hasDim (u : UnitAssignment) (p j r : Fin (n + 1)) :
    Expr.HasDim (atomUnits u) (chartVectorTree p j r) 1 := by
  apply Expr.HasDim.sub (.scalar _)
  simpa only [atomUnits, div_self', one_mul] using
    (Expr.HasDim.mul
      (Expr.HasDim.div (Expr.HasDim.atom (u := atomUnits u) (.wavevector j))
        (Expr.HasDim.atom (u := atomUnits u) (.wavevector p)))
      (Expr.HasDim.scalar (if r = p then 1 else 0)))

end
end S10Audit
