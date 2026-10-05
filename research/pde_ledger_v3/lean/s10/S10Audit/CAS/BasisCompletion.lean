import S10Audit.CAS.Support

set_option backward.isDefEq.respectTransparency false

/-! Completeness of the displayed D3 chart. Exceptional directions are handled
by the earlier spectral theorems, and are not assertions of this generic chart. -/

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section

theorem referenceBasis_independent (sigma : ℝ) (k : Vec 3) (r : Fin 3)
    (h : BasisDomain sigma k r) :
    LinearIndependent ℝ (fun _ : Fin 1 => referenceBasis sigma k r) := by
  rw [linearIndependent_unique_iff]
  intro heq
  have hv := congrFun heq 2
  rw [referenceBasis_last sigma k r h] at hv
  norm_num at hv

theorem staticBasis_span (sigma : ℝ) (k : Vec 3) (h : k 2 ≠ 0) :
    Submodule.span ℝ {referenceBasis sigma k 0} = longitudinalSpace k :=
  normalized_span k (k 2) h

theorem ordinaryBasis_span (sigma : ℝ) (k : Vec 3) (h : k 1 ≠ 0) :
    Submodule.span ℝ {referenceBasis sigma k 1} = ordinarySpace 0 k := by
  apply le_antisymm
  · apply Submodule.span_le.mpr
    intro a ha
    rcases Set.mem_singleton_iff.mp ha with rfl
    exact chartVector_ordinary 0 k 1 2 h (by decide) (by decide)
  · intro a ha
    have heq := ordinary_chart_reconstruct 0 k a 1 h ha
    have hf : ((Finset.univ.erase (0 : Fin 3)).erase 1) = {2} := by decide
    rw [hf, Finset.sum_singleton] at heq
    apply Submodule.mem_span_singleton.mpr
    exact ⟨a 2, heq.symm⟩

theorem extraBasis_span (sigma : ℝ) (k : Vec 3)
    (hs : sigma ≠ 0) (h0 : k 0 ≠ 0) (h2 : k 2 ≠ 0) :
    Submodule.span ℝ {referenceBasis sigma k 2} = longitudinalSpace (extraVector 0 sigma k) :=
  normalized_span _ _ (neg_ne_zero.mpr (mul_ne_zero (mul_ne_zero hs h0) h2))

def GenericChart (sigma : ℝ) (k : Vec 3) : Prop :=
  0 < sigma ∧ sigma ≠ 1 ∧ k 0 ≠ 0 ∧ k 1 ≠ 0 ∧ k 2 ≠ 0

theorem genericChart_basisDomain (sigma : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : BasisDomain sigma k r := by
  rcases h with ⟨hs, _, h0, h1, h2⟩
  fin_cases r
  · exact h2
  · exact h1
  · exact ⟨ne_of_gt hs, h0, h2⟩

theorem referenceRoot_normalized (rho mu sigma : ℝ) (k : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) :
    rho * referenceRoot rho mu sigma k r / mu =
      (if r = 0 then 0 else if r = 1 then normSq k else extraValue 0 sigma k) := by
  fin_cases r <;> simp [referenceRoot, coneValue, extraConeValue] <;> field_simp

theorem referenceMatrix_basis_complete (rho mu sigma : ℝ) (k a : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    (referenceMatrix rho mu sigma (referenceRoot rho mu sigma k r) k).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ {referenceBasis sigma k r} := by
  rcases h with ⟨hs, hs1, h0, h1, h2⟩
  have hk : k ≠ 0 := by intro heq; exact h0 (congrFun heq 0)
  have hq : perpSq 0 k ≠ 0 := by
    have hpos := sq_pos_of_ne_zero h1
    simp only [perpSq, normSq, dot, Fin.sum_univ_three]
    nlinarith [sq_nonneg (k 2)]
  rw [referenceMatrix_normalized rho mu sigma _ k a hm, smul_eq_zero,
    or_iff_right (show ¬ mu = 0 from hm), referenceRoot_normalized rho mu sigma k r hr hm]
  change a ∈ modeSpace 0 sigma _ k ↔ _
  fin_cases r
  · change a ∈ modeSpace 0 sigma 0 k ↔ a ∈ Submodule.span ℝ {referenceBasis sigma k 0}
    rw [zero_modeSpace hk, staticBasis_span sigma k h2]
  · change a ∈ modeSpace 0 sigma (normSq k) k ↔ a ∈ Submodule.span ℝ {referenceBasis sigma k 1}
    rw [ordinary_modeSpace hs1 hq hk, ordinaryBasis_span sigma k h1]
  · change a ∈ modeSpace 0 sigma (extraValue 0 sigma k) k ↔
      a ∈ Submodule.span ℝ {referenceBasis sigma k 2}
    rw [extra_modeSpace hs hs1 hq hk, extraBasis_span sigma k (ne_of_gt hs) h0 h2]

end
end S10Audit.CAS
