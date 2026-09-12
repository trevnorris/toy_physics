import S10Audit.CAS.CountSupport
import S10Audit.CAS.Bindings

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section

/-- Counts on the explicit chart used by the printed generic D3 bases.
The extra root has no transverse kernel on this chart. -/
theorem genericChart_mode_dimensions (sigma : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) :
    Module.finrank ℝ (rootMode sigma k r) = 1 ∧
      Module.finrank ℝ (rootMode sigma k r ⊓ transverseSpace k : Submodule ℝ (Vec 3)) =
        (if r = 1 then 1 else 0) := by
  rcases h with ⟨hs, hs1, h0, h1, _h2⟩
  have hk : k ≠ 0 := by intro heq; exact h0 (congrFun heq 0)
  have hq : perpSq 0 k ≠ 0 := by
    have hpos := sq_pos_of_ne_zero h1
    simp only [perpSq, normSq, dot, Fin.sum_univ_three]
    nlinarith [sq_nonneg (k 2)]
  have ho := oblique_counts hs hs1 hk hq h0
  fin_cases r
  · exact zero_counts hk
  · exact ⟨ho.1, ho.2.1⟩
  · exact ⟨ho.2.2.1, ho.2.2.2⟩

end
end S10Audit.CAS
