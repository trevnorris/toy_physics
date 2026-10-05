import S10Controls.FiniteAction
import S10Controls.Spectrum
import S10Pilot.VariationalCertificate

/-! Integrated variational certificates and basis-independent counts for both controls. -/

namespace S10Controls
open S10Pilot
noncomputable section
variable {D : ℕ}

theorem actionStationary_iff_mem_amplitudeSpace (form : Form) (rho mu omega : ℝ) (k a : Vec D) :
    ActionStationary form rho mu (planeWave omega k a) ↔
      a ∈ amplitudeSpace form rho mu omega k := by
  rw [actionStationary_planeWave_iff, mem_amplitudeSpace]

theorem full_variational_iff {rho mu omega : ℝ} {k a : Vec D} (hrho : rho ≠ 0) :
    ActionStationary .fullGradient rho mu (planeWave omega k a) ↔
      a = 0 ∨ omega ^ 2 = coneValue rho mu k := by
  rw [actionStationary_planeWave_iff, full_stationary_iff hrho]

theorem div_propagating_variational_iff {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0) (ha : a ≠ 0) :
    ActionStationary .divergenceOnly rho mu (planeWave omega k a) ↔
      a ∈ longitudinalSpace k ∧ omega ^ 2 = coneValue rho mu k := by
  rw [actionStationary_planeWave_iff, div_propagating_iff hrho hw hk ha]

theorem full_mode_counts {rho mu omega : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0)
    (hf : omega ^ 2 = coneValue rho mu k) :
    Module.finrank ℝ (amplitudeSpace .fullGradient rho mu omega k) = D ∧
    Module.finrank ℝ (amplitudeSpace .fullGradient rho mu omega k ⊓ transverseSpace k :
      Submodule ℝ (Vec D)) = D - 1 ∧
    Module.finrank ℝ (amplitudeSpace .fullGradient rho mu 0 k) = 0 := by
  rw [full_on_cone_space hrho hf, full_zero_space hmu hk, top_inf_eq, transverseSpace_finrank hk]
  simp

theorem div_mode_counts {rho mu omega : ℝ} {k : Vec D}
    (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0)
    (hf : omega ^ 2 = coneValue rho mu k) :
    Module.finrank ℝ (amplitudeSpace .divergenceOnly rho mu omega k) = 1 ∧
    Module.finrank ℝ (amplitudeSpace .divergenceOnly rho mu omega k ⊓ transverseSpace k :
      Submodule ℝ (Vec D)) = 0 ∧
    Module.finrank ℝ (amplitudeSpace .divergenceOnly rho mu 0 k) = D - 1 := by
  rw [div_on_cone_space hrho hw hk hf, div_zero_space hmu hk]
  simp [longitudinal_inf_transverse hk, longitudinalSpace_finrank hk, transverseSpace_finrank hk]

/-- At one common positive cone frequency, all three actions are compared using
stationarity against every smooth compact vector variation. Their zero-frequency
amplitude spaces are included, so a surviving static sector is not discarded. -/
theorem stiffness_comparison_certificate {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue rho mu k ∧
      (∀ a : Vec D, S10Pilot.ActionStationary rho mu (planeWave omega k a) ↔
        a ∈ transverseSpace k) ∧
      (∀ a : Vec D, ActionStationary .fullGradient rho mu (planeWave omega k a)) ∧
      (∀ a : Vec D, ActionStationary .divergenceOnly rho mu (planeWave omega k a) ↔
        a ∈ longitudinalSpace k) ∧
      (∀ a : Vec D, S10Pilot.ActionStationary rho mu (planeWave 0 k a) ↔
        a ∈ longitudinalSpace k) ∧
      (∀ a : Vec D, ActionStationary .fullGradient rho mu (planeWave 0 k a) ↔ a = 0) ∧
      (∀ a : Vec D, ActionStationary .divergenceOnly rho mu (planeWave 0 k a) ↔
        a ∈ transverseSpace k) := by
  obtain ⟨w, hw, hf, ht, _, hl, _⟩ := s10_variational_certificate hrho hmu hk
  refine ⟨w, hw, hf, ht, ?_, ?_, hl, ?_, ?_⟩
  · intro a
    exact (full_variational_iff (ne_of_gt hrho)).mpr (Or.inr hf)
  · intro a
    rw [actionStationary_planeWave_iff,
      div_on_cone_iff (ne_of_gt hrho) (ne_of_gt hw) hk hf]
  · intro a
    rw [actionStationary_iff_mem_amplitudeSpace, full_zero_space (ne_of_gt hmu) hk]
    rfl
  · intro a
    rw [actionStationary_planeWave_iff, div_zero_frequency_iff (ne_of_gt hmu) hk]

/-- A nonzero longitudinal amplitude propagates under both controls but fails
stationarity under the curl-only baseline, for every nonzero wavevector. -/
theorem longitudinal_control_witness {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧
      ActionStationary .fullGradient rho mu (planeWave omega k k) ∧
      ActionStationary .divergenceOnly rho mu (planeWave omega k k) ∧
      ¬ S10Pilot.ActionStationary rho mu (planeWave omega k k) := by
  obtain ⟨w, hw, _, ht, hfull, hdiv, _, _, _⟩ := stiffness_comparison_certificate hrho hmu hk
  refine ⟨w, hw, hfull k, (hdiv k).mpr (Submodule.mem_span_singleton.mpr ⟨1, by simp⟩), ?_⟩
  intro h
  have hz : dot k k = 0 := (ht k).mp h
  exact (ne_of_gt (dot_self_pos hk)) hz

/-- A nonzero transverse amplitude on the cone passes FULLGRAD, while DIVONLY
has a compact test variation whose first variation is nonzero. -/
theorem transverse_control_witness {rho mu omega : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hw : omega ≠ 0) (hk : k ≠ 0) (ha : a ≠ 0)
    (ht : a ∈ transverseSpace k) (hf : omega ^ 2 = coneValue rho mu k) :
    ActionStationary .fullGradient rho mu (planeWave omega k a) ∧
    ∃ h : Point D → Vec D, TestField h ∧
      deriv (relativeAction .divergenceOnly rho mu (planeWave omega k a) h) 0 ≠ 0 := by
  refine ⟨(full_variational_iff hrho).mpr (Or.inr hf), ?_⟩
  have hn : ¬ ActionStationary .divergenceOnly rho mu (planeWave omega k a) := by
    intro h
    have hl := div_stationary_longitudinal hrho hw
      ((actionStationary_planeWave_iff _ _ _ _ _ _).mp h)
    have hb : a ∈ longitudinalSpace k ⊓ transverseSpace k := ⟨hl, ht⟩
    rw [longitudinal_inf_transverse hk] at hb
    exact ha hb
  simpa [ActionStationary] using hn

end
end S10Controls
