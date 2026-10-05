import S11Homogeneous.Spectrum

/-! Phase matching only. The sign of the normal wavevector square does not
establish interface overlap, bound-state existence or a radiation rate. -/

namespace S11Homogeneous
open S10Pilot
noncomputable section
variable {D : ℕ}

def normalWaveSq (rho B cs : ℝ) (k : Vec D) : ℝ :=
  normSq k * (B / (rho * cs ^ 2) - 1)

theorem phase_matching {rho B cs omega : ℝ} (hrho : rho ≠ 0) (hcs : cs ≠ 0)
    (k : Vec D) (hf : omega ^ 2 = coneValue rho B k) :
    normalWaveSq rho B cs k = omega ^ 2 / cs ^ 2 - normSq k := by
  rw [hf]
  unfold normalWaveSq coneValue
  field_simp

theorem threshold_classification {rho B cs : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hcs : 0 < cs) (hk : k ≠ 0) :
    (normalWaveSq rho B cs k < 0 ↔ B < rho * cs ^ 2) ∧
    (normalWaveSq rho B cs k = 0 ↔ B = rho * cs ^ 2) ∧
    (0 < normalWaveSq rho B cs k ↔ rho * cs ^ 2 < B) := by
  have hn := dot_self_pos hk
  have hp : 0 < rho * cs ^ 2 := mul_pos hrho (sq_pos_of_pos hcs)
  have he : normalWaveSq rho B cs k = (normSq k / (rho * cs ^ 2)) * (B - rho * cs ^ 2) := by
    unfold normalWaveSq
    field_simp
  have hq : 0 < normSq k / (rho * cs ^ 2) := div_pos hn hp
  rw [he]
  constructor
  · constructor <;> intro h <;> nlinarith
  · constructor
    · constructor <;> intro h <;> nlinarith
    · constructor <;> intro h <;> nlinarith

end
end S11Homogeneous
