import S10Pilot.VariationalCertificate
import S10Pilot.PhaseAverage
import S10Controls.Spectrum

/-! The supplied coefficient-scale and sign-flip controls enter at the density.
An arbitrary real squared-frequency variable retains the negative root. -/

namespace S10ScalarControls
open S10Pilot
open MeasureTheory
noncomputable section
variable {D : ℕ}

def lagrangian (c rho mu : ℝ) (J : Jet D) : ℝ :=
  rho / 2 * normSq (J 0) - c * mu / 2 * stiffness J

theorem lagrangian_eq (c rho mu : ℝ) (J : Jet D) :
    lagrangian c rho mu J = S10Pilot.lagrangian rho (c * mu) J := rfl

theorem lagrangian_variation (c rho mu : ℝ) (J H : Jet D) :
    HasDerivAt (fun s : ℝ => lagrangian c rho mu (J + s • H))
      (S10Pilot.variationDensity rho (c * mu) J H) 0 :=
  S10Pilot.lagrangian_variation rho (c * mu) J H

def relativeAction (c rho mu : ℝ) (u h : Point D → Vec D) (s : ℝ) : ℝ :=
  ∫ x, lagrangian c rho mu (fieldJet (u + s • h) x) - lagrangian c rho mu (fieldJet u x)

theorem relativeAction_eq (c rho mu : ℝ) (u h : Point D → Vec D) :
    relativeAction c rho mu u h = S10Pilot.relativeAction rho (c * mu) u h := rfl

theorem relativeAction_hasDerivAt (c rho mu : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) :
    HasDerivAt (relativeAction c rho mu u h)
      (∫ x, S10Pilot.linearDensity rho (c * mu) (fieldJet u x) (fieldJet h x)) 0 :=
  S10Pilot.relativeAction_hasDerivAt hu hh rho (c * mu)

def ActionStationary (c rho mu : ℝ) (u : Point D → Vec D) : Prop :=
  ∀ h, TestField h → deriv (relativeAction c rho mu u h) 0 = 0

theorem actionStationary_eq (c rho mu : ℝ) (u : Point D → Vec D) :
    ActionStationary c rho mu u ↔ S10Pilot.ActionStationary rho (c * mu) u := Iff.rfl

def spectralAction (c rho mu z : ℝ) (k a : Vec D) : ℝ :=
  rho / 2 * z * normSq a - c * mu / 2 * (normSq k * normSq a - dot k a ^ 2)

theorem modalAction_eq (c rho mu omega : ℝ) (k a : Vec D) :
    lagrangian c rho mu (modeJet omega k a) = spectralAction c rho mu (omega ^ 2) k a :=
  S10Pilot.modalAction_eq rho (c * mu) omega k a

theorem phaseAverage_eq (c rho mu omega : ℝ) (k a : Vec D) :
    (1 / (2 * Real.pi)) * (∫ phi in (0 : ℝ)..(2 * Real.pi),
      lagrangian c rho mu ((-Real.sin phi) • modeJet (-omega) k a)) =
      (1 / 2 : ℝ) * spectralAction c rho mu (omega ^ 2) k a :=
  (S10Pilot.phaseAverage_eq rho (c * mu) omega k a).trans
    (congrArg ((1 / 2 : ℝ) * ·) (modalAction_eq c rho mu omega k a))

def spectralOperator (c rho mu z : ℝ) (k a : Vec D) : Vec D :=
  (rho * z - c * mu * normSq k) • a + (c * mu * dot k a) • k

theorem spectralAction_variation (c rho mu z : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => spectralAction c rho mu z k (a + s • b))
      (dot (spectralOperator c rho mu z k a) b) 0 := by
  have heq : ∀ s : ℝ, spectralAction c rho mu z k (a + s • b) =
      spectralAction c rho mu z k a + s * dot (spectralOperator c rho mu z k a) b +
      s ^ 2 * spectralAction c rho mu z k b := by
    intro s
    unfold spectralAction
    simp only [normSq, dot_add_left, dot_add_right, dot_smul_left, dot_smul_right,
      spectralOperator]
    rw [dot_comm b a]
    ring
  have h := ((hasDerivAt_const (0 : ℝ) (spectralAction c rho mu z k a)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (dot (spectralOperator c rho mu z k a) b))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (spectralAction c rho mu z k b))
  convert! h using 1
  · exact funext heq
  · simp

theorem modalOperator_eq (c rho mu omega : ℝ) (k a : Vec D) :
    S10Pilot.modalOperator rho (c * mu) omega k a = spectralOperator c rho mu (omega ^ 2) k a := by
  rw [S10Pilot.modalOperator_eq]
  rfl

theorem variational_planeWave_iff (c rho mu omega : ℝ) (k a : Vec D) :
    ActionStationary c rho mu (planeWave omega k a) ↔
      spectralOperator c rho mu (omega ^ 2) k a = 0 := by
  rw [actionStationary_eq, S10Pilot.actionStationary_planeWave_iff,
    S10Pilot.modal_stationary_iff, modalOperator_eq]

def coneValue (c rho mu : ℝ) (k : Vec D) : ℝ := (c * mu / rho) * normSq k

/-- Exhausts all nonzero real squared-frequency roots, including negative ones. -/
theorem nonzero_root_iff {c rho mu z : ℝ} {k a : Vec D}
    (hrho : rho ≠ 0) (hz : z ≠ 0) (ha : a ≠ 0) :
    spectralOperator c rho mu z k a = 0 ↔ dot k a = 0 ∧ z = coneValue c rho mu k := by
  constructor
  · intro h
    have hd := congrArg (dot k) h
    simp only [spectralOperator, dot_add_right, dot_smul_right, dot_zero] at hd
    have he : rho * z * dot k a = 0 := by
      change _ + _ * normSq k = 0 at hd
      nlinarith
    have ht := (mul_eq_zero.mp he).resolve_left (mul_ne_zero hrho hz)
    have hm : (rho * z - c * mu * normSq k) • a = 0 := by simpa [spectralOperator, ht] using h
    have hf := (smul_eq_zero.mp hm).resolve_right ha
    refine ⟨ht, ?_⟩
    unfold coneValue
    field_simp
    nlinarith
  · rintro ⟨ht, hf⟩
    have he : rho * z - c * mu * normSq k = 0 := by
      rw [hf]
      unfold coneValue
      field_simp
      ring
    simp [spectralOperator, ht, he]

theorem coneValue_pos {c rho mu : ℝ} {k : Vec D}
    (hc : 0 < c) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    0 < coneValue c rho mu k :=
  mul_pos (div_pos (mul_pos hc hmu) hrho) (dot_self_pos hk)

theorem coneValue_neg {c rho mu : ℝ} {k : Vec D}
    (hc : c < 0) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    coneValue c rho mu k < 0 :=
  mul_neg_of_neg_of_pos (div_neg_of_neg_of_pos (mul_neg_of_neg_of_pos hc hmu) hrho) (dot_self_pos hk)

def spectralMap (c rho mu z : ℝ) (k : Vec D) : Vec D →ₗ[ℝ] Vec D where
  toFun := spectralOperator c rho mu z k
  map_add' a b := by ext i; simp [spectralOperator, dot_add_right]; ring
  map_smul' s a := by ext i; simp [spectralOperator, dot_smul_right]; ring

def modeSpace (c rho mu z : ℝ) (k : Vec D) : Submodule ℝ (Vec D) :=
  (spectralMap c rho mu z k).ker

theorem cone_modeSpace {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 0) (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) :
    modeSpace c rho mu (coneValue c rho mu k) k = transverseSpace k := by
  have hz : coneValue c rho mu k ≠ 0 :=
    mul_ne_zero (div_ne_zero (mul_ne_zero hc hmu) hrho) (ne_of_gt (dot_self_pos hk))
  ext a
  change spectralOperator c rho mu (coneValue c rho mu k) k a = 0 ↔ dot k a = 0
  by_cases ha : a = 0
  · simp [ha, spectralOperator]
  · rw [nonzero_root_iff hrho hz ha]
    simp

theorem zero_modeSpace {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) :
    modeSpace c rho mu 0 k = longitudinalSpace k := by
  ext a
  change spectralOperator c rho mu 0 k a = 0 ↔ _
  have h := S10Pilot.zero_frequency_iff (rho := rho) (a := a) (mul_ne_zero hc hmu) hk
  rw [S10Pilot.modal_stationary_iff, modalOperator_eq] at h
  simpa using h

theorem cone_counts {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 0) (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) :
    Module.finrank ℝ (modeSpace c rho mu (coneValue c rho mu k) k) = D - 1 ∧
    Module.finrank ℝ (modeSpace c rho mu (coneValue c rho mu k) k ⊓ transverseSpace k :
      Submodule ℝ (Vec D)) = D - 1 := by
  rw [cone_modeSpace hc hrho hmu hk, inf_idem, transverseSpace_finrank hk]
  exact ⟨rfl, rfl⟩

theorem positive_variational_certificate {c rho mu : ℝ} {k : Vec D}
    (hc : 0 < c) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ omega : ℝ, 0 < omega ∧ omega ^ 2 = coneValue c rho mu k ∧
      (∀ a, ActionStationary c rho mu (planeWave omega k a) ↔ a ∈ transverseSpace k) ∧
      Module.finrank ℝ (transverseSpace k) = D - 1 := by
  obtain ⟨w, hw, hf, ht, hd, _, _⟩ := s10_variational_certificate hrho (mul_pos hc hmu) hk
  exact ⟨w, hw, hf, ht, hd⟩

theorem negative_control_no_real_wave {c rho mu omega : ℝ} {k a : Vec D}
    (hc : c < 0) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0)
    (hw : omega ≠ 0) (ha : a ≠ 0) : ¬ ActionStationary c rho mu (planeWave omega k a) := by
  intro h
  rw [variational_planeWave_iff, nonzero_root_iff (ne_of_gt hrho) (pow_ne_zero 2 hw) ha] at h
  have hn := coneValue_neg hc hrho hmu hk
  nlinarith [sq_nonneg omega]

theorem coneValue_scaling (c rho mu lambda : ℝ) (k : Vec D) :
    coneValue c rho mu (lambda • k) = lambda ^ 2 * coneValue c rho mu k := by
  simp only [coneValue, normSq_smul]
  ring

theorem coneValue_scaling_ratio {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 0) (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) (lambda : ℝ) :
    coneValue c rho mu (lambda • k) / coneValue c rho mu k = lambda ^ 2 := by
  rw [coneValue_scaling, mul_div_cancel_right₀]
  exact mul_ne_zero (div_ne_zero (mul_ne_zero hc hmu) hrho) (ne_of_gt (dot_self_pos hk))

theorem signflip_counts {rho mu : ℝ} {k : Vec D}
    (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    coneValue (-1) rho mu k < 0 ∧
    Module.finrank ℝ (modeSpace (-1) rho mu (coneValue (-1) rho mu k) k) = D - 1 ∧
    Module.finrank ℝ (modeSpace (-1) rho mu (coneValue (-1) rho mu k) k ⊓ transverseSpace k :
      Submodule ℝ (Vec D)) = D - 1 :=
  ⟨coneValue_neg (by norm_num) hrho hmu hk,
    cone_counts (by norm_num) (ne_of_gt hrho) (ne_of_gt hmu) hk⟩

theorem zero_counts {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) :
    Module.finrank ℝ (modeSpace c rho mu 0 k) = 1 ∧
    Module.finrank ℝ (modeSpace c rho mu 0 k ⊓ transverseSpace k : Submodule ℝ (Vec D)) = 0 := by
  rw [zero_modeSpace hc hmu hk, longitudinalSpace_finrank hk,
    S10Controls.longitudinal_inf_transverse hk]
  simp

theorem negative_root_exists {c rho mu : ℝ} {k : Vec D}
    (hD : 2 ≤ D) (hc : c < 0) (hrho : 0 < rho) (hmu : 0 < mu) (hk : k ≠ 0) :
    ∃ z : ℝ, ∃ a : Vec D, z < 0 ∧ a ≠ 0 ∧ spectralOperator c rho mu z k a = 0 := by
  have hp : 0 < Module.finrank ℝ (transverseSpace k) := by
    rw [transverseSpace_finrank hk]
    omega
  obtain ⟨a, ha⟩ := Module.finrank_pos_iff_exists_ne_zero.mp hp
  refine ⟨coneValue c rho mu k, a, coneValue_neg hc hrho hmu hk,
    (fun h => ha (Subtype.ext h)), ?_⟩
  have hm : (a : Vec D) ∈ modeSpace c rho mu (coneValue c rho mu k) k := by
    rw [cone_modeSpace (ne_of_lt hc) (ne_of_gt hrho) (ne_of_gt hmu) hk]
    exact a.property
  exact hm

theorem coefficient_changes_frequency {c rho mu : ℝ} {k : Vec D}
    (hc : c ≠ 1) (hrho : rho ≠ 0) (hmu : mu ≠ 0) (hk : k ≠ 0) :
    coneValue c rho mu k ≠ S10Pilot.coneValue rho mu k := by
  have he : coneValue c rho mu k = c * S10Pilot.coneValue rho mu k := by
    unfold coneValue S10Pilot.coneValue
    ring
  have hn : S10Pilot.coneValue rho mu k ≠ 0 :=
    mul_ne_zero (div_ne_zero hmu hrho) (ne_of_gt (dot_self_pos hk))
  rw [he]
  intro h
  have heq : c * S10Pilot.coneValue rho mu k = 1 * S10Pilot.coneValue rho mu k := by simpa using h
  exact hc (mul_right_cancel₀ hn heq)

end
end S10ScalarControls
