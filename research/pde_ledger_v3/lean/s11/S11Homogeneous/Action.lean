import S10Pilot.Variation
import S10Controls.Variation
import S10Pilot.PhaseAverage
import S10Controls.PhaseAverage

/-! The supplied homogeneous action, combined from the existing curl and
divergence densities. The compression summand has zero inertia: kinetic energy
is counted once. Integrated variation is reused, not assumed. -/

namespace S11Homogeneous
open S10Pilot MeasureTheory
noncomputable section
variable {D : ℕ}

def lagrangian (rho mu B : ℝ) (J : Jet D) : ℝ :=
  rho / 2 * normSq (J 0) - mu / 2 * stiffness J -
    B / 2 * S10Controls.divergenceOnlyStiffness J

theorem lagrangian_split (rho mu B : ℝ) (J : Jet D) :
    lagrangian rho mu B J = S10Pilot.lagrangian rho mu J +
      S10Controls.lagrangian .divergenceOnly 0 B J := by
  simp [lagrangian, S10Pilot.lagrangian, S10Controls.lagrangian, S10Controls.stiffness]
  ring

theorem lagrangian_variation (rho mu B : ℝ) (J H : Jet D) :
    HasDerivAt (fun s : ℝ => lagrangian rho mu B (J + s • H))
      (S10Pilot.variationDensity rho mu J H +
        S10Controls.variationDensity .divergenceOnly 0 B J H) 0 := by
  simp_rw [lagrangian_split]
  exact (S10Pilot.lagrangian_variation rho mu J H).add
    (S10Controls.lagrangian_variation .divergenceOnly 0 B J H)

def modalOperator (rho mu B omega : ℝ) (k a : Vec D) : Vec D :=
  (rho * omega ^ 2 - mu * normSq k) • a + ((mu - B) * dot k a) • k

theorem modalOperator_split (rho mu B omega : ℝ) (k a : Vec D) :
    modalOperator rho mu B omega k a = S10Pilot.modalOperator rho mu omega k a +
      S10Controls.modalOperator .divergenceOnly 0 B omega k a := by
  ext i
  simp [modalOperator, S10Pilot.modalOperator, S10Controls.modalOperator]
  ring

def modalAction (rho mu B omega : ℝ) (k a : Vec D) : ℝ :=
  lagrangian rho mu B (modeJet omega k a)

theorem modalAction_variation (rho mu B omega : ℝ) (k a b : Vec D) :
    HasDerivAt (fun s : ℝ => modalAction rho mu B omega k (a + s • b))
      (dot (modalOperator rho mu B omega k a) b) 0 := by
  simp only [modalAction, lagrangian_split, modalOperator_split, dot_add_left]
  exact (S10Pilot.modalAction_variation rho mu omega k a b).add
    (S10Controls.modalAction_variation .divergenceOnly 0 B omega k a b)

theorem phaseDensity_eq (rho mu B omega : ℝ) (k a : Vec D) (phi : ℝ) :
    lagrangian rho mu B ((-Real.sin phi) • modeJet (-omega) k a) =
      Real.sin phi ^ 2 * modalAction rho mu B omega k a := by
  rw [lagrangian_split]
  change S10Pilot.phaseDensity rho mu omega k a phi +
    S10Controls.phaseDensity .divergenceOnly 0 B omega k a phi = _
  rw [S10Pilot.phaseDensity_eq, S10Controls.phaseDensity_eq]
  simp only [modalAction, lagrangian_split, S10Pilot.modalAction, S10Controls.modalAction]
  ring

def phaseAverage (rho mu B omega : ℝ) (k a : Vec D) : ℝ :=
  (1 / (2 * Real.pi)) * ∫ phi in (0 : ℝ)..(2 * Real.pi),
    lagrangian rho mu B ((-Real.sin phi) • modeJet (-omega) k a)

theorem phaseAverage_eq (rho mu B omega : ℝ) (k a : Vec D) :
    phaseAverage rho mu B omega k a = (1 / 2 : ℝ) * modalAction rho mu B omega k a := by
  unfold phaseAverage
  simp_rw [phaseDensity_eq]
  rw [intervalIntegral.integral_mul_const, integral_sin_sq]
  simp only [Real.sin_zero, zero_mul, Real.sin_two_pi, sub_zero, zero_add]
  field_simp

def ModalStationary (rho mu B omega : ℝ) (k a : Vec D) : Prop :=
  ∀ b : Vec D, deriv (fun s : ℝ => modalAction rho mu B omega k (a + s • b)) 0 = 0

theorem modal_stationary_iff (rho mu B omega : ℝ) (k a : Vec D) :
    ModalStationary rho mu B omega k a ↔ modalOperator rho mu B omega k a = 0 := by
  constructor
  · intro h
    have hs := h (modalOperator rho mu B omega k a)
    rw [(modalAction_variation _ _ _ _ _ _ _).deriv] at hs
    by_contra hn
    exact (ne_of_gt (dot_self_pos hn)) hs
  · intro h b
    rw [(modalAction_variation _ _ _ _ _ _ _).deriv, h, zero_dot]

def relativeAction (rho mu B : ℝ) (u h : Point D → Vec D) (s : ℝ) : ℝ :=
  ∫ x, lagrangian rho mu B (fieldJet (u + s • h) x) - lagrangian rho mu B (fieldJet u x)

theorem relativeAction_split {u h : Point D → Vec D} (hu : SmoothField u)
    (hh : TestField h) (rho mu B s : ℝ) :
    relativeAction rho mu B u h s = S10Pilot.relativeAction rho mu u h s +
      S10Controls.relativeAction .divergenceOnly 0 B u h s := by
  unfold relativeAction S10Pilot.relativeAction S10Controls.relativeAction
  simp_rw [lagrangian_split]
  have rearrange : ∀ x : Point D,
      (S10Pilot.lagrangian rho mu (fieldJet (u + s • h) x) +
        S10Controls.lagrangian .divergenceOnly 0 B (fieldJet (u + s • h) x)) -
      (S10Pilot.lagrangian rho mu (fieldJet u x) +
        S10Controls.lagrangian .divergenceOnly 0 B (fieldJet u x)) =
      (S10Pilot.lagrangian rho mu (fieldJet (u + s • h) x) -
        S10Pilot.lagrangian rho mu (fieldJet u x)) +
      (S10Controls.lagrangian .divergenceOnly 0 B (fieldJet (u + s • h) x) -
        S10Controls.lagrangian .divergenceOnly 0 B (fieldJet u x)) := by
    intro x
    ring
  simp_rw [rearrange]
  exact integral_add (S10Pilot.relative_density_integrable hu hh rho mu s)
    (S10Controls.relative_density_integrable .divergenceOnly hu hh 0 B s)

/-- Sum of two already action-derived local Euler–Lagrange expressions. -/
def eulerLagrange (rho mu B : ℝ) (u : Point D → Vec D) (x : Point D) : Vec D :=
  S10Pilot.eulerLagrange rho mu u x + S10Controls.eulerLagrange .divergenceOnly 0 B u x

theorem relativeAction_deriv_eq_eulerLagrange {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu B : ℝ) :
    deriv (relativeAction rho mu B u h) 0 =
      ∫ x, ∑ i : Fin D, eulerLagrange rho mu B u x i * h x i := by
  have heq : relativeAction rho mu B u h =
      S10Pilot.relativeAction rho mu u h +
        S10Controls.relativeAction .divergenceOnly 0 B u h :=
    funext (relativeAction_split hu hh rho mu B)
  rw [heq, ((S10Pilot.relativeAction_hasDerivAt hu hh rho mu).add
    (S10Controls.relativeAction_hasDerivAt .divergenceOnly hu hh 0 B)).deriv]
  rw [S10Pilot.integrated_variation_by_parts hu hh,
    S10Controls.integrated_variation_by_parts .divergenceOnly hu hh]
  have hi₁ : Integrable (fun x => ∑ i : Fin D, S10Pilot.eulerLagrange rho mu u x i * h x i) :=
    integrable_finsetSum _ (fun i _ => integrable_mul_compact
      (S10Pilot.smooth_eulerLagrange hu rho mu i).continuous (hh.1 i).continuous (hh.2 i))
  have hi₂ : Integrable (fun x => ∑ i : Fin D,
      S10Controls.eulerLagrange .divergenceOnly 0 B u x i * h x i) :=
    integrable_finsetSum _ (fun i _ => integrable_mul_compact
      (S10Controls.smooth_eulerLagrange .divergenceOnly hu 0 B i).continuous
      (hh.1 i).continuous (hh.2 i))
  rw [← integral_add hi₁ hi₂]
  congr 1
  funext x
  simp [eulerLagrange, add_mul, Finset.sum_add_distrib]

def ActionStationary (rho mu B : ℝ) (u : Point D → Vec D) : Prop :=
  ∀ h : Point D → Vec D, TestField h → deriv (relativeAction rho mu B u h) 0 = 0

theorem actionStationary_iff_eulerLagrange {u : Point D → Vec D}
    (hu : SmoothField u) (rho mu B : ℝ) :
    ActionStationary rho mu B u ↔ ∀ x, eulerLagrange rho mu B u x = 0 := by
  have hc : ∀ i, Continuous (fun x => eulerLagrange rho mu B u x i) := fun i =>
    (S10Pilot.smooth_eulerLagrange hu rho mu i).continuous.add
      (S10Controls.smooth_eulerLagrange .divergenceOnly hu 0 B i).continuous
  constructor
  · intro hs
    have hz : ∀ i : Fin D, (fun x => eulerLagrange rho mu B u x i) = fun _ => 0 := by
      intro i
      have hae := ae_eq_zero_of_integral_contDiff_smul_eq_zero
        (μ := volume) (hc i).locallyIntegrable (fun f hf hcompact => ?_)
      · exact MeasureTheory.Measure.eq_of_ae_eq hae (hc i) continuous_const
      · have hh := single_testField hf hcompact i
        have h := hs (fun x => Pi.single i (f x)) hh
        rw [relativeAction_deriv_eq_eulerLagrange hu hh] at h
        simpa [Pi.single_apply, mul_comm] using h
    intro x
    ext i
    exact congrFun (hz i) x
  · intro he h hh
    rw [relativeAction_deriv_eq_eulerLagrange hu hh]
    simp [he]

theorem eulerLagrange_planeWave (rho mu B omega : ℝ) (k a : Vec D) (x : Point D) :
    eulerLagrange rho mu B (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator rho mu B omega k a := by
  rw [eulerLagrange, S10Pilot.eulerLagrange_planeWave,
    S10Controls.eulerLagrange_planeWave, modalOperator_split, smul_add]

theorem actionStationary_planeWave_iff (rho mu B omega : ℝ) (k a : Vec D) :
    ActionStationary rho mu B (planeWave omega k a) ↔ ModalStationary rho mu B omega k a := by
  rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a), modal_stationary_iff]
  constructor
  · intro h
    simpa [eulerLagrange_planeWave, phase, dot] using h 0
  · intro h x
    simp [eulerLagrange_planeWave, h]

theorem zero_compression_action (rho mu : ℝ) (J : Jet D) :
    lagrangian rho mu 0 J = S10Pilot.lagrangian rho mu J := by
  simp [lagrangian, S10Pilot.lagrangian]

theorem zero_compression_operator (rho mu omega : ℝ) (k a : Vec D) :
    modalOperator rho mu 0 omega k a = S10Pilot.modalOperator rho mu omega k a := by
  simp [modalOperator, S10Pilot.modalOperator_eq]

end
end S11Homogeneous
