import S10Anisotropic.Analytic

/-! The finite action change and its integrated first variation for the anisotropic action.
The background is smooth; variations are smooth and compactly supported. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
open MeasureTheory
open scoped ContDiff
variable {D : ℕ}

def linearDensity (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (J H : Jet D) : ℝ :=
  ∑ i : Fin D, ∑ j : Fin (D + 1), momentum e sigma rho mu J j i * H j i

theorem linearDensity_eq (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (J H : Jet D) :
    linearDensity e sigma rho mu J H = variationDensity e sigma rho mu J H := by
  unfold linearDensity
  simp only [momentum_eq, add_mul, Finset.sum_add_distrib]
  change S10Pilot.linearDensity rho mu J H + _ = _
  rw [S10Pilot.linearDensity_eq]
  simp [variationDensity, ite_mul, mul_ite]

theorem lagrangian_change (e : Fin D) (sigma : ℝ) (rho mu s : ℝ) (J H : Jet D) :
    lagrangian e sigma rho mu (J + s • H) - lagrangian e sigma rho mu J =
      s * linearDensity e sigma rho mu J H + s ^ 2 * lagrangian e sigma rho mu H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (J : Jet D) :
    linearDensity e sigma rho mu J J = 2 * lagrangian e sigma rho mu J := by
  rw [linearDensity_eq, variationDensity, lagrangian_eq,
    ← S10Pilot.linearDensity_eq, S10Pilot.linearDensity_self]
  ring

theorem integrable_linear_term (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) (i : Fin D) (j : Fin (D + 1)) :
    Integrable (fun x => momentum e sigma rho mu (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum e sigma hu rho mu j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    Integrable (fun x => linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term e sigma hu hh rho mu i j)

theorem integrable_test_lagrangian (e : Fin D) (sigma : ℝ) {h : Point D → Vec D} (hh : TestField h) (rho mu : ℝ) :
    Integrable (fun x => lagrangian e sigma rho mu (fieldJet h x)) := by
  have hi := (integrable_linearDensity e sigma hh.1 hh rho mu).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self e sigma]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (u h : Point D → Vec D) (s : ℝ) : ℝ :=
  ∫ x : Point D, lagrangian e sigma rho mu (fieldJet (u + s • h) x) -
    lagrangian e sigma rho mu (fieldJet u x)

theorem relative_density_integrable (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ) :
    Integrable (fun x => lagrangian e sigma rho mu (fieldJet (u + s • h) x) -
      lagrangian e sigma rho mu (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change e sigma]
  exact ((integrable_linearDensity e sigma hu hh rho mu).const_mul s).add
    ((integrable_test_lagrangian e sigma hh rho mu).const_mul (s ^ 2))

theorem relativeAction_expansion (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ) :
    relativeAction e sigma rho mu u h s =
      s * (∫ x, linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian e sigma rho mu (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change e sigma]
  rw [integral_add ((integrable_linearDensity e sigma hu hh rho mu).const_mul s)
    ((integrable_test_lagrangian e sigma hh rho mu).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    HasDerivAt (relativeAction e sigma rho mu u h)
      (∫ x, linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x)) 0 := by
  have heq : relativeAction e sigma rho mu u h = fun s =>
      s * (∫ x, linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian e sigma rho mu (fieldJet h x)) :=
    funext (relativeAction_expansion e sigma hu hh rho mu)
  rw [heq]
  have hd := ((hasDerivAt_id (0 : ℝ)).mul_const
    (∫ x, linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (∫ x, lagrangian e sigma rho mu (fieldJet h x)))
  convert! hd using 1
  simp

theorem integrable_el_term (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D} (hu : SmoothField u) (hh : TestField h)
    (rho mu : ℝ) (i : Fin D) (j : Fin (D + 1)) :
    Integrable (fun x => coordDeriv j (fun y => momentum e sigma rho mu (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum e sigma hu rho mu j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    (∫ x, linearDensity e sigma rho mu (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin D, eulerLagrange e sigma rho mu u x i * h x i := by
  have hleft : ∀ i : Fin D, Integrable (fun x =>
      ∑ j : Fin (D + 1), momentum e sigma rho mu (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term e sigma hu hh rho mu i j)
  have hright : ∀ i : Fin D, Integrable (fun x => eulerLagrange e sigma rho mu u x i * h x i) :=
    fun i => integrable_mul_compact (smooth_eulerLagrange e sigma hu rho mu i).continuous
      (hh.1 i).continuous (hh.2 i)
  unfold linearDensity
  rw [integral_finsetSum _ (fun i _ => hleft i),
    integral_finsetSum _ (fun i _ => hright i)]
  apply Finset.sum_congr rfl
  intro i _
  rw [integral_finsetSum _ (fun j _ => integrable_linear_term e sigma hu hh rho mu i j)]
  simp_rw [fieldJet, coord_integration_by_parts (smooth_momentum e sigma hu rho mu _ i)
    (hh.1 i) (hh.2 i)]
  simp only [eulerLagrange, neg_mul, Finset.sum_mul]
  rw [integral_neg, integral_finsetSum _ (fun j _ => integrable_el_term e sigma hu hh rho mu i j),
    Finset.sum_neg_distrib]

/-- The first variation equals the integral against the previously defined local EL expression. -/
theorem relativeAction_deriv_eq_eulerLagrange (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    deriv (relativeAction e sigma rho mu u h) 0 =
      ∫ x, ∑ i : Fin D, eulerLagrange e sigma rho mu u x i * h x i := by
  rw [(relativeAction_hasDerivAt e sigma hu hh rho mu).deriv]
  exact integrated_variation_by_parts e sigma hu hh rho mu

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (u : Point D → Vec D) : Prop :=
  ∀ h : Point D → Vec D, TestField h → deriv (relativeAction e sigma rho mu u h) 0 = 0

/-- The integrated variational principle is equivalent to the pointwise local PDE.
Compact test functions give a.e. vanishing; smoothness and full-support Lebesgue
measure upgrade that conclusion to pointwise vanishing. -/
theorem actionStationary_iff_eulerLagrange (e : Fin D) (sigma : ℝ) {u : Point D → Vec D}
    (hu : SmoothField u) (rho mu : ℝ) :
    ActionStationary e sigma rho mu u ↔ ∀ x, eulerLagrange e sigma rho mu u x = 0 := by
  constructor
  · intro hs
    have hzero : ∀ i : Fin D, (fun x => eulerLagrange e sigma rho mu u x i) = fun _ => 0 := by
      intro i
      have hc := (smooth_eulerLagrange e sigma hu rho mu i).continuous
      have hae := ae_eq_zero_of_integral_contDiff_smul_eq_zero
        (μ := volume) hc.locallyIntegrable (fun f hf hcompact => ?_)
      · exact MeasureTheory.Measure.eq_of_ae_eq hae hc continuous_const
      · have hh := single_testField hf hcompact i
        have h := hs (fun x => Pi.single i (f x)) hh
        rw [relativeAction_deriv_eq_eulerLagrange e sigma hu hh] at h
        simpa [Pi.single_apply, mul_comm] using h
    intro x
    ext i
    exact congrFun (hzero i) x
  · intro he h hh
    rw [relativeAction_deriv_eq_eulerLagrange e sigma hu hh]
    simp [he]

/-- The earlier modal condition now characterizes stationarity against arbitrary
compact test fields, including variations outside the plane-wave ansatz. -/
theorem actionStationary_planeWave_iff (e : Fin D) (sigma : ℝ) (rho mu omega : ℝ) (k a : Vec D) :
    ActionStationary e sigma rho mu (planeWave omega k a) ↔ ModalStationary e sigma rho mu omega k a := by
  rw [actionStationary_iff_eulerLagrange e sigma (smooth_planeWave omega k a)]
  exact planeWave_solves_iff e sigma rho mu omega k a

end
end S10Anisotropic
