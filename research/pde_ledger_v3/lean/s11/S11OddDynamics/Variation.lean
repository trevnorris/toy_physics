import S11OddDynamics.Action

/-! E2: actual finite relative action on smooth backgrounds, with compact test
fields. Instantiates the existing S10 analytic lemmas for the odd density. -/
namespace S11OddDynamics
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

def linearDensity (beta : ℝ) (J H : Jet 2) : ℝ :=
  ∑ i : Fin 2, ∑ j : Fin 3, momentum beta J j i * H j i

theorem linearDensity_eq (beta : ℝ) (J H : Jet 2) :
    linearDensity beta J H = variationDensity beta J H := by
  simp [linearDensity, momentum_eq, variationDensity, divergence, skew,
    Fin.sum_univ_succ]
  ring

theorem lagrangian_change (beta s : ℝ) (J H : Jet 2) :
    lagrangian beta (J + s • H) - lagrangian beta J =
      s * linearDensity beta J H + s ^ 2 * lagrangian beta H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (beta : ℝ) (J : Jet 2) :
    linearDensity beta J J = 2 * lagrangian beta J := by
  rw [linearDensity_eq]
  unfold variationDensity lagrangian
  ring

theorem integrable_linear_term {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) (i : Fin 2) (j : Fin 3) :
    Integrable (fun x => momentum beta (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum hu beta j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    Integrable (fun x => linearDensity beta (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh beta i j)

theorem integrable_test_lagrangian {h : Point 2 → Vec 2} (hh : TestField h) (beta : ℝ) :
    Integrable (fun x => lagrangian beta (fieldJet h x)) := by
  have hi := (integrable_linearDensity hh.1 hh beta).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (beta : ℝ) (u h : Point 2 → Vec 2) (s : ℝ) : ℝ :=
  ∫ x : Point 2, lagrangian beta (fieldJet (u + s • h) x) -
    lagrangian beta (fieldJet u x)

theorem relative_density_integrable {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta s : ℝ) :
    Integrable (fun x => lagrangian beta (fieldJet (u + s • h) x) -
      lagrangian beta (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  exact ((integrable_linearDensity hu hh beta).const_mul s).add
    ((integrable_test_lagrangian hh beta).const_mul (s ^ 2))

theorem relativeAction_expansion {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta s : ℝ) :
    relativeAction beta u h s =
      s * (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian beta (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  rw [integral_add ((integrable_linearDensity hu hh beta).const_mul s)
    ((integrable_test_lagrangian hh beta).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    HasDerivAt (relativeAction beta u h)
      (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) 0 := by
  have heq : relativeAction beta u h = fun s =>
      s * (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian beta (fieldJet h x)) :=
    funext (relativeAction_expansion hu hh beta)
  rw [heq]
  have hd := ((hasDerivAt_id (0 : ℝ)).mul_const
    (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (∫ x, lagrangian beta (fieldJet h x)))
  convert! hd using 1
  simp

theorem integrable_el_term {u h : Point 2 → Vec 2} (hu : SmoothField u) (hh : TestField h)
    (beta : ℝ) (i : Fin 2) (j : Fin 3) :
    Integrable (fun x => coordDeriv j (fun y => momentum beta (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum hu beta j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin 2, eulerLagrange beta u x i * h x i := by
  have hleft : ∀ i : Fin 2, Integrable (fun x =>
      ∑ j : Fin 3, momentum beta (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh beta i j)
  have hright : ∀ i : Fin 2, Integrable (fun x => eulerLagrange beta u x i * h x i) :=
    fun i => integrable_mul_compact (smooth_eulerLagrange hu beta i).continuous
      (hh.1 i).continuous (hh.2 i)
  unfold linearDensity
  rw [integral_finsetSum _ (fun i _ => hleft i),
    integral_finsetSum _ (fun i _ => hright i)]
  apply Finset.sum_congr rfl
  intro i _
  rw [integral_finsetSum _ (fun j _ => integrable_linear_term hu hh beta i j)]
  simp_rw [fieldJet, coord_integration_by_parts (smooth_momentum hu beta _ i)
    (hh.1 i) (hh.2 i)]
  simp only [eulerLagrange, neg_mul, Finset.sum_mul]
  rw [integral_neg, integral_finsetSum _ (fun j _ => integrable_el_term hu hh beta i j),
    Finset.sum_neg_distrib]

/-- The first variation equals the integral against the previously defined local EL expression. -/
theorem relativeAction_deriv_eq_eulerLagrange {u h : Point 2 → Vec 2}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    deriv (relativeAction beta u h) 0 =
      ∫ x, ∑ i : Fin 2, eulerLagrange beta u x i * h x i := by
  rw [(relativeAction_hasDerivAt hu hh beta).deriv]
  exact integrated_variation_by_parts hu hh beta

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (beta : ℝ) (u : Point 2 → Vec 2) : Prop :=
  ∀ h : Point 2 → Vec 2, TestField h → deriv (relativeAction beta u h) 0 = 0

/-- The integrated variational principle is equivalent to the pointwise local PDE.
Compact test functions give a.e. vanishing; smoothness and full-support Lebesgue
measure upgrade that conclusion to pointwise vanishing. -/
theorem actionStationary_iff_eulerLagrange {u : Point 2 → Vec 2}
    (hu : SmoothField u) (beta : ℝ) :
    ActionStationary beta u ↔ ∀ x, eulerLagrange beta u x = 0 := by
  constructor
  · intro hs
    have hzero : ∀ i : Fin 2, (fun x => eulerLagrange beta u x i) = fun _ => 0 := by
      intro i
      have hc := (smooth_eulerLagrange hu beta i).continuous
      have hae := ae_eq_zero_of_integral_contDiff_smul_eq_zero
        (μ := volume) hc.locallyIntegrable (fun f hf hcompact => ?_)
      · exact MeasureTheory.Measure.eq_of_ae_eq hae hc continuous_const
      · have hh := single_testField hf hcompact i
        have h := hs (fun x => Pi.single i (f x)) hh
        rw [relativeAction_deriv_eq_eulerLagrange hu hh] at h
        simpa [Pi.single_apply, mul_comm] using h
    intro x
    ext i
    exact congrFun (hzero i) x
  · intro he h hh
    rw [relativeAction_deriv_eq_eulerLagrange hu hh]
    simp [he]

end
end S11OddDynamics
