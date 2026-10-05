import S11D5Bulk.Action

/-! D5B.1: actual finite relative action on smooth backgrounds, with compact test
fields. Instantiates the existing S10 analytic lemmas for the D5 invariant family. -/
namespace S11D5Bulk
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

def linearDensity (v : Coeff) (J H : Jet 5) : ℝ :=
  ∑ i : Fin 5, ∑ j : Fin 6, momentum v J j i * H j i

theorem linearDensity_eq (v : Coeff) (J H : Jet 5) :
    linearDensity v J H = variationDensity v J H := by
  simp [linearDensity, momentum_eq, variationDensity, divergence,
    Fin.sum_univ_succ]
  ring

theorem lagrangian_change (v : Coeff) (s : ℝ) (J H : Jet 5) :
    lagrangian v (J + s • H) - lagrangian v J =
      s * linearDensity v J H + s ^ 2 * lagrangian v H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (v : Coeff) (J : Jet 5) :
    linearDensity v J J = 2 * lagrangian v J := by
  rw [linearDensity_eq]
  simp only [variationDensity, lagrangian, divergence, transposePair, gradientSquare, S11D5Invariants.sum_five]
  ring

theorem integrable_linear_term {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) (i : Fin 5) (j : Fin 6) :
    Integrable (fun x => momentum v (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum hu v j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) :
    Integrable (fun x => linearDensity v (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh v i j)

theorem integrable_test_lagrangian {h : Point 5 → Vec 5} (hh : TestField h) (v : Coeff) :
    Integrable (fun x => lagrangian v (fieldJet h x)) := by
  have hi := (integrable_linearDensity hh.1 hh v).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (v : Coeff) (u h : Point 5 → Vec 5) (s : ℝ) : ℝ :=
  ∫ x : Point 5, lagrangian v (fieldJet (u + s • h) x) -
    lagrangian v (fieldJet u x)

theorem relative_density_integrable {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) (s : ℝ) :
    Integrable (fun x => lagrangian v (fieldJet (u + s • h) x) -
      lagrangian v (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  exact ((integrable_linearDensity hu hh v).const_mul s).add
    ((integrable_test_lagrangian hh v).const_mul (s ^ 2))

theorem relativeAction_expansion {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) (s : ℝ) :
    relativeAction v u h s =
      s * (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian v (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  rw [integral_add ((integrable_linearDensity hu hh v).const_mul s)
    ((integrable_test_lagrangian hh v).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) :
    HasDerivAt (relativeAction v u h)
      (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) 0 := by
  have heq : relativeAction v u h = fun s =>
      s * (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian v (fieldJet h x)) :=
    funext (relativeAction_expansion hu hh v)
  rw [heq]
  have hd := ((hasDerivAt_id (0 : ℝ)).mul_const
    (∫ x, linearDensity v (fieldJet u x) (fieldJet h x))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (∫ x, lagrangian v (fieldJet h x)))
  convert! hd using 1
  simp

theorem integrable_el_term {u h : Point 5 → Vec 5} (hu : SmoothField u) (hh : TestField h)
    (v : Coeff) (i : Fin 5) (j : Fin 6) :
    Integrable (fun x => coordDeriv j (fun y => momentum v (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum hu v j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) :
    (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin 5, eulerLagrange v u x i * h x i := by
  have hleft : ∀ i : Fin 5, Integrable (fun x =>
      ∑ j : Fin 6, momentum v (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh v i j)
  have hright : ∀ i : Fin 5, Integrable (fun x => eulerLagrange v u x i * h x i) :=
    fun i => integrable_mul_compact (smooth_eulerLagrange hu v i).continuous
      (hh.1 i).continuous (hh.2 i)
  unfold linearDensity
  rw [integral_finsetSum _ (fun i _ => hleft i),
    integral_finsetSum _ (fun i _ => hright i)]
  apply Finset.sum_congr rfl
  intro i _
  rw [integral_finsetSum _ (fun j _ => integrable_linear_term hu hh v i j)]
  simp_rw [fieldJet, coord_integration_by_parts (smooth_momentum hu v _ i)
    (hh.1 i) (hh.2 i)]
  simp only [eulerLagrange, neg_mul, Finset.sum_mul]
  rw [integral_neg, integral_finsetSum _ (fun j _ => integrable_el_term hu hh v i j),
    Finset.sum_neg_distrib]

/-- The first variation equals the integral against the previously defined local EL expression. -/
theorem relativeAction_deriv_eq_eulerLagrange {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (v : Coeff) :
    deriv (relativeAction v u h) 0 =
      ∫ x, ∑ i : Fin 5, eulerLagrange v u x i * h x i := by
  rw [(relativeAction_hasDerivAt hu hh v).deriv]
  exact integrated_variation_by_parts hu hh v

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (v : Coeff) (u : Point 5 → Vec 5) : Prop :=
  ∀ h : Point 5 → Vec 5, TestField h → deriv (relativeAction v u h) 0 = 0

/-- The integrated variational principle is equivalent to the pointwise local PDE.
Compact test functions give a.e. vanishing; smoothness and full-support Lebesgue
measure upgrade that conclusion to pointwise vanishing. -/
theorem actionStationary_iff_eulerLagrange {u : Point 5 → Vec 5}
    (hu : SmoothField u) (v : Coeff) :
    ActionStationary v u ↔ ∀ x, eulerLagrange v u x = 0 := by
  constructor
  · intro hs
    have hzero : ∀ i : Fin 5, (fun x => eulerLagrange v u x i) = fun _ => 0 := by
      intro i
      have hc := (smooth_eulerLagrange hu v i).continuous
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
end S11D5Bulk
