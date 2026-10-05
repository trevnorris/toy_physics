import S11D4Bulk.Action

/-! D4C.1: actual finite relative action on smooth backgrounds, with compact test
fields. Instantiates the existing S10 analytic lemmas for the full D4 invariant family. -/
namespace S11D4Bulk
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

def linearDensity (v : Vec 4) (J H : Jet 4) : ℝ :=
  ∑ i : Fin 4, ∑ j : Fin 5, momentum v J j i * H j i

theorem linearDensity_eq (v : Vec 4) (J H : Jet 4) :
    linearDensity v J H = variationDensity v J H := by
  unfold linearDensity
  simp_rw [momentum_eq, add_mul]
  simp only [Finset.sum_add_distrib]
  change (∑ i : Fin 4, ∑ j : Fin 5, Even.momentum v J j i * H j i) +
    S11D4Odd.linearDensity (v 3) J H = _
  rw [S11D4Odd.linearDensity_eq]
  simp only [variationDensity, Even.momentum_eq, Even.variationDensity,
    Even.divergence, Fin.sum_univ_succ, Fin.sum_univ_zero,
    Fin.cases_zero, Fin.cases_succ]
  norm_num [Fin.succ]
  ring!


theorem lagrangian_change (v : Vec 4) (s : ℝ) (J H : Jet 4) :
    lagrangian v (J + s • H) - lagrangian v J =
      s * linearDensity v J H + s ^ 2 * lagrangian v H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (v : Vec 4) (J : Jet 4) :
    linearDensity v J J = 2 * lagrangian v J := by
  rw [linearDensity_eq]
  simp only [variationDensity, lagrangian, Even.variationDensity, Even.lagrangian,
    Even.divergence, Even.transposePair, Even.gradientSquare,
    S11D4Odd.variationDensity, S11D4Odd.lagrangian, S11D4Odd.orientationJet_eq,
    Fin.sum_univ_four]
  ring

theorem integrable_linear_term {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) (i : Fin 4) (j : Fin 5) :
    Integrable (fun x => momentum v (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum hu v j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) :
    Integrable (fun x => linearDensity v (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh v i j)

theorem integrable_test_lagrangian {h : Point 4 → Vec 4} (hh : TestField h) (v : Vec 4) :
    Integrable (fun x => lagrangian v (fieldJet h x)) := by
  have hi := (integrable_linearDensity hh.1 hh v).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (v : Vec 4) (u h : Point 4 → Vec 4) (s : ℝ) : ℝ :=
  ∫ x : Point 4, lagrangian v (fieldJet (u + s • h) x) -
    lagrangian v (fieldJet u x)

theorem relative_density_integrable {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) (s : ℝ) :
    Integrable (fun x => lagrangian v (fieldJet (u + s • h) x) -
      lagrangian v (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  exact ((integrable_linearDensity hu hh v).const_mul s).add
    ((integrable_test_lagrangian hh v).const_mul (s ^ 2))

theorem relativeAction_expansion {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) (s : ℝ) :
    relativeAction v u h s =
      s * (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian v (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  rw [integral_add ((integrable_linearDensity hu hh v).const_mul s)
    ((integrable_test_lagrangian hh v).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) :
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

theorem integrable_el_term {u h : Point 4 → Vec 4} (hu : SmoothField u) (hh : TestField h)
    (v : Vec 4) (i : Fin 4) (j : Fin 5) :
    Integrable (fun x => coordDeriv j (fun y => momentum v (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum hu v j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) :
    (∫ x, linearDensity v (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin 4, eulerLagrange v u x i * h x i := by
  have hleft : ∀ i : Fin 4, Integrable (fun x =>
      ∑ j : Fin 5, momentum v (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh v i j)
  have hright : ∀ i : Fin 4, Integrable (fun x => eulerLagrange v u x i * h x i) :=
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
theorem relativeAction_deriv_eq_eulerLagrange {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (v : Vec 4) :
    deriv (relativeAction v u h) 0 =
      ∫ x, ∑ i : Fin 4, eulerLagrange v u x i * h x i := by
  rw [(relativeAction_hasDerivAt hu hh v).deriv]
  exact integrated_variation_by_parts hu hh v

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (v : Vec 4) (u : Point 4 → Vec 4) : Prop :=
  ∀ h : Point 4 → Vec 4, TestField h → deriv (relativeAction v u h) 0 = 0

/-- The integrated variational principle is equivalent to the pointwise local PDE.
Compact test functions give a.e. vanishing; smoothness and full-support Lebesgue
measure upgrade that conclusion to pointwise vanishing. -/
theorem actionStationary_iff_eulerLagrange {u : Point 4 → Vec 4}
    (hu : SmoothField u) (v : Vec 4) :
    ActionStationary v u ↔ ∀ x, eulerLagrange v u x = 0 := by
  constructor
  · intro hs
    have hzero : ∀ i : Fin 4, (fun x => eulerLagrange v u x i) = fun _ => 0 := by
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
end S11D4Bulk
