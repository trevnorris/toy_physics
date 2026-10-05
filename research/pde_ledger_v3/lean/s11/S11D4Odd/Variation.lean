import S11D4Odd.Action

/-! D4B.1: actual relative action with smooth compact variations for the odd density. -/
namespace S11D4Odd
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

def linearDensity (beta : ℝ) (J H : Jet 4) : ℝ :=
  ∑ i : Fin 4, ∑ j : Fin 5, momentum beta J j i * H j i

theorem linearDensity_eq (beta : ℝ) (J H : Jet 4) :
    linearDensity beta J H = variationDensity beta J H := by
  simp only [linearDensity, momentum_eq, Fin.sum_univ_succ, Fin.sum_univ_zero,
    Fin.cases_zero, Fin.cases_succ, dualCurl, Matrix.of_apply, Matrix.cons_val,
    variationDensity, antisym]
  norm_num [Fin.succ, Matrix.cons_val_two, Matrix.cons_val_three]
  ring!

theorem lagrangian_change (beta : ℝ) (s : ℝ) (J H : Jet 4) :
    lagrangian beta (J + s • H) - lagrangian beta J =
      s * linearDensity beta J H + s ^ 2 * lagrangian beta H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (beta : ℝ) (J : Jet 4) :
    linearDensity beta J J = 2 * lagrangian beta J := by
  rw [linearDensity_eq]
  simp only [variationDensity, lagrangian, orientationJet_eq]
  ring

theorem integrable_linear_term {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) (i : Fin 4) (j : Fin 5) :
    Integrable (fun x => momentum beta (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum hu beta j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    Integrable (fun x => linearDensity beta (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh beta i j)

theorem integrable_test_lagrangian {h : Point 4 → Vec 4} (hh : TestField h) (beta : ℝ) :
    Integrable (fun x => lagrangian beta (fieldJet h x)) := by
  have hi := (integrable_linearDensity hh.1 hh beta).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (beta : ℝ) (u h : Point 4 → Vec 4) (s : ℝ) : ℝ :=
  ∫ x : Point 4, lagrangian beta (fieldJet (u + s • h) x) -
    lagrangian beta (fieldJet u x)

theorem relative_density_integrable {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) (s : ℝ) :
    Integrable (fun x => lagrangian beta (fieldJet (u + s • h) x) -
      lagrangian beta (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  exact ((integrable_linearDensity hu hh beta).const_mul s).add
    ((integrable_test_lagrangian hh beta).const_mul (s ^ 2))

theorem relativeAction_expansion {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) (s : ℝ) :
    relativeAction beta u h s =
      s * (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian beta (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  rw [integral_add ((integrable_linearDensity hu hh beta).const_mul s)
    ((integrable_test_lagrangian hh beta).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt {u h : Point 4 → Vec 4}
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

theorem integrable_el_term {u h : Point 4 → Vec 4} (hu : SmoothField u) (hh : TestField h)
    (beta : ℝ) (i : Fin 4) (j : Fin 5) :
    Integrable (fun x => coordDeriv j (fun y => momentum beta (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum hu beta j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    (∫ x, linearDensity beta (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin 4, eulerLagrange beta u x i * h x i := by
  have hleft : ∀ i : Fin 4, Integrable (fun x =>
      ∑ j : Fin 5, momentum beta (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh beta i j)
  have hright : ∀ i : Fin 4, Integrable (fun x => eulerLagrange beta u x i * h x i) :=
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
theorem relativeAction_deriv_eq_eulerLagrange {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (beta : ℝ) :
    deriv (relativeAction beta u h) 0 =
      ∫ x, ∑ i : Fin 4, eulerLagrange beta u x i * h x i := by
  rw [(relativeAction_hasDerivAt hu hh beta).deriv]
  exact integrated_variation_by_parts hu hh beta

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (beta : ℝ) (u : Point 4 → Vec 4) : Prop :=
  ∀ h : Point 4 → Vec 4, TestField h → deriv (relativeAction beta u h) 0 = 0

end
end S11D4Odd
