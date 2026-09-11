import S10Pilot.Analytic

/-!
Integrated variation of the S10 action in arbitrary dimension. A smooth background need not decay.
Only the variation has compact support; the integral is of the pointwise change
in density, not a subtraction of two potentially infinite action integrals.
-/

namespace S10Pilot
noncomputable section
open MeasureTheory
open scoped ContDiff

variable {D : ℕ}

def linearDensity (rho mu : ℝ) (J H : Jet D) : ℝ :=
  ∑ i : Fin D, ∑ j : Fin (D + 1), momentum rho mu J j i * H j i

theorem linearDensity_eq (rho mu : ℝ) (J H : Jet D) :
    linearDensity rho mu J H = variationDensity rho mu J H := by
  have hs : (∑ i, ∑ j, antisym J i j * H j.succ i) =
      -(∑ i, ∑ j, antisym J i j * H i.succ j) := by
    rw [Finset.sum_comm]
    simp only [← Finset.sum_neg_distrib]
    apply Finset.sum_congr rfl
    intro i _
    apply Finset.sum_congr rfl
    intro j _
    unfold antisym
    ring
  have hp : (∑ i, ∑ j, antisym J i j * antisym H i j) =
      2 * ∑ i, ∑ j, antisym J i j * H i.succ j := by
    change (∑ i, ∑ j, antisym J i j * (H i.succ j - H j.succ i)) = _
    simp only [mul_sub, Finset.sum_sub_distrib]
    rw [hs]
    ring
  have ht : ∀ i : Fin D, (∑ j : Fin (D + 1), momentum rho mu J j i * H j i) =
      rho * J 0 i * H 0 i - mu * ∑ j, antisym J j i * H j.succ i := by
    intro i
    simp only [Fin.sum_univ_succ, momentum_eq, Fin.cases_zero, Fin.cases_succ,
      neg_mul, mul_assoc, Finset.sum_neg_distrib, ← Finset.mul_sum]
    ring
  unfold linearDensity variationDensity
  simp_rw [ht]
  rw [hp]
  simp only [Finset.sum_sub_distrib, mul_assoc, ← Finset.mul_sum]
  rw [Finset.sum_comm (f := fun i j : Fin D => antisym J j i * H j.succ i)]
  unfold dot
  ring

theorem lagrangian_change (rho mu s : ℝ) (J H : Jet D) :
    lagrangian rho mu (J + s • H) - lagrangian rho mu J =
      s * linearDensity rho mu J H + s ^ 2 * lagrangian rho mu H := by
  rw [lagrangian_increment, linearDensity_eq]
  ring

theorem linearDensity_self (rho mu : ℝ) (J : Jet D) :
    linearDensity rho mu J J = 2 * lagrangian rho mu J := by
  simp only [linearDensity_eq, variationDensity, lagrangian, normSq, stiffness, sq]
  ring

theorem integrable_linear_term {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) (i : Fin D) (j : Fin (D + 1)) :
    Integrable (fun x => momentum rho mu (fieldJet u x) j i * fieldJet h x j i) :=
  integrable_mul_compact (smooth_momentum hu rho mu j i).continuous
    (smooth_coordDeriv (hh.1 i) j).continuous (compact_coordDeriv (hh.1 i) (hh.2 i) j)

theorem integrable_linearDensity {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    Integrable (fun x => linearDensity rho mu (fieldJet u x) (fieldJet h x)) := by
  apply integrable_finsetSum
  intro i _
  exact integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh rho mu i j)

theorem integrable_test_lagrangian {h : Point D → Vec D} (hh : TestField h) (rho mu : ℝ) :
    Integrable (fun x => lagrangian rho mu (fieldJet h x)) := by
  have hi := (integrable_linearDensity hh.1 hh rho mu).const_mul (1 / 2 : ℝ)
  convert! hi using 1
  funext x
  rw [linearDensity_self]
  ring

/-- Finite change of action under a compactly supported variation. -/
def relativeAction (rho mu : ℝ) (u h : Point D → Vec D) (s : ℝ) : ℝ :=
  ∫ x : Point D, lagrangian rho mu (fieldJet (u + s • h) x) -
    lagrangian rho mu (fieldJet u x)

theorem relative_density_integrable {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ) :
    Integrable (fun x => lagrangian rho mu (fieldJet (u + s • h) x) -
      lagrangian rho mu (fieldJet u x)) := by
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  exact ((integrable_linearDensity hu hh rho mu).const_mul s).add
    ((integrable_test_lagrangian hh rho mu).const_mul (s ^ 2))

theorem relativeAction_expansion {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ) :
    relativeAction rho mu u h s =
      s * (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian rho mu (fieldJet h x)) := by
  unfold relativeAction
  simp_rw [fieldJet_add_smul hu hh.1, lagrangian_change]
  rw [integral_add ((integrable_linearDensity hu hh rho mu).const_mul s)
    ((integrable_test_lagrangian hh rho mu).const_mul (s ^ 2))]
  simp only [integral_const_mul]

/-- Differentiate the actual integral using its proved quadratic expansion. -/
theorem relativeAction_hasDerivAt {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    HasDerivAt (relativeAction rho mu u h)
      (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x)) 0 := by
  have heq : relativeAction rho mu u h = fun s =>
      s * (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x)) +
      s ^ 2 * (∫ x, lagrangian rho mu (fieldJet h x)) :=
    funext (relativeAction_expansion hu hh rho mu)
  rw [heq]
  have hd := ((hasDerivAt_id (0 : ℝ)).mul_const
    (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (∫ x, lagrangian rho mu (fieldJet h x)))
  convert! hd using 1
  simp

theorem integrable_el_term {u h : Point D → Vec D} (hu : SmoothField u) (hh : TestField h)
    (rho mu : ℝ) (i : Fin D) (j : Fin (D + 1)) :
    Integrable (fun x => coordDeriv j (fun y => momentum rho mu (fieldJet u y) j i) x *
      h x i) :=
  integrable_mul_compact
    (smooth_coordDeriv (smooth_momentum hu rho mu j i) j).continuous
    (hh.1 i).continuous (hh.2 i)

theorem integrated_variation_by_parts {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x)) =
      ∫ x, ∑ i : Fin D, eulerLagrange rho mu u x i * h x i := by
  have hleft : ∀ i : Fin D, Integrable (fun x =>
      ∑ j : Fin (D + 1), momentum rho mu (fieldJet u x) j i * fieldJet h x j i) :=
    fun i => integrable_finsetSum _ (fun j _ => integrable_linear_term hu hh rho mu i j)
  have hright : ∀ i : Fin D, Integrable (fun x => eulerLagrange rho mu u x i * h x i) :=
    fun i => integrable_mul_compact (smooth_eulerLagrange hu rho mu i).continuous
      (hh.1 i).continuous (hh.2 i)
  unfold linearDensity
  rw [integral_finsetSum _ (fun i _ => hleft i),
    integral_finsetSum _ (fun i _ => hright i)]
  apply Finset.sum_congr rfl
  intro i _
  rw [integral_finsetSum _ (fun j _ => integrable_linear_term hu hh rho mu i j)]
  simp_rw [fieldJet, coord_integration_by_parts (smooth_momentum hu rho mu _ i)
    (hh.1 i) (hh.2 i)]
  simp only [eulerLagrange, neg_mul, Finset.sum_mul]
  rw [integral_neg, integral_finsetSum _ (fun j _ => integrable_el_term hu hh rho mu i j),
    Finset.sum_neg_distrib]

/-- The first variation equals the integral against the previously defined local EL expression. -/
theorem relativeAction_deriv_eq_eulerLagrange {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ) :
    deriv (relativeAction rho mu u h) 0 =
      ∫ x, ∑ i : Fin D, eulerLagrange rho mu u x i * h x i := by
  rw [(relativeAction_hasDerivAt hu hh rho mu).deriv]
  exact integrated_variation_by_parts hu hh rho mu

/-- Stationarity against all smooth, compactly supported vector variations. -/
def ActionStationary (rho mu : ℝ) (u : Point D → Vec D) : Prop :=
  ∀ h : Point D → Vec D, TestField h → deriv (relativeAction rho mu u h) 0 = 0

theorem single_testField {f : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (i : Fin D) :
    TestField (fun x => Pi.single i (f x)) := by
  constructor
  · intro j
    by_cases hij : i = j
    · subst j
      simpa using hf
    · simpa [Pi.single_apply, hij] using
        (contDiff_const : ContDiff ℝ ∞ (fun _ : Point D => (0 : ℝ)))
  · intro j
    by_cases hij : i = j
    · subst j
      simpa using hc
    · have heq : (fun x : Point D => (Pi.single i (f x) : Vec D) j) = (0 : Point D → ℝ) := by
        funext x
        simp [hij]
      rw [heq]
      exact HasCompactSupport.zero

/-- The integrated variational principle is equivalent to the pointwise local PDE.
Compact test functions give a.e. vanishing; smoothness and full-support Lebesgue
measure upgrade that conclusion to pointwise vanishing. -/
theorem actionStationary_iff_eulerLagrange {u : Point D → Vec D}
    (hu : SmoothField u) (rho mu : ℝ) :
    ActionStationary rho mu u ↔ ∀ x, eulerLagrange rho mu u x = 0 := by
  constructor
  · intro hs
    have hzero : ∀ i : Fin D, (fun x => eulerLagrange rho mu u x i) = fun _ => 0 := by
      intro i
      have hc := (smooth_eulerLagrange hu rho mu i).continuous
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

theorem smooth_planeWave (omega : ℝ) (k a : Vec D) : SmoothField (planeWave omega k a) := by
  intro i
  unfold planeWave phase dot
  fun_prop

/-- The earlier modal condition now characterizes stationarity against arbitrary
compact test fields, including variations outside the plane-wave ansatz. -/
theorem actionStationary_planeWave_iff (rho mu omega : ℝ) (k a : Vec D) :
    ActionStationary rho mu (planeWave omega k a) ↔ ModalStationary rho mu omega k a := by
  rw [actionStationary_iff_eulerLagrange (smooth_planeWave omega k a)]
  exact planeWave_solves_iff rho mu omega k a

end
end S10Pilot
