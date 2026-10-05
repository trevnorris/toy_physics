import S11VariableCoefficients.Common
import S11D3Bulk.Bulk

/-! VC1: prescribed profiles in the same derivative-defined D3 action. -/
namespace S11VariableCoefficients.D3
noncomputable section
open S10Pilot
open scoped ContDiff

def lagrangian (v : Point 3 → Vec 3) (x : Point 3) (J : Jet 3) : ℝ :=
  S11D3Bulk.lagrangian (v x) J

def momentum (v : Point 3 → Vec 3) (x : Point 3) (J : Jet 3) (j : Fin 4) (i : Fin 3) : ℝ :=
  deriv (fun s : ℝ => lagrangian v x (J + s • basisJet j i)) 0

theorem momentum_identity (v : Point 3 → Vec 3) (x : Point 3) (J : Jet 3) (j : Fin 4) (i : Fin 3) :
    momentum v x J j i = S11D3Bulk.momentum (v x) J j i := rfl

theorem pointwise_variation (v : Point 3 → Vec 3) (x : Point 3) (J H : Jet 3) :
    HasDerivAt (fun s : ℝ => lagrangian v x (J + s • H))
      (S11D3Bulk.variationDensity (v x) J H) 0 := S11D3Bulk.lagrangian_variation _ _ _

def eulerLagrange (v : Point 3 → Vec 3) (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun i => -∑ j : Fin 4, coordDeriv j (fun y => momentum v y (fieldJet u y) j i) x

def correction (v : Point 3 → Vec 3) (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun j => coordDeriv j.succ (fun y => v y 0) x * S11D3Bulk.divergence (fieldJet u x) +
    ∑ i : Fin 3, (coordDeriv i.succ (fun y => v y 1) x * fieldJet u x j.succ i +
      coordDeriv i.succ (fun y => v y 2) x * fieldJet u x i.succ j)

theorem eulerLagrange_product {v : Point 3 → Vec 3} {u : Point 3 → Vec 3}
    (hv : ∀ i, ContDiff ℝ ∞ (fun x => v x i)) (hu : SmoothField u) (x : Point 3) :
    eulerLagrange v u x = S11D3Bulk.eulerLagrange (v x) u x + correction v u x := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  ext i
  fin_cases i
  all_goals simp only [eulerLagrange, momentum_identity, S11D3Bulk.eulerLagrange,
    S11D3Bulk.momentum_eq, correction, S11D3Bulk.divergence,
    Fin.sum_univ_succ, Fin.sum_univ_zero, Fin.cases_zero, Fin.cases_succ, Pi.add_apply]
  all_goals norm_num [Fin.succ] at *
  all_goals simp (disch := fun_prop) only [S11D4Odd.partial_neg, S11D4Odd.partial_add,
    S11D4Odd.partial_const_mul, S11D4Odd.partial_mul, S11D4Odd.partial_zero]
  all_goals ring!

theorem eulerLagrange_eq {v : Point 3 → Vec 3} {u : Point 3 → Vec 3}
    (hv : ∀ i, ContDiff ℝ ∞ (fun x => v x i)) (hu : SmoothField u) (x : Point 3) :
    eulerLagrange v u x =
      (v x 0 + v x 1) • S11D3Bulk.gradDiv u x +
      v x 2 • S11D3Bulk.laplacian u x + correction v u x := by
  rw [eulerLagrange_product hv hu, S11D3Bulk.eulerLagrange_eq hu]

theorem constant_profile (v : Vec 3) (u : Point 3 → Vec 3) (x : Point 3) :
    eulerLagrange (fun _ => v) u x = S11D3Bulk.eulerLagrange v u x := rfl

theorem null_profile_residual {a : Point 3 → ℝ} {u : Point 3 → Vec 3}
    (ha : ContDiff ℝ ∞ a) (hu : SmoothField u) (x : Point 3) (j : Fin 3) :
    eulerLagrange (fun y => ![a y, -a y, 0]) u x j =
      coordDeriv j.succ a x * S11D3Bulk.divergence (fieldJet u x) -
        ∑ i : Fin 3, coordDeriv i.succ a x * fieldJet u x j.succ i := by
  have hv : ∀ i : Fin 3, ContDiff ℝ ∞ (fun y => (![a y, -a y, 0] : Vec 3) i) := by
    intro i
    fin_cases i
    · exact ha
    · exact ha.neg
    · exact contDiff_const
  rw [eulerLagrange_eq hv hu]
  simp only [Matrix.cons_val, add_neg_cancel, zero_smul, Pi.add_apply, Pi.zero_apply, zero_add,
    correction, S11D4Odd.partial_neg ha, S11D4Odd.partial_zero,
    zero_mul, add_zero, neg_mul, Finset.sum_neg_distrib]
  ring

theorem weighted_null_density {a : Point 3 → ℝ} {u : Point 3 → Vec 3}
    (ha : ContDiff ℝ ∞ a) (hu : SmoothField u) (x : Point 3) :
    lagrangian (fun y => ![a y,-a y,0]) x (fieldJet u x) =
      -(1 / 2 : ℝ) * spatialDiv (fun y i => a y * S11D3Bulk.boundaryCurrent u y i) x +
        (1 / 2 : ℝ) * gradientPair a (S11D3Bulk.boundaryCurrent u) x := by
  have hJ : ∀ i, ContDiff ℝ ∞ (fun y => S11D3Bulk.boundaryCurrent u y i) := by
    have hU : ∀ i, ContDiff ℝ ∞ (fun x => u x i) := hu
    have hG : ∀ i j, ContDiff ℝ ∞ (fun x => fieldJet u x i j) :=
      fun i j => smooth_coordDeriv (hu j) i
    intro i
    unfold S11D3Bulk.boundaryCurrent
    fun_prop
  rw [weighted_divergence ha hJ]
  have he := S11D3Bulk.null_density_is_divergence hu (a x) x
  change S11D3Bulk.lagrangian _ _ = _
  rw [he]
  unfold spatialDiv
  ring

def traction (v : Vec 3) (J : Jet 3) (n : Vec 3) : Vec 3 :=
  normalFlux n (fun i j => S11D3Bulk.momentum v J i.succ j)

theorem traction_eq (v : Vec 3) (J : Jet 3) (n : Vec 3) (j : Fin 3) :
    traction v J n j = -∑ i : Fin 3, n i *
      (v 0 * (if i = j then S11D3Bulk.divergence J else 0) + v 1 * J j.succ i + v 2 * J i.succ j) := by
  simp only [traction, normalFlux, S11D3Bulk.momentum_eq, Fin.cases_succ]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

end
end S11VariableCoefficients.D3
