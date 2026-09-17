import S11VariableCoefficients.Common
import S11D4Odd.Boundary

/-! VC2: the formerly bulk-null D4 odd term gains a beta-gradient response. -/
namespace S11VariableCoefficients.D4
noncomputable section
open S10Pilot
open scoped ContDiff

def lagrangian (beta : Point 4 → ℝ) (x : Point 4) (J : Jet 4) : ℝ :=
  S11D4Odd.lagrangian (beta x) J

def momentum (beta : Point 4 → ℝ) (x : Point 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) : ℝ :=
  deriv (fun s : ℝ => lagrangian beta x (J + s • basisJet j i)) 0

theorem momentum_identity (beta : Point 4 → ℝ) (x : Point 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) :
    momentum beta x J j i = S11D4Odd.momentum (beta x) J j i := rfl

theorem pointwise_variation (beta : Point 4 → ℝ) (x : Point 4) (J H : Jet 4) :
    HasDerivAt (fun s : ℝ => lagrangian beta x (J + s • H))
      (S11D4Odd.variationDensity (beta x) J H) 0 := S11D4Odd.lagrangian_variation _ _ _

def eulerLagrange (beta : Point 4 → ℝ) (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun j => -∑ i : Fin 5, coordDeriv i (fun y => momentum beta y (fieldJet u y) i j) x

theorem eulerLagrange_eq {beta : Point 4 → ℝ} {u : Point 4 → Vec 4}
    (hb : ContDiff ℝ ∞ beta) (hu : SmoothField u) (x : Point 4) (j : Fin 4) :
    eulerLagrange beta u x j = (1 / 2 : ℝ) *
      ∑ i : Fin 4, coordDeriv i.succ beta x * S11D4Odd.dualCurl (fieldJet u x) i j := by
  have hm : ∀ i, ContDiff ℝ ∞ (fun y => S11D4Odd.dualCurl (fieldJet u y) i j) :=
    fun i => S11D4Odd.smooth_dualCurl hu i j
  have he (i : Fin 4) :
      (fun y => -(beta y / 2) * S11D4Odd.dualCurl (fieldJet u y) i j) =
      (fun y => -(1 / 2 : ℝ) * (beta y * S11D4Odd.dualCurl (fieldJet u y) i j)) := by
    funext y
    ring
  unfold eulerLagrange
  rw [Fin.sum_univ_succ]
  simp only [momentum_identity, S11D4Odd.momentum_eq, Fin.cases_zero, Fin.cases_succ,
    S11D4Odd.partial_zero, zero_add]
  simp_rw [he, S11D4Odd.partial_const_mul (hb.mul (hm _)), S11D4Odd.partial_mul hb (hm _)]
  simp_rw [← Finset.mul_sum, Finset.sum_add_distrib, ← Finset.mul_sum,
    S11D4Odd.dualCurl_divergence_zero hu]
  ring

theorem constant_profile {u : Point 4 → Vec 4} (hu : SmoothField u)
    (beta : ℝ) (x : Point 4) : eulerLagrange (fun _ => beta) u x = 0 :=
  S11D4Odd.eulerLagrange_zero hu beta x

theorem weighted_odd_density {beta : Point 4 → ℝ} {u : Point 4 → Vec 4}
    (hb : ContDiff ℝ ∞ beta) (hu : SmoothField u) (x : Point 4) :
    lagrangian beta x (fieldJet u x) =
      -(1 / 2 : ℝ) * spatialDiv (fun y i => beta y * S11D4Odd.boundaryCurrent u y i) x +
        (1 / 2 : ℝ) * gradientPair beta (S11D4Odd.boundaryCurrent u) x := by
  have hk : ∀ i, ContDiff ℝ ∞ (fun y => S11D4Odd.boundaryCurrent u y i) := by
    have hU : ∀ j, ContDiff ℝ ∞ (fun y => u y j) := hu
    have hM : ∀ i j, ContDiff ℝ ∞ (fun y => S11D4Odd.dualCurl (fieldJet u y) i j) :=
      S11D4Odd.smooth_dualCurl hu
    intro i
    unfold S11D4Odd.boundaryCurrent S11D4Odd.currentJet
    fun_prop
  rw [weighted_divergence hb hk]
  have he : spatialDiv (S11D4Odd.boundaryCurrent u) x = S11D4Odd.orientationJet (fieldJet u x) :=
    S11D4Odd.boundary_identity hu x
  rw [he]
  unfold lagrangian S11D4Odd.lagrangian
  ring

def traction (beta : ℝ) (J : Jet 4) (n : Vec 4) : Vec 4 :=
  normalFlux n (fun i j => S11D4Odd.momentum beta J i.succ j)

theorem traction_eq (beta : ℝ) (J : Jet 4) (n : Vec 4) (j : Fin 4) :
    traction beta J n j = -(beta / 2) * ∑ i : Fin 4, n i * S11D4Odd.dualCurl J i j := by
  simp only [traction, normalFlux, S11D4Odd.momentum_eq, Fin.cases_succ, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

end
end S11VariableCoefficients.D4
