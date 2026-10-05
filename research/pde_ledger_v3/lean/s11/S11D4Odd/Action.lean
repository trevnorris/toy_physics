import S10Pilot.Analytic
import S11D4Invariants.Classification

/-! D4B.1: the supplied constant-coefficient odd density and its actual derivative.
The spatial row is the derivative index. No even density or inertia is added. -/
namespace S11D4Odd
noncomputable section
open S10Pilot
open scoped ContDiff

def spatialGradient (J : Jet 4) : S11D4Invariants.Mat := fun i j => J i.succ j
def orientationJet (J : Jet 4) : ℝ := S11D4Invariants.orientation (spatialGradient J)

theorem orientationJet_eq (J : Jet 4) :
    orientationJet J = antisym J 0 1 * antisym J 2 3 -
      antisym J 0 2 * antisym J 1 3 + antisym J 0 3 * antisym J 1 2 := rfl

def dualCurl (J : Jet 4) : S11D4Invariants.Mat := Matrix.of
  ![![0, antisym J 2 3, -antisym J 1 3, antisym J 1 2],
    ![-antisym J 2 3, 0, antisym J 0 3, -antisym J 0 2],
    ![antisym J 1 3, -antisym J 0 3, 0, antisym J 0 1],
    ![-antisym J 1 2, antisym J 0 2, -antisym J 0 1, 0]]

def lagrangian (beta : ℝ) (J : Jet 4) : ℝ := -(beta / 2) * orientationJet J

theorem density_identity (beta : ℝ) (J : Jet 4) :
    lagrangian beta J = -(1 / 2 : ℝ) *
      S11D4Invariants.invariantForm ![0,0,0,beta] (spatialGradient J) := by
  rw [S11D4Invariants.invariantForm_apply]
  simp only [Matrix.cons_val, lagrangian, orientationJet, zero_mul, zero_add]
  ring

theorem all_odd_densities (Q : S11D4Invariants.Quad) :
    S11D4Invariants.SOInvariant Q ∧ S11D4Invariants.ReflectionOdd Q ↔
      ∃ beta : ℝ, Q = S11D4Invariants.invariantForm ![0,0,0,beta] :=
  S11D4Invariants.odd_classification Q

def variationDensity (beta : ℝ) (J H : Jet 4) : ℝ := -(beta / 2) *
  (antisym H 0 1 * antisym J 2 3 + antisym J 0 1 * antisym H 2 3 -
   antisym H 0 2 * antisym J 1 3 - antisym J 0 2 * antisym H 1 3 +
   antisym H 0 3 * antisym J 1 2 + antisym J 0 3 * antisym H 1 2)

theorem lagrangian_increment (beta s : ℝ) (J H : Jet 4) :
    lagrangian beta (J + s • H) = lagrangian beta J +
      s * variationDensity beta J H + s ^ 2 * lagrangian beta H := by
  simp only [lagrangian, orientationJet_eq, variationDensity, antisym,
    Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (beta : ℝ) (J H : Jet 4) :
    HasDerivAt (fun s : ℝ => lagrangian beta (J + s • H))
      (variationDensity beta J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian beta J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity beta J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian beta H))
  convert! h using 1
  · funext s
    exact lagrangian_increment beta s J H
  · simp

def momentum (beta : ℝ) (J : Jet 4) (j : Fin 5) (i : Fin 4) : ℝ :=
  deriv (fun s : ℝ => lagrangian beta (J + s • basisJet j i)) 0

theorem momentum_eq (beta : ℝ) (J : Jet 4) (j : Fin 5) (i : Fin 4) :
    momentum beta J j i = Fin.cases 0 (fun r => -(beta / 2) * dualCurl J r i) j := by
  unfold momentum
  rw [(lagrangian_variation beta J (basisJet j i)).deriv]
  cases j using Fin.cases
  · fin_cases i <;> simp [variationDensity, antisym, basisJet_apply]
  · rename_i r
    simp only [Fin.cases_succ]
    fin_cases r <;> fin_cases i <;>
      simp [variationDensity, dualCurl, antisym, basisJet_apply,
        Matrix.cons_val_two, Matrix.cons_val_three]

theorem dualCurl_contraction (J : Jet 4) :
    (∑ i : Fin 4, ∑ j : Fin 4, J i.succ j * dualCurl J i j) = 2 * orientationJet J := by
  simp only [Fin.sum_univ_four, dualCurl, Matrix.of_apply, Matrix.cons_val,
    orientationJet_eq, antisym]
  ring

def eulerLagrange (beta : ℝ) (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => -∑ j : Fin 5, coordDeriv j (fun y => momentum beta (fieldJet u y) j i) x

theorem smooth_dualCurl {u : Point 4 → Vec 4} (hu : SmoothField u) (i j : Fin 4) :
    ContDiff ℝ ∞ (fun x => dualCurl (fieldJet u x) i j) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  fin_cases i <;> fin_cases j <;>
    simp [dualCurl, antisym, Matrix.cons_val_two, Matrix.cons_val_three]
  all_goals fun_prop

theorem smooth_momentum {u : Point 4 → Vec 4} (hu : SmoothField u)
    (beta : ℝ) (j : Fin 5) (i : Fin 4) :
    ContDiff ℝ ∞ (fun x => momentum beta (fieldJet u x) j i) := by
  cases j using Fin.cases
  · simp only [momentum_eq, Fin.cases_zero]
    fun_prop
  · rename_i r
    simp only [momentum_eq, Fin.cases_succ]
    exact contDiff_const.mul (smooth_dualCurl hu r i)

theorem smooth_eulerLagrange {u : Point 4 → Vec 4} (hu : SmoothField u)
    (beta : ℝ) (i : Fin 4) : ContDiff ℝ ∞ (fun x => eulerLagrange beta u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞
      (coordDeriv j (fun y => momentum beta (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu beta j i) j
  unfold eulerLagrange
  fun_prop

end
end S11D4Odd
