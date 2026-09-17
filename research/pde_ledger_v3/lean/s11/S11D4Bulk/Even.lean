import S10Pilot.Variation
import S11D4Invariants.Classification

/-! D4C.1: the even part of the full D4 density, with constant real coefficients.
L = -Q/2. Spatial rows are derivative indices; no inertia is added. -/
namespace S11D4Bulk.Even
noncomputable section
open S10Pilot
open scoped ContDiff

def spatialGradient (J : Jet 4) : S11D4Invariants.Mat := fun j i => J j.succ i
def divergence (J : Jet 4) : ℝ := ∑ i : Fin 4, J i.succ i
def transposePair (J : Jet 4) : ℝ := ∑ i : Fin 4, ∑ j : Fin 4, J i.succ j * J j.succ i
def gradientSquare (J : Jet 4) : ℝ := ∑ i : Fin 4, ∑ j : Fin 4, J i.succ j ^ 2

def lagrangian (v : Vec 4) (J : Jet 4) : ℝ :=
  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2 + v 1 * transposePair J + v 2 * gradientSquare J)

def variationDensity (v : Vec 4) (J H : Jet 4) : ℝ :=
  -(v 0 * divergence J * divergence H +
    v 1 * (∑ i : Fin 4, ∑ j : Fin 4, J j.succ i * H i.succ j) +
    v 2 * (∑ i : Fin 4, ∑ j : Fin 4, J i.succ j * H i.succ j))

theorem lagrangian_increment (v : Vec 4) (s : ℝ) (J H : Jet 4) :
    lagrangian v (J + s • H) = lagrangian v J +
      s * variationDensity v J H + s ^ 2 * lagrangian v H := by
  simp only [lagrangian, variationDensity, divergence, transposePair, gradientSquare,
    Fin.sum_univ_four, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (v : Vec 4) (J H : Jet 4) :
    HasDerivAt (fun s : ℝ => lagrangian v (J + s • H)) (variationDensity v J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian v J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity v J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian v H))
  convert! h using 1
  · funext s
    exact lagrangian_increment v s J H
  · simp

def momentum (v : Vec 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) : ℝ :=
  deriv (fun s : ℝ => lagrangian v (J + s • basisJet j i)) 0

theorem momentum_eq (v : Vec 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) :
    momentum v J j i = Fin.cases 0
      (fun r => -(v 0 * (if r = i then divergence J else 0) +
        v 1 * J i.succ r + v 2 * J r.succ i)) j := by
  unfold momentum
  rw [(lagrangian_variation v J (basisJet j i)).deriv]
  cases j using Fin.cases
  · fin_cases i <;> simp [variationDensity, divergence, basisJet_apply, Fin.sum_univ_four]
  · rename_i r
    simp only [Fin.cases_succ]
    fin_cases r <;> fin_cases i <;>
      simp [variationDensity, divergence, basisJet_apply, Fin.sum_univ_four]

theorem momentum_smul (v : Vec 4) (s : ℝ) (J : Jet 4) (j : Fin 5) (i : Fin 4) :
    momentum v (s • J) j i = s * momentum v J j i := by
  rw [momentum_eq, momentum_eq]
  cases j using Fin.cases
  · simp
  · rename_i r
    by_cases h : r = i <;> simp [h, divergence, Fin.sum_univ_four] <;> ring

def eulerLagrange (v : Vec 4) (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => -∑ j : Fin 5, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x

theorem smooth_momentum {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (j : Fin 5) (i : Fin 4) :
    ContDiff ℝ ∞ (fun x => momentum v (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  cases j using Fin.cases
  · simp only [momentum_eq, Fin.cases_zero]
    fun_prop
  · rename_i r
    by_cases h : r = i <;> simp only [momentum_eq, Fin.cases_succ, h, if_true, if_false, divergence]
    all_goals fun_prop

theorem smooth_eulerLagrange {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (i : Fin 4) : ContDiff ℝ ∞ (fun x => eulerLagrange v u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞ (coordDeriv j (fun y => momentum v (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu v j i) j
  unfold eulerLagrange
  fun_prop

end
end S11D4Bulk.Even
