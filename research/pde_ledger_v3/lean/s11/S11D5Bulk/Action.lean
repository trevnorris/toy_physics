import S10Pilot.Variation
import S11D5Invariants.Classification

/-! D5B.1: the full classified D5 density, with constant real coefficients.
L = -Q/2. Spatial rows are derivative indices; no inertia is added. -/
namespace S11D5Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

/-- The three coefficients are independent of the five field components. -/
abbrev Coeff := Vec 3

def spatialGradient (J : Jet 5) : S11D5Invariants.Mat := fun j i => J j.succ i
def divergence (J : Jet 5) : ℝ := ∑ i : Fin 5, J i.succ i
def transposePair (J : Jet 5) : ℝ := ∑ i : Fin 5, ∑ j : Fin 5, J i.succ j * J j.succ i
def gradientSquare (J : Jet 5) : ℝ := ∑ i : Fin 5, ∑ j : Fin 5, J i.succ j ^ 2

def lagrangian (v : Coeff) (J : Jet 5) : ℝ :=
  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2 + v 1 * transposePair J + v 2 * gradientSquare J)

theorem density_identity (v : Coeff) (J : Jet 5) :
    lagrangian v J = -(1 / 2 : ℝ) * S11D5Invariants.invariantForm v (spatialGradient J) := by
  rw [S11D5Invariants.invariantForm_apply]
  simp only [lagrangian, divergence, transposePair, gradientSquare, spatialGradient,
    Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
    S11D5Invariants.sum_five]
  ring

theorem all_invariant_densities (Q : S11D5Invariants.Quad)
    (hQ : S11D5Invariants.SOInvariant Q) :
    ∃! v : Coeff, Q = S11D5Invariants.invariantForm v :=
  S11D5Invariants.SO_unique Q hQ

def variationDensity (v : Coeff) (J H : Jet 5) : ℝ :=
  -(v 0 * divergence J * divergence H +
    v 1 * (∑ i : Fin 5, ∑ j : Fin 5, J j.succ i * H i.succ j) +
    v 2 * (∑ i : Fin 5, ∑ j : Fin 5, J i.succ j * H i.succ j))

theorem lagrangian_increment (v : Coeff) (s : ℝ) (J H : Jet 5) :
    lagrangian v (J + s • H) = lagrangian v J +
      s * variationDensity v J H + s ^ 2 * lagrangian v H := by
  simp only [lagrangian, variationDensity, divergence, transposePair, gradientSquare,
    S11D5Invariants.sum_five, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (v : Coeff) (J H : Jet 5) :
    HasDerivAt (fun s : ℝ => lagrangian v (J + s • H)) (variationDensity v J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian v J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity v J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian v H))
  convert! h using 1
  · funext s
    exact lagrangian_increment v s J H
  · simp

def momentum (v : Coeff) (J : Jet 5) (j : Fin 6) (i : Fin 5) : ℝ :=
  deriv (fun s : ℝ => lagrangian v (J + s • basisJet j i)) 0

theorem momentum_eq (v : Coeff) (J : Jet 5) (j : Fin 6) (i : Fin 5) :
    momentum v J j i = Fin.cases 0
      (fun r => -(v 0 * (if r = i then divergence J else 0) +
        v 1 * J i.succ r + v 2 * J r.succ i)) j := by
  unfold momentum
  rw [(lagrangian_variation v J (basisJet j i)).deriv]
  cases j using Fin.cases
  · fin_cases i <;> simp [variationDensity, divergence, basisJet_apply, S11D5Invariants.sum_five]
  · rename_i r
    simp only [Fin.cases_succ]
    fin_cases r <;> fin_cases i <;>
      simp [variationDensity, divergence, basisJet_apply, S11D5Invariants.sum_five]

theorem momentum_time (v : Coeff) (J : Jet 5) (i : Fin 5) : momentum v J 0 i = 0 := by
  simp [momentum_eq]

theorem momentum_smul (v : Coeff) (s : ℝ) (J : Jet 5) (j : Fin 6) (i : Fin 5) :
    momentum v (s • J) j i = s * momentum v J j i := by
  rw [momentum_eq, momentum_eq]
  cases j using Fin.cases
  · simp
  · rename_i r
    by_cases h : r = i <;> simp [h, divergence, S11D5Invariants.sum_five] <;> ring

def eulerLagrange (v : Coeff) (u : Point 5 → Vec 5) (x : Point 5) : Vec 5 :=
  fun i => -∑ j : Fin 6, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x

theorem smooth_momentum {u : Point 5 → Vec 5} (hu : SmoothField u)
    (v : Coeff) (j : Fin 6) (i : Fin 5) :
    ContDiff ℝ ∞ (fun x => momentum v (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  cases j using Fin.cases
  · simp only [momentum_eq, Fin.cases_zero]
    fun_prop
  · rename_i r
    by_cases h : r = i <;> simp only [momentum_eq, Fin.cases_succ, h, if_true, if_false, divergence]
    all_goals fun_prop

theorem smooth_eulerLagrange {u : Point 5 → Vec 5} (hu : SmoothField u)
    (v : Coeff) (i : Fin 5) : ContDiff ℝ ∞ (fun x => eulerLagrange v u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞ (coordDeriv j (fun y => momentum v (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu v j i) j
  unfold eulerLagrange
  fun_prop

end
end S11D5Bulk
