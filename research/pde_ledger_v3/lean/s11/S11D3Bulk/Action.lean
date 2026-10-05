import S10Pilot.Variation
import S11D3Invariants.Classification

/-! K1: the full classified D3 density, with constant real coefficients.
L = -Q/2. Spatial rows are derivative indices; no inertia is added. -/
namespace S11D3Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

def spatialGradient (J : Jet 3) : S11D3Invariants.Mat := fun j i => J j.succ i
def divergence (J : Jet 3) : ℝ := ∑ i : Fin 3, J i.succ i
def transposePair (J : Jet 3) : ℝ := ∑ i : Fin 3, ∑ j : Fin 3, J i.succ j * J j.succ i
def gradientSquare (J : Jet 3) : ℝ := ∑ i : Fin 3, ∑ j : Fin 3, J i.succ j ^ 2

def lagrangian (v : Vec 3) (J : Jet 3) : ℝ :=
  -(1 / 2 : ℝ) * (v 0 * divergence J ^ 2 + v 1 * transposePair J + v 2 * gradientSquare J)

theorem density_identity (v : Vec 3) (J : Jet 3) :
    lagrangian v J = -(1 / 2 : ℝ) * S11D3Invariants.invariantForm v (spatialGradient J) := by
  rw [S11D3Invariants.invariantForm_apply]
  simp only [lagrangian, divergence, transposePair, gradientSquare, spatialGradient,
    Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Fin.sum_univ_three]
  ring

theorem all_invariant_densities (Q : S11D3Invariants.Quad)
    (hQ : S11D3Invariants.SOInvariant Q) :
    ∃! v : Vec 3, Q = S11D3Invariants.invariantForm v :=
  S11D3Invariants.SO_unique Q hQ

def variationDensity (v : Vec 3) (J H : Jet 3) : ℝ :=
  -(v 0 * divergence J * divergence H +
    v 1 * (∑ i : Fin 3, ∑ j : Fin 3, J j.succ i * H i.succ j) +
    v 2 * (∑ i : Fin 3, ∑ j : Fin 3, J i.succ j * H i.succ j))

theorem lagrangian_increment (v : Vec 3) (s : ℝ) (J H : Jet 3) :
    lagrangian v (J + s • H) = lagrangian v J +
      s * variationDensity v J H + s ^ 2 * lagrangian v H := by
  simp only [lagrangian, variationDensity, divergence, transposePair, gradientSquare,
    Fin.sum_univ_three, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (v : Vec 3) (J H : Jet 3) :
    HasDerivAt (fun s : ℝ => lagrangian v (J + s • H)) (variationDensity v J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian v J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity v J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian v H))
  convert! h using 1
  · funext s
    exact lagrangian_increment v s J H
  · simp

def momentum (v : Vec 3) (J : Jet 3) (j : Fin 4) (i : Fin 3) : ℝ :=
  deriv (fun s : ℝ => lagrangian v (J + s • basisJet j i)) 0

theorem momentum_eq (v : Vec 3) (J : Jet 3) (j : Fin 4) (i : Fin 3) :
    momentum v J j i = Fin.cases 0
      (fun r => -(v 0 * (if r = i then divergence J else 0) +
        v 1 * J i.succ r + v 2 * J r.succ i)) j := by
  unfold momentum
  rw [(lagrangian_variation v J (basisJet j i)).deriv]
  cases j using Fin.cases
  · fin_cases i <;> simp [variationDensity, divergence, basisJet_apply, Fin.sum_univ_three]
  · rename_i r
    simp only [Fin.cases_succ]
    fin_cases r <;> fin_cases i <;>
      simp [variationDensity, divergence, basisJet_apply, Fin.sum_univ_three]

theorem momentum_smul (v : Vec 3) (s : ℝ) (J : Jet 3) (j : Fin 4) (i : Fin 3) :
    momentum v (s • J) j i = s * momentum v J j i := by
  rw [momentum_eq, momentum_eq]
  cases j using Fin.cases
  · simp
  · rename_i r
    by_cases h : r = i <;> simp [h, divergence, Fin.sum_univ_three] <;> ring

def eulerLagrange (v : Vec 3) (u : Point 3 → Vec 3) (x : Point 3) : Vec 3 :=
  fun i => -∑ j : Fin 4, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x

theorem smooth_momentum {u : Point 3 → Vec 3} (hu : SmoothField u)
    (v : Vec 3) (j : Fin 4) (i : Fin 3) :
    ContDiff ℝ ∞ (fun x => momentum v (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  cases j using Fin.cases
  · simp only [momentum_eq, Fin.cases_zero]
    fun_prop
  · rename_i r
    by_cases h : r = i <;> simp only [momentum_eq, Fin.cases_succ, h, if_true, if_false, divergence]
    all_goals fun_prop

theorem smooth_eulerLagrange {u : Point 3 → Vec 3} (hu : SmoothField u)
    (v : Vec 3) (i : Fin 3) : ContDiff ℝ ∞ (fun x => eulerLagrange v u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞ (coordDeriv j (fun y => momentum v (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu v j i) j
  unfold eulerLagrange
  fun_prop

end
end S11D3Bulk
