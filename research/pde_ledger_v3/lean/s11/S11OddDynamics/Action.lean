import S10Pilot.Variation
import S11Invariants.Controls

/-! E1 and the local E2/E3 identities. Spatial rows are derivative indices;
coordinates are (t,x1,x2). This is only XFORM_EXTRA minus MAIN, with constant beta. -/
namespace S11OddDynamics
noncomputable section
open S10Pilot
open scoped ContDiff

def spatialGradient (J : Jet 2) : S11Invariants.Mat := fun j i => J j.succ i

def divergence (J : Jet 2) : ℝ := J (0 : Fin 2).succ 0 + J (1 : Fin 2).succ 1
def skew (J : Jet 2) : ℝ := J (0 : Fin 2).succ 1 - J (1 : Fin 2).succ 0

def lagrangian (beta : ℝ) (J : Jet 2) : ℝ :=
  -beta / 2 * divergence J * skew J

theorem density_identity (beta : ℝ) (J : Jet 2) :
    lagrangian beta J = -beta / 2 * S11Invariants.oddPairing (spatialGradient J) := by
  simp only [S11Invariants.oddPairing, S11Invariants.invariantForm_apply,
    S11Invariants.coordinates, Matrix.cons_val, spatialGradient]
  simp [lagrangian, divergence, skew]
  ring

def variationDensity (beta : ℝ) (J H : Jet 2) : ℝ :=
  -beta / 2 * (divergence H * skew J + divergence J * skew H)

theorem lagrangian_increment (beta s : ℝ) (J H : Jet 2) :
    lagrangian beta (J + s • H) = lagrangian beta J +
      s * variationDensity beta J H + s ^ 2 * lagrangian beta H := by
  simp only [lagrangian, variationDensity, divergence, skew, Pi.add_apply, Pi.smul_apply,
    smul_eq_mul]
  ring

theorem lagrangian_smul (beta s : ℝ) (J : Jet 2) :
    lagrangian beta (s • J) = s ^ 2 * lagrangian beta J := by
  simp only [lagrangian, divergence, skew, Pi.smul_apply, smul_eq_mul]
  ring

theorem lagrangian_variation (beta : ℝ) (J H : Jet 2) :
    HasDerivAt (fun s : ℝ => lagrangian beta (J + s • H))
      (variationDensity beta J H) 0 := by
  have h := ((hasDerivAt_const (0 : ℝ) (lagrangian beta J)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (variationDensity beta J H))).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).mul_const (lagrangian beta H))
  convert! h using 1
  · funext s
    exact lagrangian_increment beta s J H
  · simp

/-- An actual derivative of the supplied density, not an asserted constitutive law. -/
def momentum (beta : ℝ) (J : Jet 2) (j : Fin 3) (i : Fin 2) : ℝ :=
  deriv (fun s : ℝ => lagrangian beta (J + s • basisJet j i)) 0

theorem momentum_eq (beta : ℝ) (J : Jet 2) (j : Fin 3) (i : Fin 2) :
    momentum beta J j i =
      !![0, 0; -beta / 2 * skew J, -beta / 2 * divergence J;
        beta / 2 * divergence J, -beta / 2 * skew J] j i := by
  unfold momentum
  rw [(lagrangian_variation beta J (basisJet j i)).deriv]
  fin_cases j <;> fin_cases i <;>
    simp [variationDensity, divergence, skew, basisJet_apply]
  all_goals ring

theorem momentum_smul (beta s : ℝ) (J : Jet 2) (j : Fin 3) (i : Fin 2) :
    momentum beta (s • J) j i = s * momentum beta J j i := by
  rw [momentum_eq, momentum_eq]
  fin_cases j <;> fin_cases i <;> simp [divergence, skew]
  all_goals ring

/-- Variational sign dL/du - sum_j partial_j(dL/dJ_ji). -/
def eulerLagrange (beta : ℝ) (u : Point 2 → Vec 2) (x : Point 2) : Vec 2 :=
  fun i => -∑ j : Fin 3, coordDeriv j (fun y => momentum beta (fieldJet u y) j i) x

theorem smooth_momentum {u : Point 2 → Vec 2} (hu : SmoothField u)
    (beta : ℝ) (j : Fin 3) (i : Fin 2) :
    ContDiff ℝ ∞ (fun x => momentum beta (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  fin_cases j <;> fin_cases i <;> simp [momentum_eq, divergence, skew]
  all_goals fun_prop

theorem smooth_eulerLagrange {u : Point 2 → Vec 2} (hu : SmoothField u) (beta : ℝ)
    (i : Fin 2) : ContDiff ℝ ∞ (fun x => eulerLagrange beta u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞ (coordDeriv j (fun y => momentum beta (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu beta j i) j
  unfold eulerLagrange
  fun_prop

theorem coordDeriv_const_mul {f : Point 2 → ℝ} (hf : Differentiable ℝ f)
    (c : ℝ) (j : Fin 3) (x : Point 2) :
    coordDeriv j (fun y => c * f y) x = c * coordDeriv j f x := by
  rw [coordDeriv_eq_fderiv (hf.const_mul c), coordDeriv_eq_fderiv hf]
  simp [fderiv_const_mul (hf x)]

theorem coordDeriv_zero (j : Fin 3) (x : Point 2) :
    coordDeriv j (fun _ => (0 : ℝ)) x = 0 := by simp [coordDeriv]

/-- Compact local formula; no commutation of mixed derivatives is required. -/
theorem eulerLagrange_eq {u : Point 2 → Vec 2} (hu : SmoothField u) (beta : ℝ)
    (x : Point 2) : eulerLagrange beta u x =
    ![beta / 2 * (coordDeriv 1 (fun y => skew (fieldJet u y)) x -
        coordDeriv 2 (fun y => divergence (fieldJet u y)) x),
      beta / 2 * (coordDeriv 1 (fun y => divergence (fieldJet u y)) x +
        coordDeriv 2 (fun y => skew (fieldJet u y)) x)] := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun y => fieldJet u y a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  have hd : ContDiff ℝ ∞ (fun y => divergence (fieldJet u y)) := by unfold divergence; fun_prop
  have hs : ContDiff ℝ ∞ (fun y => skew (fieldJet u y)) := by unfold skew; fun_prop
  ext i
  fin_cases i <;> simp [eulerLagrange, momentum_eq, Fin.sum_univ_succ,
    coordDeriv_const_mul (hs.differentiable (by simp)),
    coordDeriv_const_mul (hd.differentiable (by simp)), coordDeriv_zero] <;> ring

/-- The positively oriented perpendicular direction for row-major curl convention. -/
def turn (k : Vec 2) : Vec 2 := ![-k 1, k 0]

def modalAction (beta : ℝ) (k a : Vec 2) : ℝ :=
  lagrangian beta (modeJet 0 k a)

/-- Identified with both the action derivative and the local PDE below. -/
def modalOperator (beta : ℝ) (k a : Vec 2) : Vec 2 :=
  (-beta / 2) • ((dot (turn k) a) • k + (dot k a) • turn k)

theorem modalAction_eq (beta : ℝ) (k a : Vec 2) :
    modalAction beta k a = -beta / 2 * dot k a * dot (turn k) a := by
  change -beta / 2 * (k 0 * a 0 + k 1 * a 1) * (k 0 * a 1 - k 1 * a 0) = _
  simp only [dot, turn, Fin.sum_univ_two, Matrix.cons_val]
  ring

theorem modal_variation (beta : ℝ) (k a b : Vec 2) :
    HasDerivAt (fun s : ℝ => modalAction beta k (a + s • b))
      (dot (modalOperator beta k a) b) 0 := by
  simp_rw [modalAction_eq, dot_add_right, dot_smul_right]
  have h := (((hasDerivAt_const (0 : ℝ) (dot k a)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (dot k b))).const_mul (-beta / 2)).mul
    ((hasDerivAt_const (0 : ℝ) (dot (turn k) a)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (dot (turn k) b)))
  convert! h using 1
  simp [modalOperator, dot_smul_left, dot_add_left]
  ring

theorem momentum_contraction (beta omega : ℝ) (k a : Vec 2) (i : Fin 2) :
    (∑ j : Fin 3, waveCovector omega k j * momentum beta (modeJet (-omega) k a) j i) =
      modalOperator beta k a i := by
  simp only [Fin.sum_univ_succ, waveCovector, Fin.cases_zero, Fin.cases_succ,
    momentum_eq, Fin.sum_univ_zero]
  fin_cases i
  · change -omega * 0 + (k 0 * (-beta / 2 * (k 0 * a 1 - k 1 * a 0)) +
        (k 1 * (beta / 2 * (k 0 * a 0 + k 1 * a 1)) + 0)) = _
    simp [modalOperator, dot, turn, Fin.sum_univ_two, Matrix.vecHead, Matrix.vecTail]
    ring
  · change -omega * 0 + (k 0 * (-beta / 2 * (k 0 * a 0 + k 1 * a 1)) +
        (k 1 * (-beta / 2 * (k 0 * a 1 - k 1 * a 0)) + 0)) = _
    simp [modalOperator, dot, turn, Fin.sum_univ_two, Matrix.vecHead, Matrix.vecTail]
    ring

theorem eulerLagrange_planeWave (beta omega : ℝ) (k a : Vec 2) (x : Point 2) :
    eulerLagrange beta (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator beta k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have hrearrange : ∀ j : Fin 3,
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum beta (modeJet (-omega) k a) j i) =
      (fun y => -momentum beta (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [hrearrange, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction beta omega]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

end
end S11OddDynamics
