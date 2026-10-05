import S9Pilot.Action

/-!
The position-space route uses actual scalar coordinate derivatives. The
Euler–Lagrange differential expression is defined by differentiating the supplied
density with respect to its derivative coordinates. The Variation module proves
equivalence with stationarity of the integrated action change under smooth,
compactly supported variations. No boundary or decay assumptions are needed to
evaluate this local expression on a smooth plane wave.
-/

namespace S9Pilot
noncomputable section

def axis (j : Fin 4) : Point := Pi.single j 1
def basisJet (j : Fin 4) (i : Fin 3) : Jet := Pi.single j (Pi.single i 1)
def phase (q x : Point) : ℝ := q 0 * x 0 + q 1 * x 1 + q 2 * x 2 + q 3 * x 3
def coordDeriv (j : Fin 4) (f : Point → ℝ) (x : Point) : ℝ :=
  deriv (fun s : ℝ => f (x + s • axis j)) 0
def fieldJet (u : Point → Vec) (x : Point) : Jet :=
  fun j i => coordDeriv j (fun x => u x i) x
def waveCovector (omega : ℝ) (k : Vec) : Point := ![-omega, k 0, k 1, k 2]
def planeWave (omega : ℝ) (k a : Vec) (x : Point) : Vec :=
  fun i => a i * Real.cos (phase (waveCovector omega k) x)

theorem phase_shift (q x : Point) (j : Fin 4) (s : ℝ) :
    phase q (x + s • axis j) = phase q x + s * q j := by
  fin_cases j <;> simp [phase, axis] <;> ring

theorem partial_const_cos (q x : Point) (j : Fin 4) (c : ℝ) :
    coordDeriv j (fun y => c * Real.cos (phase q y)) x =
      -c * Real.sin (phase q x) * q j := by
  unfold coordDeriv
  simp_rw [phase_shift]
  have h := (((hasDerivAt_const (0 : ℝ) (phase q x)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (q j))).cos).const_mul c
  convert! h.deriv using 1
  simp
  ring

theorem partial_const_sin (q x : Point) (j : Fin 4) (c : ℝ) :
    coordDeriv j (fun y => c * Real.sin (phase q y)) x =
      c * Real.cos (phase q x) * q j := by
  unfold coordDeriv
  simp_rw [phase_shift]
  have h := (((hasDerivAt_const (0 : ℝ) (phase q x)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (q j))).sin).const_mul c
  convert! h.deriv using 1
  simp
  ring

theorem fieldJet_planeWave (omega : ℝ) (k a : Vec) (x : Point) :
    fieldJet (planeWave omega k a) x =
      (-Real.sin (phase (waveCovector omega k) x)) • modeJet (-omega) k a := by
  ext j i
  simp only [fieldJet, planeWave, partial_const_cos]
  fin_cases j <;> simp [waveCovector, modeJet] <;> ring

/-- Exact pointwise density identity; no infinite-volume action integral is taken. -/
theorem planeWave_lagrangian (rho mu omega : ℝ) (k a : Vec) (x : Point) :
    lagrangian rho mu (fieldJet (planeWave omega k a) x) =
      Real.sin (phase (waveCovector omega k) x) ^ 2 * modalAction rho mu omega k a := by
  rw [fieldJet_planeWave]
  simp [lagrangian, modalAction, modeJet, normSq, dot, jetCurl]
  ring

/-- Conjugate derivative obtained from the density by an actual derivative. -/
def momentum (rho mu : ℝ) (J : Jet) (j : Fin 4) (i : Fin 3) : ℝ :=
  deriv (fun s : ℝ => lagrangian rho mu (J + s • basisJet j i)) 0

theorem momentum_eq (rho mu : ℝ) (J : Jet) (j : Fin 4) (i : Fin 3) :
    momentum rho mu J j i =
      rho * dot (J 0) (basisJet j i 0) - mu * dot (jetCurl J) (jetCurl (basisJet j i)) :=
  (lagrangian_variation rho mu J (basisJet j i)).deriv

theorem momentum_smul (rho mu s : ℝ) (J : Jet) (j : Fin 4) (i : Fin 3) :
    momentum rho mu (s • J) j i = s * momentum rho mu J j i := by
  simp [momentum_eq, dot, jetCurl]
  ring

/-- First-order, field-independent density: dL/du - sum_j d_j(dL/dJ_ji). -/
def eulerLagrange (rho mu : ℝ) (u : Point → Vec) (x : Point) : Vec :=
  fun i => -∑ j : Fin 4, coordDeriv j (fun y => momentum rho mu (fieldJet u y) j i) x

theorem momentum_contraction (rho mu omega : ℝ) (k a : Vec) (i : Fin 3) :
    (∑ j : Fin 4, waveCovector omega k j * momentum rho mu (modeJet (-omega) k a) j i) =
      modalOperator rho mu omega k a i := by
  fin_cases i <;>
    simp [Fin.sum_univ_succ, waveCovector, momentum_eq, modeJet, basisJet, jetCurl,
      dot, normSq, modalOperator] <;> ring

/-- The independently differentiated position-space expression yields the same matrix. -/
theorem eulerLagrange_planeWave (rho mu omega : ℝ) (k a : Vec) (x : Point) :
    eulerLagrange rho mu (planeWave omega k a) x =
      Real.cos (phase (waveCovector omega k) x) • modalOperator rho mu omega k a := by
  ext i
  unfold eulerLagrange
  simp_rw [fieldJet_planeWave, momentum_smul]
  have hrearrange : ∀ j : Fin 4,
      (fun y => -Real.sin (phase (waveCovector omega k) y) *
        momentum rho mu (modeJet (-omega) k a) j i) =
      (fun y => -momentum rho mu (modeJet (-omega) k a) j i *
        Real.sin (phase (waveCovector omega k) y)) := by
    intro j
    funext y
    ring
  simp_rw [hrearrange, partial_const_sin]
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [← momentum_contraction]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- A global plane wave solves the local EL expression iff its amplitude is stationary. -/
theorem planeWave_solves_iff (rho mu omega : ℝ) (k a : Vec) :
    (∀ x, eulerLagrange rho mu (planeWave omega k a) x = 0) ↔
      ModalStationary rho mu omega k a := by
  rw [modal_stationary_iff]
  constructor
  · intro h
    simpa [eulerLagrange_planeWave, phase] using h 0
  · intro h x
    simp [eulerLagrange_planeWave, h]

end
end S9Pilot
