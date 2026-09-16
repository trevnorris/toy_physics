import S9Pilot.PlaneWave
import S9Pilot.Spectrum

/-!
Narrow scalar-phase obstruction. The supplied Madelung velocity law is
`v = (hbar / mass) grad(theta)`. We differentiate both real plane-wave phase
quadratures about a constant phase, and prove that their complete velocity
amplitude family is longitudinal. This is a kinematic statement in a smooth,
single-valued phase chart (physically: nonzero density and no vortices), not a
derivation of GNLS dynamics, a branch count, or a general no-photon theorem.
-/

namespace S9Pilot.Madelung
noncomputable section

def phasePerturbation (theta0 epsilon omega A B : ℝ) (k : Vec) (x : Point) : ℝ :=
  theta0 + epsilon * (A * Real.cos (phase (waveCovector omega k) x) +
    B * Real.sin (phase (waveCovector omega k) x))

/-- The actual spatial derivative in the supplied velocity law. -/
def velocity (hbar mass : ℝ) (theta : Point → ℝ) (x : Point) : Vec :=
  fun i => (hbar / mass) * coordDeriv i.succ theta x

theorem phasePerturbation_coordDeriv (theta0 epsilon omega A B : ℝ)
    (k : Vec) (x : Point) (j : Fin 4) :
    coordDeriv j (phasePerturbation theta0 epsilon omega A B k) x =
      epsilon * (-A * Real.sin (phase (waveCovector omega k) x) +
        B * Real.cos (phase (waveCovector omega k) x)) * waveCovector omega k j := by
  unfold coordDeriv phasePerturbation
  simp_rw [phase_shift]
  have hp := (hasDerivAt_const (0 : ℝ) (phase (waveCovector omega k) x)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (waveCovector omega k j))
  have h := (hasDerivAt_const (0 : ℝ) theta0).add
    (((hp.cos.const_mul A).add (hp.sin.const_mul B)).const_mul epsilon)
  convert! h.deriv using 1
  simp
  ring

theorem velocity_phasePerturbation (hbar mass theta0 epsilon omega A B : ℝ)
    (k : Vec) (x : Point) :
    velocity hbar mass (phasePerturbation theta0 epsilon omega A B k) x =
      epsilon • (((hbar / mass) *
        (-A * Real.sin (phase (waveCovector omega k) x) +
          B * Real.cos (phase (waveCovector omega k) x))) • k) := by
  ext i
  simp only [velocity, phasePerturbation_coordDeriv, Pi.smul_apply, smul_eq_mul]
  have hi : waveCovector omega k i.succ = k i := by
    fin_cases i <;> rfl
  rw [hi]
  ring

/-- First variation with respect to perturbation size, not an assumed amplitude. -/
def linearVelocity (hbar mass theta0 omega A B : ℝ) (k : Vec) (x : Point) : Vec :=
  fun i => deriv (fun epsilon : ℝ =>
    velocity hbar mass (phasePerturbation theta0 epsilon omega A B k) x i) 0

def cosAmplitude (hbar mass B : ℝ) (k : Vec) : Vec := ((hbar / mass) * B) • k
def sinAmplitude (hbar mass A : ℝ) (k : Vec) : Vec := (-(hbar / mass) * A) • k

theorem linearVelocity_eq (hbar mass theta0 omega A B : ℝ) (k : Vec) (x : Point) :
    linearVelocity hbar mass theta0 omega A B k x =
      Real.cos (phase (waveCovector omega k) x) • cosAmplitude hbar mass B k +
      Real.sin (phase (waveCovector omega k) x) • sinAmplitude hbar mass A k := by
  ext i
  simp only [linearVelocity, velocity_phasePerturbation, Pi.smul_apply, smul_eq_mul]
  have h := (hasDerivAt_id (0 : ℝ)).mul_const
    ((hbar / mass * (-A * Real.sin (phase (waveCovector omega k) x) +
      B * Real.cos (phase (waveCovector omega k) x))) * k i)
  convert! h.deriv using 1
  simp [cosAmplitude, sinAmplitude]
  ring

theorem cosAmplitude_mem (hbar mass B : ℝ) (k : Vec) :
    cosAmplitude hbar mass B k ∈ longitudinalSpace k :=
  Submodule.mem_span_singleton.mpr ⟨(hbar / mass) * B, rfl⟩

theorem sinAmplitude_mem (hbar mass A : ℝ) (k : Vec) :
    sinAmplitude hbar mass A k ∈ longitudinalSpace k :=
  Submodule.mem_span_singleton.mpr ⟨-(hbar / mass) * A, rfl⟩

/-- All longitudinal amplitudes occur; the admissible family is not empty. -/
theorem cosAmplitude_range {hbar mass : ℝ} (hhbar : hbar ≠ 0) (hmass : mass ≠ 0)
    (k a : Vec) : (∃ B, cosAmplitude hbar mass B k = a) ↔ a ∈ longitudinalSpace k := by
  constructor
  · rintro ⟨B, rfl⟩
    exact cosAmplitude_mem _ _ _ _
  · intro ha
    obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp ha
    refine ⟨s / (hbar / mass), ?_⟩
    unfold cosAmplitude
    rw [mul_div_cancel₀ _ (div_ne_zero hhbar hmass)]

theorem longitudinal_transverse_eq_zero {k a : Vec} (hk : k ≠ 0)
    (hl : a ∈ longitudinalSpace k) (ht : a ∈ transverseSpace k) : a = 0 := by
  obtain ⟨s, rfl⟩ := Submodule.mem_span_singleton.mp hl
  change dot k (s • k) = 0 at ht
  have heq : dot k (s • k) = s * normSq k := by simp [dot, normSq]; ring
  rw [heq] at ht
  have hs := (mul_eq_zero.mp ht).resolve_right (ne_of_gt (dot_self_pos hk))
  simp [hs]

/-- No nonzero transverse velocity quadrature in this entire scalar phase family. -/
theorem no_transverse_velocity {k : Vec} (hk : k ≠ 0)
    (hbar mass A B : ℝ)
    (hc : cosAmplitude hbar mass B k ∈ transverseSpace k)
    (hs : sinAmplitude hbar mass A k ∈ transverseSpace k) :
    cosAmplitude hbar mass B k = 0 ∧ sinAmplitude hbar mass A k = 0 :=
  ⟨longitudinal_transverse_eq_zero hk (cosAmplitude_mem _ _ _ _) hc,
    longitudinal_transverse_eq_zero hk (sinAmplitude_mem _ _ _ _) hs⟩

theorem transverse_linearVelocity_zero {k : Vec} (hk : k ≠ 0)
    (hbar mass theta0 omega A B : ℝ)
    (hc : cosAmplitude hbar mass B k ∈ transverseSpace k)
    (hs : sinAmplitude hbar mass A k ∈ transverseSpace k) (x : Point) :
    linearVelocity hbar mass theta0 omega A B k x = 0 := by
  obtain ⟨hc0, hs0⟩ := no_transverse_velocity hk hbar mass A B hc hs
  simp [linearVelocity_eq, hc0, hs0]

theorem zero_wavevector (hbar mass theta0 omega A B : ℝ) (x : Point) :
    linearVelocity hbar mass theta0 omega A B 0 x = 0 := by
  simp [linearVelocity_eq, cosAmplitude, sinAmplitude]

theorem concrete_longitudinal_velocity :
    linearVelocity 1 1 0 1 0 1 ![0, 0, 1] 0 = ![0, 0, 1] ∧
      linearVelocity 1 1 0 1 0 1 ![0, 0, 1] 0 ≠ 0 := by
  norm_num [linearVelocity_eq, cosAmplitude, sinAmplitude, phase]

end
end S9Pilot.Madelung
