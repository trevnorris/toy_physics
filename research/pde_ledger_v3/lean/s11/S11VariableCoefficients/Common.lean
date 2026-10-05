import S11D4Odd.Calculus

/-! VC2: a weighted divergence retains the coefficient-gradient term. -/
namespace S11VariableCoefficients
noncomputable section
open S10Pilot
open scoped ContDiff

variable {D : ℕ}

def spatialDiv (J : Point D → Vec D) (x : Point D) : ℝ :=
  ∑ i : Fin D, coordDeriv i.succ (fun y => J y i) x

def gradientPair (a : Point D → ℝ) (J : Point D → Vec D) (x : Point D) : ℝ :=
  ∑ i : Fin D, coordDeriv i.succ a x * J x i

theorem weighted_divergence {a : Point D → ℝ} {J : Point D → Vec D}
    (ha : ContDiff ℝ ∞ a) (hJ : ∀ i, ContDiff ℝ ∞ (fun x => J x i)) (x : Point D) :
    spatialDiv (fun y i => a y * J y i) x = gradientPair a J x + a x * spatialDiv J x := by
  unfold spatialDiv gradientPair
  simp_rw [S11D4Odd.partial_mul ha (hJ _) _ _]
  rw [Finset.sum_add_distrib, Finset.mul_sum]

theorem weighted_density {a : Point D → ℝ} {J : Point D → Vec D} {q : Point D → ℝ}
    (ha : ContDiff ℝ ∞ a) (hJ : ∀ i, ContDiff ℝ ∞ (fun x => J x i))
    (hq : ∀ x, spatialDiv J x = q x) (x : Point D) :
    a x * q x = spatialDiv (fun y i => a y * J y i) x - gradientPair a J x := by
  rw [weighted_divergence ha hJ, hq]
  ring

def normalFlux (n : Vec D) (p : Fin D → Fin D → ℝ) : Vec D :=
  fun j => ∑ i, n i * p i j

def jumpPair (minus plus h : Vec D) : ℝ := ∑ j, (minus j - plus j) * h j

theorem jumpPair_zero_iff (minus plus : Vec D) :
    (∀ h, jumpPair minus plus h = 0) ↔ minus = plus := by
  constructor
  · intro h
    funext j
    have he := h (Pi.single j 1)
    simpa [jumpPair, Pi.single_apply, sub_eq_zero] using he
  · rintro rfl h
    simp [jumpPair]

theorem jumpPair_swap (minus plus h : Vec D) :
    jumpPair plus minus h = -jumpPair minus plus h := by
  simp only [jumpPair, ← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro j _
  ring

end
end S11VariableCoefficients
