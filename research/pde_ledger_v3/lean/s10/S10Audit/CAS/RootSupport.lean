import S10Audit.CAS.GenericCountSupport
import Mathlib.Algebra.Polynomial.Roots

namespace S10Audit.CAS
open S10Pilot S10Anisotropic Polynomial
noncomputable section
set_option backward.isDefEq.respectTransparency false

/-- Syntactic counterpart of the emitted FreeQ[expression, omegaSquared] filter. -/
def frequencyFree : Expr Symbol → Bool
  | .atom a => decide (a ≠ .z)
  | .scalar _ => true
  | .add a b | .sub a b | .mul a b | .div a b => frequencyFree a && frequencyFree b
  | .pow a _ => frequencyFree a
  | .sum _ terms => (List.ofFn (fun i => frequencyFree (terms i))).all id

/-- The cubic determinant polynomial, including its non-monic coefficient. -/
def rootPolynomial (rho mu sigma : ℝ) (k : Vec 3) : Polynomial ℝ :=
  C (rho ^ 3 * sigma / 8) *
    ((X - C 0) * (X - C (referenceRoot rho mu sigma k 1)) *
      (X - C (referenceRoot rho mu sigma k 2)))

theorem rootPolynomial_eval (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (rootPolynomial rho mu sigma k).eval z =
      Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z k) := by
  simp only [rootPolynomial, eval_mul, eval_sub, eval_X, eval_C]
  cas_equal

theorem rootPolynomial_ne_zero (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) : rootPolynomial rho mu sigma k ≠ 0 := by
  simp [rootPolynomial, hr, hs, X_sub_C_ne_zero]

theorem rootPolynomial_roots (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (rootPolynomial rho mu sigma k).roots =
      ([0, referenceRoot rho mu sigma k 1, referenceRoot rho mu sigma k 2] : List ℝ) := by
  rw [rootPolynomial, roots_C_mul _ (by positivity)]
  rw [roots_mul (mul_ne_zero (mul_ne_zero (X_sub_C_ne_zero _) (X_sub_C_ne_zero _)) (X_sub_C_ne_zero _)),
    roots_mul (mul_ne_zero (X_sub_C_ne_zero _) (X_sub_C_ne_zero _))]
  simp only [roots_X_sub_C]
  rfl

theorem rootPolynomial_complete (rho mu sigma z : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z k) = 0 ↔
      z ∈ (rootPolynomial rho mu sigma k).roots := by
  rw [mem_roots (rootPolynomial_ne_zero rho mu sigma k hr hs)]
  change _ = 0 ↔ (rootPolynomial rho mu sigma k).eval z = 0
  rw [rootPolynomial_eval rho mu sigma z k hr hs]

theorem rootPolynomial_card (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (rootPolynomial rho mu sigma k).roots.card = 3 := by
  rw [rootPolynomial_roots rho mu sigma k hr hs]
  rfl

theorem root_nonzero (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hk : k ≠ 0) :
    referenceRoot rho mu sigma k 1 ≠ 0 ∧ referenceRoot rho mu sigma k 2 ≠ 0 := by
  exact ⟨mul_ne_zero (div_ne_zero hm hr) (ne_of_gt (dot_self_pos hk)),
    mul_ne_zero (div_ne_zero hm hr) (ne_of_gt (extraValue_pos hs hk))⟩

theorem root_coincidence (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    referenceRoot rho mu sigma k 2 = referenceRoot rho mu sigma k 1 ↔ perpSq 0 k = 0 := by
  change (mu / rho) * extraValue 0 sigma k = (mu / rho) * normSq k ↔ _
  rw [mul_right_inj' (div_ne_zero hm hr), frequency_coincidence_iff hs hs1]

theorem split_root_nodup (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)
    (hk : k ≠ 0) (hq : perpSq 0 k ≠ 0) :
    ([0, referenceRoot rho mu sigma k 1, referenceRoot rho mu sigma k 2] : List ℝ).Nodup := by
  have hn := root_nonzero rho mu sigma k hr hm hs hk
  have hne : referenceRoot rho mu sigma k 2 ≠ referenceRoot rho mu sigma k 1 :=
    fun he => hq ((root_coincidence rho mu sigma k hr hm (ne_of_gt hs) hs1).mp he)
  simp [hn.1.symm, hn.2.symm, hne.symm]

theorem parallel_root_nodup (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hk : k ≠ 0) :
    ([0, referenceRoot rho mu sigma k 1] : List ℝ).Nodup := by
  simp [(root_nonzero rho mu sigma k hr hm hs hk).1.symm]

theorem split_root_multiplicities (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)
    (hk : k ≠ 0) (hq : perpSq 0 k ≠ 0) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma k 1) (rootPolynomial rho mu sigma k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma k 2) (rootPolynomial rho mu sigma k) = 1 := by
  classical
  have hn := root_nonzero rho mu sigma k hr hm hs hk
  have hne : referenceRoot rho mu sigma k 2 ≠ referenceRoot rho mu sigma k 1 :=
    fun he => hq ((root_coincidence rho mu sigma k hr hm (ne_of_gt hs) hs1).mp he)
  simp [← count_roots, rootPolynomial_roots rho mu sigma k hr (ne_of_gt hs),
    hn.1, hn.2, hn.1.symm, hn.2.symm, hne, hne.symm]

theorem parallel_root_multiplicities (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)
    (hk : k ≠ 0) (hq : perpSq 0 k = 0) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma k 1) (rootPolynomial rho mu sigma k) = 2 := by
  classical
  have hn := (root_nonzero rho mu sigma k hr hm hs hk).1
  have he := (root_coincidence rho mu sigma k hr hm (ne_of_gt hs) hs1).mpr hq
  simp [← count_roots, rootPolynomial_roots rho mu sigma k hr (ne_of_gt hs), he, hn, hn.symm]

theorem genericChart_split (sigma : ℝ) (k : Vec 3) (h : GenericChart sigma k) :
    k ≠ 0 ∧ perpSq 0 k ≠ 0 := by
  rcases h with ⟨_hs, _hs1, h0, h1, _h2⟩
  refine ⟨by intro heq; exact h0 (congrFun heq 0), ?_⟩
  have hpos := sq_pos_of_ne_zero h1
  simp only [perpSq, normSq, dot, Fin.sum_univ_three]
  nlinarith [sq_nonneg (k 2)]

end
end S10Audit.CAS
