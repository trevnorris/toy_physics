import S10Pilot.Analytic
import Mathlib.Analysis.Calculus.FDeriv.Symmetric

/-! Coordinate rules needed for D4B.2–D4B.3, on the existing spacetime. -/
namespace S11D4Odd
noncomputable section
open S10Pilot
open scoped ContDiff

variable {D : ℕ}

theorem partial_add {f g : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (hg : ContDiff ℝ ∞ g) (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => f y + g y) x = coordDeriv j f x + coordDeriv j g x := by
  rw [coordDeriv_eq_fderiv ((hf.add hg).differentiable (by simp)),
    coordDeriv_eq_fderiv (hf.differentiable (by simp)),
    coordDeriv_eq_fderiv (hg.differentiable (by simp))]
  simp [fderiv_fun_add (hf.differentiable (by simp) x) (hg.differentiable (by simp) x)]

theorem partial_sub {f g : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (hg : ContDiff ℝ ∞ g) (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => f y - g y) x = coordDeriv j f x - coordDeriv j g x := by
  rw [coordDeriv_eq_fderiv ((hf.sub hg).differentiable (by simp)),
    coordDeriv_eq_fderiv (hf.differentiable (by simp)),
    coordDeriv_eq_fderiv (hg.differentiable (by simp))]
  simp [fderiv_fun_sub (hf.differentiable (by simp) x) (hg.differentiable (by simp) x)]

theorem partial_neg {f : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => -f y) x = -coordDeriv j f x := by
  rw [coordDeriv_eq_fderiv (hf.neg.differentiable (by simp)),
    coordDeriv_eq_fderiv (hf.differentiable (by simp))]
  simp [fderiv_fun_neg]

theorem partial_const_mul {f : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (c : ℝ) (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => c * f y) x = c * coordDeriv j f x := by
  rw [coordDeriv_eq_fderiv ((contDiff_const.mul hf).differentiable (by simp)),
    coordDeriv_eq_fderiv (hf.differentiable (by simp))]
  simp [fderiv_const_mul (hf.differentiable (by simp) x)]

theorem partial_mul {f g : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (hg : ContDiff ℝ ∞ g) (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => f y * g y) x =
      coordDeriv j f x * g x + f x * coordDeriv j g x := by
  rw [coordDeriv_eq_fderiv ((hf.mul hg).differentiable (by simp)),
    coordDeriv_eq_fderiv (hf.differentiable (by simp)),
    coordDeriv_eq_fderiv (hg.differentiable (by simp))]
  simp [fderiv_fun_mul (hf.differentiable (by simp) x) (hg.differentiable (by simp) x)]
  ring

theorem partial_zero (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun _ => (0 : ℝ)) x = 0 := by simp [coordDeriv]

theorem partial_sum {n : ℕ} (f : Fin n → Point D → ℝ)
    (hf : ∀ i, ContDiff ℝ ∞ (f i)) (j : Fin (D + 1)) (x : Point D) :
    coordDeriv j (fun y => ∑ i, f i y) x = ∑ i, coordDeriv j (f i) x := by
  have hs : ContDiff ℝ ∞ (fun y => ∑ i, f i y) := by fun_prop
  rw [coordDeriv_eq_fderiv (hs.differentiable (by simp)),
    fderiv_fun_sum (fun i _ => (hf i).differentiable (by simp) x)]
  simp only [sum_apply]
  apply Finset.sum_congr rfl
  intro i _
  exact (coordDeriv_eq_fderiv ((hf i).differentiable (by simp)) j x).symm

theorem partial_second {f : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (i j : Fin (D + 1)) (x : Point D) :
    coordDeriv i (coordDeriv j f) x = fderiv ℝ (fderiv ℝ f) x (axis i) (axis j) := by
  rw [coordDeriv_eq_fderiv ((smooth_coordDeriv hf j).differentiable (by simp))]
  have he : coordDeriv j f = fun y => fderiv ℝ f y (axis j) :=
    funext (coordDeriv_eq_fderiv (hf.differentiable (by simp)) j)
  have hdf : ContDiff ℝ ∞ (fderiv ℝ f) := hf.fderiv_right (by simp)
  rw [he, fderiv_clm_apply (hdf.differentiable (by simp) x)
    (differentiableAt_const (axis j))]
  simp

theorem partial_commute {f : Point D → ℝ} (hf : ContDiff ℝ ∞ f)
    (i j : Fin (D + 1)) (x : Point D) :
    coordDeriv i (coordDeriv j f) x = coordDeriv j (coordDeriv i f) x := by
  rw [partial_second hf, partial_second hf]
  exact (hf.contDiffAt.isSymmSndFDerivAt (by
    simp only [minSmoothness_of_isRCLikeNormedField]
    exact WithTop.coe_le_coe.mpr le_top)).eq (axis i) (axis j)

end
end S11D4Odd
