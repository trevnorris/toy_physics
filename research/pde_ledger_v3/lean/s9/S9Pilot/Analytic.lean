import S9Pilot.PlaneWave
import Mathlib.Analysis.Calculus.LineDeriv.IntegrationByParts
import Mathlib.Analysis.Distribution.AEEqOfIntegralContDiff

/-! Analytic preliminaries for compactly supported variations on spacetime R^4. -/

namespace S9Pilot
noncomputable section
open MeasureTheory
open scoped ContDiff

def SmoothField (u : Point → Vec) : Prop :=
  ∀ i, ContDiff ℝ ∞ (fun x => u x i)

/-- With finitely many components this is equivalent to compact support of the vector field. -/
def TestField (h : Point → Vec) : Prop :=
  SmoothField h ∧ ∀ i, HasCompactSupport (fun x => h x i)

theorem coordDeriv_eq_fderiv {f : Point → ℝ} (hf : Differentiable ℝ f)
    (j : Fin 4) (x : Point) : coordDeriv j f x = fderiv ℝ f x (axis j) :=
  (hf x).lineDeriv_eq_fderiv

theorem smooth_coordDeriv {f : Point → ℝ} (hf : ContDiff ℝ ∞ f) (j : Fin 4) :
    ContDiff ℝ ∞ (coordDeriv j f) := by
  have hd := hf.differentiable (by simp)
  have heq : coordDeriv j f = fun x => fderiv ℝ f x (axis j) :=
    funext (coordDeriv_eq_fderiv hd j)
  rw [heq]
  exact (hf.fderiv_right (by simp)).clm_apply contDiff_const

theorem compact_coordDeriv {f : Point → ℝ} (hf : ContDiff ℝ ∞ f)
    (hc : HasCompactSupport f) (j : Fin 4) : HasCompactSupport (coordDeriv j f) := by
  have heq : coordDeriv j f = fun x => fderiv ℝ f x (axis j) :=
    funext (coordDeriv_eq_fderiv (hf.differentiable (by simp)) j)
  rw [heq]
  exact hc.fderiv_apply ℝ (axis j)

theorem coordDeriv_add_smul {f g : Point → ℝ}
    (hf : Differentiable ℝ f) (hg : Differentiable ℝ g) (s : ℝ) (j : Fin 4) (x : Point) :
    coordDeriv j (fun y => f y + s * g y) x = coordDeriv j f x + s * coordDeriv j g x := by
  rw [coordDeriv_eq_fderiv (hf.fun_add (hg.const_mul s)),
    coordDeriv_eq_fderiv hf, coordDeriv_eq_fderiv hg]
  simp [fderiv_fun_add (hf x) ((hg.const_mul s) x), fderiv_const_mul (hg x)]

theorem fieldJet_add_smul {u h : Point → Vec}
    (hu : SmoothField u) (hh : SmoothField h) (s : ℝ) (x : Point) :
    fieldJet (u + s • h) x = fieldJet u x + s • fieldJet h x := by
  ext j i
  exact coordDeriv_add_smul ((hu i).differentiable (by simp))
    ((hh i).differentiable (by simp)) s j x

theorem smooth_momentum {u : Point → Vec} (hu : SmoothField u)
    (rho mu : ℝ) (j : Fin 4) (i : Fin 3) :
    ContDiff ℝ ∞ (fun x => momentum rho mu (fieldJet u x) j i) := by
  have hJ : ∀ a b, ContDiff ℝ ∞ (fun x => fieldJet u x a b) :=
    fun a b => smooth_coordDeriv (hu b) a
  simp only [momentum_eq, dot, jetCurl, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  fun_prop

theorem smooth_eulerLagrange {u : Point → Vec} (hu : SmoothField u)
    (rho mu : ℝ) (i : Fin 3) :
    ContDiff ℝ ∞ (fun x => eulerLagrange rho mu u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞
      (coordDeriv j (fun y => momentum rho mu (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu rho mu j i) j
  unfold eulerLagrange
  fun_prop

theorem integrable_mul_compact {f g : Point → ℝ}
    (hf : Continuous f) (hg : Continuous g) (hc : HasCompactSupport g) :
    Integrable (fun x => f x * g x) := by
  exact (hf.mul hg).integrable_of_hasCompactSupport hc.mul_left

/-- Boundary terms vanish by actual integration by parts and compact support. -/
theorem coord_integration_by_parts {f g : Point → ℝ}
    (hf : ContDiff ℝ ∞ f) (hg : ContDiff ℝ ∞ g) (hc : HasCompactSupport g) (j : Fin 4) :
    (∫ x, f x * coordDeriv j g x) = -(∫ x, coordDeriv j f x * g x) := by
  have hfd := hf.differentiable (by simp)
  have hgd := hg.differentiable (by simp)
  have A := integrable_mul_compact (smooth_coordDeriv hf j).continuous hg.continuous hc
  have B := integrable_mul_compact hf.continuous (smooth_coordDeriv hg j).continuous
    (compact_coordDeriv hg hc j)
  have C := integrable_mul_compact hf.continuous hg.continuous hc
  simp_rw [coordDeriv_eq_fderiv hfd, coordDeriv_eq_fderiv hgd] at A B ⊢
  exact integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable A B C
    (fun x _ => hfd x) (fun x _ => hgd x)

end
end S9Pilot
