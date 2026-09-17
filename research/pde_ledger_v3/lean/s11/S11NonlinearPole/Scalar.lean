import S11NonlinearPole.Moments

/-! NP3: scalar pencils, their actual derivatives, and contour integrals. -/
namespace S11NonlinearPole
noncomputable section
open Complex Metric

def squarePencil (z : ℂ) : ℂ := z ^ 2

def squareInverse (z : ℂ) : ℂ := (squarePencil z)⁻¹

def squareLog (z : ℂ) : ℂ := squareInverse z * deriv squarePencil z

theorem square_derivative (z : ℂ) : deriv squarePencil z = 2 * z := by
  unfold squarePencil
  simp

theorem square_inverse_identity {z : ℂ} (hz : z ≠ 0) :
    squarePencil z * squareInverse z = 1 := by
  exact mul_inv_cancel₀ (pow_ne_zero 2 hz)

theorem square_log_formula {z : ℂ} (hz : z ≠ 0) : squareLog z = 2 * z⁻¹ := by
  rw [squareLog, square_derivative]
  simp only [squareInverse, squarePencil]
  field_simp

theorem square_residue (R : ℝ) (hR : 0 < R) : moment R squareInverse = 0 := by
  change moment R (fun z : ℂ => (z ^ 2)⁻¹) = 0
  simpa [zpow_neg] using moment_zpow R hR (-2)

theorem square_log_moment (R : ℝ) (hR : 0 < R) : moment R squareLog = 2 := by
  rw [moment_congr hR.le (fun _ hz => square_log_formula (nonzero_on_circle hR hz))]
  change moment R (fun z : ℂ => (2 : ℂ) • z⁻¹) = _
  rw [moment_smul]
  have hi := moment_zpow R hR (-1)
  simpa using hi

theorem square_log_not_idempotent (R : ℝ) (hR : 0 < R) :
    moment R squareLog * moment R squareLog ≠ moment R squareLog := by
  rw [square_log_moment R hR]
  norm_num

theorem square_response_nonzero : squareInverse 1 = 1 := by
  norm_num [squareInverse, squarePencil]

theorem square_higher_coefficient (R : ℝ) (hR : 0 < R) :
    moment R (fun z => z * squareInverse z) = 1 := by
  have heq : ∀ z : ℂ, z ≠ 0 → z * squareInverse z = z ^ (-1 : ℤ) := by
    intro z hz
    simp only [squareInverse, squarePencil, zpow_neg_one]
    field_simp
  rw [moment_congr hR.le (fun z hz => heq z (nonzero_on_circle hR hz))]
  simpa using moment_zpow R hR (-1)

def twoRootPencil (z : ℂ) : ℂ := z ^ 2 - 1

def twoRootLog (z : ℂ) : ℂ := (twoRootPencil z)⁻¹ * deriv twoRootPencil z

theorem twoRoot_derivative (z : ℂ) : deriv twoRootPencil z = 2 * z := by
  unfold twoRootPencil
  simp

theorem twoRoot_zeros (z : ℂ) : twoRootPencil z = 0 ↔ z = 1 ∨ z = -1 := by
  unfold twoRootPencil
  constructor
  · intro h
    have hf : (z - 1) * (z + 1) = 0 := by
      calc (z - 1) * (z + 1) = z ^ 2 - 1 := by ring
           _ = 0 := h
    simpa only [sub_eq_zero, add_eq_zero_iff_eq_neg] using mul_eq_zero.mp hf
  · rintro (rfl | rfl) <;> norm_num

theorem twoRoot_log_formula {z : ℂ} (hp : z ≠ 1) (hm : z ≠ -1) :
    twoRootLog z = (z - 1)⁻¹ + (z + 1)⁻¹ := by
  have hp' : z - 1 ≠ 0 := sub_ne_zero.mpr hp
  have hm' : z + 1 ≠ 0 := fun h => hm (add_eq_zero_iff_eq_neg.mp h)
  have hq : z ^ 2 - 1 ≠ 0 := by
    rw [show z ^ 2 - 1 = (z - 1) * (z + 1) by ring]
    exact mul_ne_zero hp' hm'
  rw [twoRootLog, twoRoot_derivative]
  simp only [twoRootPencil]
  field_simp
  ring

theorem twoRoot_log_moment : moment 2 twoRootLog = 2 := by
  have hp : (1 : ℂ) ∉ sphere 0 (2 : ℝ) := by norm_num [mem_sphere_iff_norm]
  have hm : (-1 : ℂ) ∉ sphere 0 (2 : ℝ) := by norm_num [mem_sphere_iff_norm]
  rw [moment_congr (by norm_num : (0 : ℝ) ≤ 2) (fun z hz =>
    twoRoot_log_formula (ne_of_mem_of_not_mem hz hp) (ne_of_mem_of_not_mem hz hm))]
  have ip : CircleIntegrable (fun z : ℂ => (z - 1)⁻¹) 0 2 :=
    circleIntegrable_sub_inv_iff.mpr (Or.inr (by norm_num [mem_sphere_iff_norm]))
  have im : CircleIntegrable (fun z : ℂ => (z + 1)⁻¹) 0 2 := by
    simpa only [sub_neg_eq_add] using
      (circleIntegrable_sub_inv_iff (c := (0 : ℂ)) (w := -1) (R := 2)).mpr
        (Or.inr (by norm_num [mem_sphere_iff_norm]))
  have jp := circleIntegral.integral_sub_inv_of_mem_ball
    (c := (0 : ℂ)) (w := 1) (R := 2) (by norm_num [mem_ball_iff_norm])
  have jm := circleIntegral.integral_sub_inv_of_mem_ball
    (c := (0 : ℂ)) (w := -1) (R := 2) (by norm_num [mem_ball_iff_norm])
  simp only [sub_neg_eq_add] at jm
  simp only [moment, circleIntegral.integral_add ip im, jp, jm, smul_eq_mul]
  field_simp [normalization_ne_zero]
  norm_num

theorem twoRoot_log_not_idempotent :
    moment 2 twoRootLog * moment 2 twoRootLog ≠ moment 2 twoRootLog := by
  rw [twoRoot_log_moment]
  norm_num

/-- Actual observation and forcing vary holomorphically with frequency. -/
def affineResponse (z : ℂ) : ℂ := (2 + 5*z) * squareInverse z * (1 + 3*z)

theorem affine_response_expansion {z : ℂ} (hz : z ≠ 0) :
    affineResponse z = ∑ k : Fin 3, z ^ (![-2,-1,0] k : ℤ) * (![2,11,15] k : ℂ) := by
  simp [affineResponse, squareInverse, squarePencil, Fin.sum_univ_succ, zpow_neg]
  field_simp
  ring

theorem affine_response_residue (R : ℝ) (hR : 0 < R) : moment R affineResponse = 11 := by
  rw [moment_congr hR.le (fun z hz => affine_response_expansion (nonzero_on_circle hR hz))]
  change moment R (fun z : ℂ => ∑ k : Fin 3, z ^ (![-2,-1,0] k : ℤ) • (![2,11,15] k : ℂ)) = _
  rw [moment_finite_laurent Finset.univ R hR]
  norm_num [Fin.sum_univ_succ]

end
end S11NonlinearPole
