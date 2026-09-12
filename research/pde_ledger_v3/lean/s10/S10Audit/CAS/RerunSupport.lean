import S10Audit.CAS.BasisCompletion
import S10Audit.CAS.LocusBindings

/-! Fixed-point transcripts use dimensionless coordinates p and a reduced
squared frequency z. Physical variables are k = κ • p and ω² = κ² z,
where κ has inverse-length units. Coefficients keep their physical units.
No nonzero literal is assigned a physical wavevector unit. -/

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section
set_option backward.isDefEq.respectTransparency false

def coordinateUnits : Symbol → Dim
  | .rho => units .rho
  | .mu => units .mu
  | .sigma => units .sigma
  | .z => dimensions 2 (-2) 0
  | .k _ => dimensions 0 0 0

theorem coordinate_frequency_units :
    dimensions (-1) 0 0 ^ (2 : ℕ) * coordinateUnits .z = units .z := by
  norm_num [coordinateUnits, units, dimensions_pow, dimensions_mul]

theorem coordinate_matrix_units :
    dimensions (-1) 0 0 ^ (2 : ℕ) * coordinateUnits .mu = dimensions (-3) (-2) 1 := by
  norm_num [coordinateUnits, units, dimensions_pow, dimensions_mul]

theorem coordinate_determinant_units :
    dimensions (-1) 0 0 ^ (6 : ℕ) * coordinateUnits .mu ^ (3 : ℕ) = dimensions (-9) (-6) 3 := by
  norm_num [coordinateUnits, units, dimensions_pow, dimensions_mul]

theorem coordinate_product_units :
    dimensions (-1) 0 0 ^ (3 : ℕ) * coordinateUnits .mu = dimensions (-4) (-2) 1 := by
  norm_num [coordinateUnits, units, dimensions_pow, dimensions_mul]

theorem coordinate_wavevector_units (i : Fin 3) :
    dimensions (-1) 0 0 * coordinateUnits (.k i) = units (.k i) := by
  norm_num [coordinateUnits, units, dimensions_mul]

theorem coordinate_residual_units :
    dimensions (-1) 0 0 ^ (2 : ℕ) * dimensions 0 0 0 = dimensions (-2) 0 0 := by
  norm_num [dimensions_pow, dimensions_mul]

theorem normSq_scale (c : ℝ) (p : Vec 3) : normSq (c • p) = c ^ 2 * normSq p := by
  simp only [normSq, dot, Fin.sum_univ_three, Pi.smul_apply, smul_eq_mul]
  ring

theorem referenceMatrix_scale (rho mu sigma z c : ℝ) (p : Vec 3) :
    referenceMatrix rho mu sigma (c ^ 2 * z) (c • p) =
      c ^ 2 • referenceMatrix rho mu sigma z p := by
  ext i j
  simp only [referenceMatrix, normSq_scale, Matrix.smul_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem referenceRoot_scale (rho mu sigma c : ℝ) (p : Vec 3) (r : Fin 3) :
    referenceRoot rho mu sigma (c • p) r = c ^ 2 * referenceRoot rho mu sigma p r := by
  fin_cases r <;>
    simp [referenceRoot, coneValue, extraConeValue, extraValue, extraNumerator,
      perpSq, normSq_scale, Pi.smul_apply, smul_eq_mul] <;> ring

theorem coordinate_determinant_scale (rho mu sigma z c : ℝ) (p : Vec 3) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma (c ^ 2 * z) (c • p)) =
      c ^ 6 * Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z p) := by
  rw [referenceMatrix_scale]
  simp only [Matrix.det_fin_three, Matrix.smul_apply, smul_eq_mul]
  ring

theorem coordinate_kernel_scale (rho mu sigma z c : ℝ) (p a : Vec 3) (hc : c ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma (c ^ 2 * z) (c • p)).mulVec a = 0 ↔
      ((1/2 : ℝ) • referenceMatrix rho mu sigma z p).mulVec a = 0 := by
  rw [referenceMatrix_scale, smul_comm (1/2 : ℝ) (c ^ 2), Matrix.smul_mulVec, smul_eq_zero]
  exact or_iff_right (pow_ne_zero 2 hc)

theorem coordinate_product_scale (rho mu sigma z c : ℝ) (p : Vec 3) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma (c ^ 2 * z) (c • p)).mulVec (c • p) =
      c ^ 3 • (((1/2 : ℝ) • referenceMatrix rho mu sigma z p).mulVec p) := by
  rw [referenceMatrix_scale]
  ext i
  simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three, Matrix.smul_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem coordinate_dot_scale (c : ℝ) (p a : Vec 3) : dot (c • p) a = c * dot p a := by
  exact dot_smul_left p a c

theorem coordinate_residual_scale (c : ℝ) (p a : Vec 3) :
    normSq (c • p) • a - dot (c • p) a • (c • p) =
      c ^ 2 • (normSq p • a - dot p a • p) := by
  rw [normSq_scale, coordinate_dot_scale]
  ext i
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  ring

def rootMode (sigma : ℝ) (p : Vec 3) (r : Fin 3) : Submodule ℝ (Vec 3) :=
  modeSpace 0 sigma (if r = 0 then 0 else if r = 1 then normSq p else extraValue 0 sigma p) p

theorem coordinate_kernel_mode (rho mu sigma : ℝ) (p a : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) :
    ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma p r) p).mulVec a = 0 ↔
      a ∈ rootMode sigma p r := by
  rw [Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),
    referenceMatrix_normalized rho mu sigma _ p a hm, smul_eq_zero, or_iff_right hm,
    referenceRoot_normalized rho mu sigma p r hr hm]
  rfl

/-- Membership plus independence and the entire kernel dimension certifies
all basis vectors and rules out an incomplete sample from a repeated root. -/
theorem complete_of_dimension {n : ℕ} (b : Fin n → Vec 3) (V : Submodule ℝ (Vec 3))
    (hli : LinearIndependent ℝ b) (hmem : ∀ j, b j ∈ V)
    (hdim : Module.finrank ℝ V = n) : Submodule.span ℝ (Set.range b) = V := by
  apply Submodule.eq_of_le_of_finrank_eq
  · exact Submodule.span_le.mpr (by rintro a ⟨j, rfl⟩; exact hmem j)
  · rw [finrank_span_eq_card hli, Fintype.card_fin, hdim]

end
end S10Audit.CAS
