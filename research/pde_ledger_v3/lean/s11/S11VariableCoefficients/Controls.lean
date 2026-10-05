import S11VariableCoefficients.D3
import S11VariableCoefficients.D4
import S11VariableCoefficients.Interface

/-! VC4: actual smooth-field and integral witnesses, with nonzero corrections. -/
namespace S11VariableCoefficients
noncomputable section
open S10Pilot MeasureTheory
open scoped ContDiff

theorem partial_coordinate {D : ℕ} (i j : Fin (D+1)) (x : Point D) :
    coordDeriv i (fun y => y j) x = if j = i then 1 else 0 := by
  change deriv (fun s : ℝ => x j + s * axis i j) 0 = _
  have h := ((hasDerivAt_const (0 : ℝ) (x j)).add ((hasDerivAt_id (0 : ℝ)).mul_const (axis i j))).deriv
  calc
    _ = 0 + 1 * axis i j := h
    _ = _ := by simp [axis, Pi.single_apply]

def d3Field (x : Point 3) : Vec 3 := ![0,x 2,0]

theorem d3Field_smooth : SmoothField d3Field := by
  intro i
  fin_cases i <;> dsimp [d3Field] <;> fun_prop

theorem d3_nonzero_response (x : Point 3) :
    D3.eulerLagrange (fun y => ![y 1,-y 1,0]) d3Field x 0 = 1 := by
  rw [D3.null_profile_residual (by fun_prop) d3Field_smooth]
  norm_num [d3Field, S11D3Bulk.divergence, fieldJet, Fin.sum_univ_three,
    Matrix.cons_val_two, partial_coordinate, S11D4Odd.partial_zero]
  decide

def d4Field (x : Point 4) : Vec 4 := ![0,0,0,x 3]

theorem d4Field_smooth : SmoothField d4Field := by
  intro i
  fin_cases i <;> dsimp [d4Field] <;> fun_prop

theorem d4_nonzero_response (x : Point 4) :
    D4.eulerLagrange (fun y => y 1) d4Field x 1 = 1 / 2 := by
  rw [D4.eulerLagrange_eq (by fun_prop) d4Field_smooth]
  norm_num [d4Field, Fin.sum_univ_four, partial_coordinate,
    S11D4Odd.dualCurl, antisym, fieldJet, S11D4Odd.partial_zero,
    Matrix.cons_val_two, Matrix.cons_val_three]
  norm_num [Fin.ext_iff]

theorem nonzero_weighted_correction (x : Point 1) :
    gradientPair (fun y => y 1) (fun _ => ![1]) x = 1 := by
  rw [gradientPair, Fin.sum_univ_one]
  norm_num [partial_coordinate, Fin.succ, Fin.ext_iff]

def interfaceWitness : ℝ :=
  (∫ x in (-1 : ℝ)..0, 2 * (-2 * x)) + (∫ x in (0 : ℝ)..1, 5 * (-2 * x))

theorem interfaceWitness_eq : interfaceWitness = -3 := by
  have hm : ∀ x ∈ Set.uIcc (-1 : ℝ) 0, HasDerivAt (fun _ : ℝ => (2 : ℝ)) 0 x :=
    fun x _ => hasDerivAt_const x 2
  have hp : ∀ x ∈ Set.uIcc (0 : ℝ) 1, HasDerivAt (fun _ : ℝ => (5 : ℝ)) 0 x :=
    fun x _ => hasDerivAt_const x 5
  have hh (x : ℝ) : HasDerivAt (fun y : ℝ => 1-y^2) (-2*x) x := by
    convert! (hasDerivAt_const x (1 : ℝ)).sub ((hasDerivAt_id x).pow 2) using 1
    simp
  have e := compact_endpoint_split hm hp (fun x _ => hh x) (fun x _ => hh x)
    (intervalIntegrable_const) (intervalIntegrable_const)
    ((by fun_prop : Continuous (fun x : ℝ => -2*x)).intervalIntegrable _ _)
    ((by fun_prop : Continuous (fun x : ℝ => -2*x)).intervalIntegrable _ _)
    (by norm_num) (by norm_num)
  norm_num [interfaceWitness] at e ⊢
  exact e

theorem d3_traction_nonzero :
    D3.traction ![1,-1,0] (Fin.cases (fun _ => 0)
      (Matrix.of ![![0,0,0],![0,1,0],![0,0,0]])) ![1,0,0] 0 = -1 := by
  rw [D3.traction_eq]
  simp only [S11D3Bulk.divergence, Fin.sum_univ_three]
  change -(1 * (1 * (0 + 1 + 0) + (-1) * 0 + 0 * 0) + 0 * _ + 0 * _) = (-1 : ℝ)
  norm_num

theorem d4_traction_nonzero :
    D4.traction 2 S11D4Odd.witnessJet ![1,0,0,0] 1 = -1 := by
  rw [D4.traction_eq]
  simp only [Fin.sum_univ_four]
  change -(2 / 2 : ℝ) * (1 * (1 - 0) + 0 * _ + 0 * _ + 0 * _) = -1
  norm_num

theorem trace_jump_nonzero : jumpPair (![2] : Vec 1) ![5] ![1] = -3 := by
  norm_num [jumpPair]

theorem matched_trace_zero (v h : Vec 3) : jumpPair v v h = 0 := by
  simp [jumpPair]

end
end S11VariableCoefficients
