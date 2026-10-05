import S11D5Bulk.Census

/-! D5B.4: nonzero action/momentum/current and modal witnesses.
Each paired external negative control changes a true value or proposition. -/
namespace S11D5Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

def evenJet : Jet 5 := Fin.cases (fun _ => 0)
  (Matrix.of ![![1,0,0,0,0],![0,1,0,0,0],![0,0,0,0,0],![0,0,0,0,0],![0,0,0,0,0]])

theorem even_null_density_nonzero : lagrangian ![1,-1,0] evenJet = -1 := by
  norm_num [lagrangian, divergence, transposePair,
    gradientSquare, evenJet, S11D5Invariants.sum_five,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem even_momentum_normalization : momentum ![1,0,0] evenJet 1 0 = -2 := by
  rw [momentum_eq]
  simp only [divergence, S11D5Invariants.sum_five]
  change (-(1 * (1+1+0+0+0) + 0 * 1 + 0 * 1) : ℝ) = -2
  norm_num


theorem partial_coordinate (i j : Fin 6) (x : Point 5) :
    coordDeriv i (fun y => y j) x = if j = i then 1 else 0 := by
  change deriv (fun s : ℝ => x j + s * axis i j) 0 = _
  have h := ((hasDerivAt_const (0 : ℝ) (x j)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (axis i j))).deriv
  calc
    _ = 0 + 1 * axis i j := h
    _ = _ := by simp [axis, Pi.single_apply]

def affineWitness (x : Point 5) : Vec 5 := ![x 1,x 2,0,0,0]

theorem affineWitness_smooth : SmoothField affineWitness := by
  intro i
  fin_cases i <;> dsimp [affineWitness] <;> fun_prop

theorem affineWitness_jet (x : Point 5) : fieldJet affineWitness x = evenJet := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [fieldJet, affineWitness, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four,
      partial_coordinate, S11D4Odd.partial_zero, Fin.ext_iff] <;> rfl


theorem even_current_normalization (x : Point 5) :
    (∑ i : Fin 5, coordDeriv i.succ (fun y => boundaryCurrent affineWitness y i) x) = 2 := by
  rw [boundary_identity affineWitness_smooth, affineWitness_jet]
  norm_num [divergence, transposePair, evenJet, S11D5Invariants.sum_five, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem transverse_response :
    eulerLagrange ![2,3,7] (planeWave 0 ![1,0,0,0,0] ![0,1,0,0,0]) 0 1 = -7 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, normSq, dot, phase, S11D5Invariants.sum_five,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem longitudinal_response :
    eulerLagrange ![2,3,7] (planeWave 0 ![1,0,0,0,0] ![1,0,0,0,0]) 0 0 = -12 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, normSq, dot, phase, S11D5Invariants.sum_five,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem fifth_transverse_response :
    eulerLagrange ![2,3,7] (planeWave 0 ![0,0,0,0,1] ![1,0,0,0,0]) 0 0 = -7 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, normSq, dot, phase, S11D5Invariants.sum_five,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem fifth_longitudinal_response :
    eulerLagrange ![2,3,7] (planeWave 0 ![0,0,0,0,1] ![0,0,0,0,1]) 0 4 = -12 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, normSq, dot, phase, S11D5Invariants.sum_five,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four, Fin.cases]

theorem null_positive (t : ℝ) : VariationallyNull ![t,-t,0] := by
  rw [variationallyNull_iff]
  simp [Matrix.cons_val_two]

theorem null_firstVariation {u h : Point 5 → Vec 5}
    (hu : SmoothField u) (hh : TestField h) (t : ℝ) :
    deriv (relativeAction ![t,-t,0] u h) 0 = 0 :=
  null_positive t u hu h hh

end
end S11D5Bulk
