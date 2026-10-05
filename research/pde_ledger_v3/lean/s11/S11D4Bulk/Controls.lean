import S11D4Bulk.Census

/-! D4C.4: nonzero action/momentum/current and modal witnesses.
Each paired external negative control changes a true value or proposition. -/
namespace S11D4Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

def evenJet : Jet 4 := Fin.cases (fun _ => 0)
  (Matrix.of ![![1,0,0,0],![0,1,0,0],![0,0,0,0],![0,0,0,0]])

theorem even_null_density_nonzero : lagrangian ![1,-1,0,0] evenJet = -1 := by
  norm_num [lagrangian, Even.lagrangian, Even.divergence, Even.transposePair,
    Even.gradientSquare, S11D4Odd.lagrangian, evenJet, Fin.sum_univ_four,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Fin.cases]

theorem even_momentum_normalization : momentum ![1,0,0,0] evenJet 1 0 = -2 := by
  rw [momentum_eq, Even.momentum_eq, S11D4Odd.momentum_eq]
  simp only [Even.divergence, Fin.sum_univ_four]
  change (-(1 * (1+1+0+0) + 0 * 1 + 0 * 1) + -(0 / 2) * 0 : ℝ) = -2
  norm_num


theorem odd_density_normalization :
    lagrangian ![0,0,0,2] S11D4Odd.witnessJet = -1 := by
  norm_num [lagrangian, Even.lagrangian, S11D4Odd.witness_lagrangian, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Fin.cases]

theorem odd_momentum_normalization :
    momentum ![0,0,0,2] S11D4Odd.witnessJet 1 1 = -1 := by
  rw [momentum_eq]
  have he : Even.momentum ![0,0,0,2] S11D4Odd.witnessJet 1 1 = 0 := by
    rw [Even.momentum_eq]
    change (-(0 * 0 + 0 * 0 + 0 * 1) : ℝ) = 0
    norm_num
  rw [he]
  simpa using S11D4Odd.nonzero_momentum


theorem partial_coordinate (i j : Fin 5) (x : Point 4) :
    coordDeriv i (fun y => y j) x = if j = i then 1 else 0 := by
  change deriv (fun s : ℝ => x j + s * axis i j) 0 = _
  have h := ((hasDerivAt_const (0 : ℝ) (x j)).add
    ((hasDerivAt_id (0 : ℝ)).mul_const (axis i j))).deriv
  calc
    _ = 0 + 1 * axis i j := h
    _ = _ := by simp [axis, Pi.single_apply]

def affineWitness (x : Point 4) : Vec 4 := ![x 1,x 2,0,0]

theorem affineWitness_smooth : SmoothField affineWitness := by
  intro i
  fin_cases i <;> dsimp [affineWitness] <;> fun_prop

theorem affineWitness_jet (x : Point 4) : fieldJet affineWitness x = evenJet := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [fieldJet, affineWitness, Matrix.cons_val_two, Matrix.cons_val_three,
      partial_coordinate, S11D4Odd.partial_zero, Fin.ext_iff] <;> rfl


theorem even_current_normalization (x : Point 4) :
    (∑ i : Fin 4, coordDeriv i.succ (fun y => boundaryCurrent affineWitness y i) x) = 2 := by
  rw [Even.boundary_identity affineWitness_smooth, affineWitness_jet]
  norm_num [Even.divergence, Even.transposePair, evenJet, Fin.sum_univ_four, Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Fin.cases]

theorem odd_current_normalization :
    S11D4Odd.currentJet ![0,2,0,0] S11D4Odd.witnessJet 0 = 1 :=
  S11D4Odd.current_normalization

theorem transverse_response :
    eulerLagrange ![2,3,7,11] (planeWave 0 ![1,0,0,0] ![0,1,0,0]) 0 1 = -7 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, Even.modalOperator, normSq, dot, phase, Fin.sum_univ_four,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Fin.cases]

theorem longitudinal_response :
    eulerLagrange ![2,3,7,11] (planeWave 0 ![1,0,0,0] ![1,0,0,0]) 0 0 = -12 := by
  rw [eulerLagrange_planeWave]
  norm_num [modalOperator, Even.modalOperator, normSq, dot, phase, Fin.sum_univ_four,
    Matrix.cons_val, Matrix.cons_val_two, Matrix.cons_val_three, Fin.cases]

theorem mixed_null_positive (t beta : ℝ) : VariationallyNull ![t,-t,0,beta] := by
  rw [variationallyNull_iff]
  simp [Matrix.cons_val_two]

theorem mixed_null_firstVariation {u h : Point 4 → Vec 4}
    (hu : SmoothField u) (hh : TestField h) (t beta : ℝ) :
    deriv (relativeAction ![t,-t,0,beta] u h) 0 = 0 :=
  mixed_null_positive t beta u hu h hh

end
end S11D4Bulk
