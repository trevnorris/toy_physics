import S11D3Invariants.Rotation

/-! Three explicit admissible rotations give necessary constraints only.
Every retained equation is derived from full SO invariance. Rational linear
combinations are checked by Lean, including the exhaustive representation. -/
namespace S11D3Invariants
noncomputable section
set_option maxRecDepth 2048
set_option maxHeartbeats 3200000

theorem polynomial_vec (c : Coefficients) (a0 a1 a2 a3 a4 a5 a6 a7 a8 : ℝ) :
    polynomial c ![a0,a1,a2,a3,a4,a5,a6,a7,a8] = c 0 * a0 * a0 + c 1 * a0 * a1 + c 2 * a0 * a2 + c 3 * a0 * a3 + c 4 * a0 * a4 + c 5 * a0 * a5 + c 6 * a0 * a6 + c 7 * a0 * a7 + c 8 * a0 * a8 + c 9 * a1 * a1 + c 10 * a1 * a2 + c 11 * a1 * a3 + c 12 * a1 * a4 + c 13 * a1 * a5 + c 14 * a1 * a6 + c 15 * a1 * a7 + c 16 * a1 * a8 + c 17 * a2 * a2 + c 18 * a2 * a3 + c 19 * a2 * a4 + c 20 * a2 * a5 + c 21 * a2 * a6 + c 22 * a2 * a7 + c 23 * a2 * a8 + c 24 * a3 * a3 + c 25 * a3 * a4 + c 26 * a3 * a5 + c 27 * a3 * a6 + c 28 * a3 * a7 + c 29 * a3 * a8 + c 30 * a4 * a4 + c 31 * a4 * a5 + c 32 * a4 * a6 + c 33 * a4 * a7 + c 34 * a4 * a8 + c 35 * a5 * a5 + c 36 * a5 * a6 + c 37 * a5 * a7 + c 38 * a5 * a8 + c 39 * a6 * a6 + c 40 * a6 * a7 + c 41 * a6 * a8 + c 42 * a7 * a7 + c 43 * a7 * a8 + c 44 * a8 * a8 := by rfl

theorem invariant_polynomial {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ∀ G, Q G = (c 4/2) * G.trace ^ 2 +
      (c 11/2) * (G*G).trace + c 9 * (G*G.transpose).trace := by
  have e0 : (-1) * c 0 + (1) * c 30 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![0,1,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,1,0,0,0,0] = polynomial c ![1,0,0,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e1 : (-1) * c 0 + (-1) * c 1 + (-1) * c 9 + (1) * c 24 + (-1) * c 25 + (1) * c 30 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,1,0],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![-1,1,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,1,0],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,1,0,0,0,0] = polynomial c ![1,1,0,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e2 : (-1) * c 0 + (-1) * c 2 + (-1) * c 17 + (1) * c 30 + (1) * c 31 + (1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,1],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![0,1,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,1],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,1,1,0,0,0] = polynomial c ![1,0,1,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e3 : (-1) * c 0 + (-1) * c 3 + (1) * c 9 + (-1) * c 12 + (-1) * c 24 + (1) * c 30 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![1,0,0],![0,0,0]])) = (Matrix.of ![![0,-1,0],![0,1,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![1,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,-1,0,0,1,0,0,0,0] = polynomial c ![1,0,0,1,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e4 : (-1) * c 0 + (-1) * c 5 + (1) * c 17 + (-1) * c 19 + (1) * c 30 + (-1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![0,0,1],![0,0,0]])) = (Matrix.of ![![0,0,-1],![0,1,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,-1,0,1,0,0,0,0] = polynomial c ![1,0,0,0,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e5 : (-1) * c 0 + (-1) * c 6 + (1) * c 30 + (1) * c 33 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![0,0,0],![1,0,0]])) = (Matrix.of ![![0,0,0],![0,1,0],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,0],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,1,0,0,1,0] = polynomial c ![1,0,0,0,0,0,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e6 : (-1) * c 0 + (-1) * c 7 + (1) * c 30 + (-1) * c 32 + (1) * c 39 + (-1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![0,0,0],![0,1,0]])) = (Matrix.of ![![0,0,0],![0,1,0],![-1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,0],![0,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,1,0,-1,0,0] = polynomial c ![1,0,0,0,0,0,0,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e7 : (-1) * c 0 + (-1) * c 8 + (1) * c 30 + (1) * c 34 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,1]])) = (Matrix.of ![![0,0,0],![0,1,0],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,1,0,0,0,1] = polynomial c ![1,0,0,0,0,0,0,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e8 : (-1) * c 9 + (1) * c 24 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![-1,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,0,0,0,0,0] = polynomial c ![0,1,0,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e9 : (-1) * c 9 + (-1) * c 10 + (-1) * c 17 + (1) * c 24 + (-1) * c 26 + (1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,1],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![-1,0,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,1],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,0,1,0,0,0] = polynomial c ![0,1,1,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e10 : (-1) * c 9 + (-1) * c 13 + (1) * c 17 + (1) * c 18 + (1) * c 24 + (-1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,0],![0,0,1],![0,0,0]])) = (Matrix.of ![![0,0,-1],![-1,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,-1,-1,0,0,0,0,0] = polynomial c ![0,1,0,0,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e11 : (-1) * c 9 + (-1) * c 14 + (1) * c 24 + (-1) * c 28 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,0],![0,0,0],![1,0,0]])) = (Matrix.of ![![0,0,0],![-1,0,0],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,0],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,0,0,0,1,0] = polynomial c ![0,1,0,0,0,0,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e12 : (-1) * c 9 + (-1) * c 15 + (1) * c 24 + (1) * c 27 + (1) * c 39 + (-1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,0],![0,0,0],![0,1,0]])) = (Matrix.of ![![0,0,0],![-1,0,0],![-1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,0],![0,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,0,0,-1,0,0] = polynomial c ![0,1,0,0,0,0,0,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e13 : (-1) * c 9 + (-1) * c 16 + (1) * c 24 + (-1) * c 29 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,1]])) = (Matrix.of ![![0,0,0],![-1,0,0],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,-1,0,0,0,0,1] = polynomial c ![0,1,0,0,0,0,0,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e14 : (-1) * c 17 + (1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,0,0],![0,0,0]])) = (Matrix.of ![![0,0,0],![0,0,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,1,0,0,0] = polynomial c ![0,0,1,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e15 : (1) * c 9 + (-1) * c 13 + (-1) * c 17 + (-1) * c 18 + (-1) * c 24 + (1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![1,0,0],![0,0,0]])) = (Matrix.of ![![0,-1,0],![0,0,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![1,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,-1,0,0,0,1,0,0,0] = polynomial c ![0,0,1,1,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e16 : (1) * c 0 + (1) * c 5 + (-1) * c 17 + (-1) * c 19 + (-1) * c 30 + (1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,1,0],![0,0,0]])) = (Matrix.of ![![1,0,0],![0,0,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,1,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,1,0,0,0] = polynomial c ![0,0,1,0,1,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e17 : (-2) * c 20 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,0,1],![0,0,0]])) = (Matrix.of ![![0,0,-1],![0,0,1],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,-1,0,0,1,0,0,0] = polynomial c ![0,0,1,0,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e18 : (-1) * c 17 + (-1) * c 21 + (1) * c 35 + (1) * c 37 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,0,0],![1,0,0]])) = (Matrix.of ![![0,0,0],![0,0,1],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,0,0],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,1,0,1,0] = polynomial c ![0,0,1,0,0,0,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e19 : (-1) * c 17 + (-1) * c 22 + (1) * c 35 + (-1) * c 36 + (1) * c 39 + (-1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,0,0],![0,1,0]])) = (Matrix.of ![![0,0,0],![0,0,1],![-1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,0,0],![0,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,1,-1,0,0] = polynomial c ![0,0,1,0,0,0,0,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e20 : (-1) * c 17 + (-1) * c 23 + (1) * c 35 + (1) * c 38 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,1],![0,0,0],![0,0,1]])) = (Matrix.of ![![0,0,0],![0,0,1],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,1],![0,0,0],![0,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,1,0,0,1] = polynomial c ![0,0,1,0,0,0,0,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e21 : (1) * c 9 + (1) * c 10 + (1) * c 17 + (-1) * c 24 + (-1) * c 26 + (-1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![1,0,1],![0,0,0]])) = (Matrix.of ![![0,-1,-1],![0,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![1,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,-1,-1,0,0,0,0,0,0] = polynomial c ![0,0,0,1,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e22 : (1) * c 9 + (-1) * c 15 + (-1) * c 24 + (-1) * c 27 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![1,0,0],![1,0,0]])) = (Matrix.of ![![0,-1,0],![0,0,0],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![1,0,0],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,-1,0,0,0,0,0,1,0] = polynomial c ![0,0,0,1,0,0,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e23 : (1) * c 9 + (1) * c 14 + (-1) * c 24 + (-1) * c 28 + (1) * c 39 + (-1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![1,0,0],![0,1,0]])) = (Matrix.of ![![0,-1,0],![0,0,0],![-1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![1,0,0],![0,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,-1,0,0,0,0,-1,0,0] = polynomial c ![0,0,0,1,0,0,0,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e24 : (1) * c 0 + (-1) * c 2 + (1) * c 17 + (-1) * c 30 + (-1) * c 31 + (-1) * c 35 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,1,1],![0,0,0]])) = (Matrix.of ![![1,0,-1],![0,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,1,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,-1,0,0,0,0,0,0] = polynomial c ![0,0,0,0,1,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e25 : (1) * c 0 + (1) * c 7 + (-1) * c 30 + (-1) * c 32 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,1,0],![1,0,0]])) = (Matrix.of ![![1,0,0],![0,0,0],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,1,0],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,0,0,1,0] = polynomial c ![0,0,0,0,1,0,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e26 : (1) * c 0 + (-1) * c 6 + (-1) * c 30 + (-1) * c 33 + (1) * c 39 + (-1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,1,0],![0,1,0]])) = (Matrix.of ![![1,0,0],![0,0,0],![-1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,1,0],![0,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,0,-1,0,0] = polynomial c ![0,0,0,0,1,0,0,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e27 : (1) * c 17 + (-1) * c 22 + (-1) * c 35 + (-1) * c 36 + (-1) * c 39 + (1) * c 42 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,0,1],![1,0,0]])) = (Matrix.of ![![0,0,-1],![0,0,0],![0,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,0,1],![1,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,-1,0,0,0,0,1,0] = polynomial c ![0,0,0,0,0,1,1,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e28 : (1) * c 17 + (-1) * c 23 + (-1) * c 35 + (-1) * c 38 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,0,1],![0,0,1]])) = (Matrix.of ![![0,0,-1],![0,0,0],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,0,1],![0,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,-1,0,0,0,0,0,1] = polynomial c ![0,0,0,0,0,1,0,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e29 : (-2) * c 40 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,0,0],![1,1,0]])) = (Matrix.of ![![0,0,0],![0,0,0],![-1,1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,0,0],![1,1,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,0,-1,1,0] = polynomial c ![0,0,0,0,0,0,1,1,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e30 : (-1) * c 39 + (-1) * c 41 + (1) * c 42 + (1) * c 43 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,0,0],![1,0,1]])) = (Matrix.of ![![0,0,0],![0,0,0],![0,1,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,0,0],![1,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,0,0,1,1] = polynomial c ![0,0,0,0,0,0,1,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e31 : (1) * c 39 + (-1) * c 41 + (-1) * c 42 + (-1) * c 43 = 0 := by
    have image : conjugate (rotationXY 0 1) ((Matrix.of ![![0,0,0],![0,0,0],![0,1,1]])) = (Matrix.of ![![0,0,0],![0,0,0],![-1,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY 0 1) (rotationXY_proper (by norm_num)) ((Matrix.of ![![0,0,0],![0,0,0],![0,1,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,0,0,0,0,-1,0,1] = polynomial c ![0,0,0,0,0,0,0,1,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e32 : (-1) * c 1 + (1) * c 2 + (-1) * c 9 + (1) * c 17 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![1,1,0],![0,0,0],![0,0,0]])) = (Matrix.of ![![1,0,1],![0,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![1,1,0],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,1,0,0,0,0,0,0] = polynomial c ![1,1,0,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e33 : (-1) * c 1 + (-1) * c 2 + (1) * c 9 + (-1) * c 17 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![1,0,1],![0,0,0],![0,0,0]])) = (Matrix.of ![![1,-1,0],![0,0,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![1,0,1],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,-1,0,0,0,0,0,0,0] = polynomial c ![1,0,1,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e34 : (-1) * c 3 + (1) * c 6 + (-1) * c 24 + (1) * c 39 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![1,0,0],![1,0,0],![0,0,0]])) = (Matrix.of ![![1,0,0],![0,0,0],![1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![1,0,0],![1,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,0,1,0,0] = polynomial c ![1,0,0,1,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e35 : (-1) * c 4 + (1) * c 8 + (-1) * c 30 + (1) * c 44 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![1,0,0],![0,1,0],![0,0,0]])) = (Matrix.of ![![1,0,0],![0,0,0],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,1,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,0,0,0,1] = polynomial c ![1,0,0,0,1,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e36 : (-1) * c 5 + (-1) * c 7 + (-1) * c 35 + (1) * c 42 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![1,0,0],![0,0,1],![0,0,0]])) = (Matrix.of ![![1,0,0],![0,0,0],![0,-1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![1,0,0,0,0,0,0,-1,0] = polynomial c ![1,0,0,0,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e37 : (-1) * c 9 + (-1) * c 11 + (1) * c 17 + (1) * c 21 + (-1) * c 24 + (1) * c 39 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![0,1,0],![1,0,0],![0,0,0]])) = (Matrix.of ![![0,0,1],![0,0,0],![1,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![0,1,0],![1,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,1,0,0,0,1,0,0] = polynomial c ![0,1,0,1,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e38 : (-1) * c 9 + (-1) * c 12 + (1) * c 17 + (1) * c 23 + (-1) * c 30 + (1) * c 44 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![0,1,0],![0,1,0],![0,0,0]])) = (Matrix.of ![![0,0,1],![0,0,0],![0,0,1]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,1,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,1,0,0,0,0,0,1] = polynomial c ![0,1,0,0,1,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e39 : (-1) * c 9 + (-1) * c 13 + (1) * c 17 + (-1) * c 22 + (-1) * c 35 + (1) * c 42 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![0,1,0],![0,0,1],![0,0,0]])) = (Matrix.of ![![0,0,1],![0,0,0],![0,-1,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,1],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,1,0,0,0,0,-1,0] = polynomial c ![0,1,0,0,0,1,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e40 : (-1) * c 9 + (-1) * c 16 + (1) * c 17 + (1) * c 19 + (1) * c 30 + (-1) * c 44 = 0 := by
    have image : conjugate (rotationYZ 0 1) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,1]])) = (Matrix.of ![![0,0,1],![0,1,0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationYZ, Matrix.cons_val_two]
    have raw := hQ (rotationYZ 0 1) (rotationYZ_proper (by norm_num)) ((Matrix.of ![![0,1,0],![0,0,0],![0,0,1]]))
    rw [image, hc, hc] at raw
    change polynomial c ![0,0,1,0,1,0,0,0,0] = polynomial c ![0,1,0,0,0,0,0,0,1] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have e41 : ((-544/625)) * c 0 + ((108/625)) * c 1 + ((108/625)) * c 3 + ((144/625)) * c 4 + ((144/625)) * c 9 + ((144/625)) * c 11 + ((192/625)) * c 12 + ((144/625)) * c 24 + ((192/625)) * c 25 + ((256/625)) * c 30 = 0 := by
    have image : conjugate (rotationXY (3/5) (4/5)) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,0]])) = (Matrix.of ![![(9/25),(12/25),0],![(12/25),(16/25),0],![0,0,0]]) := by
      ext i j
      rw [conjugate_apply]
      fin_cases i <;> fin_cases j <;>
        norm_num [Fin.sum_univ_three, rotationXY, Matrix.cons_val_two]
    have raw := hQ (rotationXY (3/5) (4/5)) (rotationXY_proper (by norm_num)) ((Matrix.of ![![1,0,0],![0,0,0],![0,0,0]]))
    rw [image, hc, hc] at raw
    change polynomial c ![(9/25),(12/25),0,(12/25),(16/25),0,0,0,0] = polynomial c ![1,0,0,0,0,0,0,0,0] at raw
    rw [polynomial_vec, polynomial_vec] at raw
    norm_num at raw
    linarith only [raw]
  have h0 : c 0 = 1 * (c 4/2) + 1 * (c 11/2) + 1 * c 9 := by
    linear_combination ((59/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e3 + ((7/48)) * e4 + ((7/48)) * e5 + ((7/48)) * e6 + ((19/24)) * e8 + ((7/12)) * e14 + ((-7/48)) * e16 + ((-7/24)) * e19 + ((-7/48)) * e24 + ((-7/48)) * e25 + ((7/48)) * e26 + ((7/24)) * e27 + ((7/24)) * e33 + ((7/24)) * e34 + ((-7/24)) * e36 + ((-625/288)) * e41
  have h1 : c 1 = 0 := by
    linear_combination ((-1/2)) * e32 + ((-1/2)) * e33
  have h2 : c 2 = 0 := by
    linear_combination ((-1/2)) * e2 + ((-1/2)) * e24
  have h3 : c 3 = 0 := by
    linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e5 + ((-1/2)) * e6 + (-1) * e8 + (-2) * e14 + ((1/2)) * e16 + (1) * e19 + ((1/2)) * e24 + ((1/2)) * e25 + ((-1/2)) * e26 + (-1) * e27 + ((1/2)) * e32 + ((-1/2)) * e33 + (-1) * e34 + (1) * e36
  have h5 : c 5 = 0 := by
    linear_combination (1) * e0 + ((-1/2)) * e4 + (-1) * e14 + ((1/2)) * e16
  have h6 : c 6 = 0 := by
    linear_combination ((-1/2)) * e5 + ((-1/2)) * e26
  have h7 : c 7 = 0 := by
    linear_combination (1) * e0 + ((-1/2)) * e6 + (-1) * e14 + ((1/2)) * e19 + ((1/2)) * e25 + ((-1/2)) * e27
  have h8 : c 8 = 2 * (c 4/2) := by
    linear_combination (1) * e0 + (1) * e2 + (1) * e3 + ((-1/2)) * e4 + ((-1/2)) * e5 + ((-1/2)) * e6 + (-2) * e14 + ((1/2)) * e16 + (1) * e19 + ((-1/2)) * e20 + (1) * e24 + ((1/2)) * e25 + ((-1/2)) * e26 + (-1) * e27 + ((-1/2)) * e28 + (1) * e32 + (-1) * e33 + (-1) * e34 + (1) * e35 + (1) * e36 + (-1) * e38
  have h10 : c 10 = 0 := by
    linear_combination (1) * e8 + ((-1/2)) * e9 + (1) * e14 + ((1/2)) * e21
  have h12 : c 12 = 0 := by
    linear_combination (-1) * e0 + ((-1/2)) * e2 + (-1) * e3 + ((1/2)) * e4 + ((1/2)) * e5 + ((1/2)) * e6 + (2) * e14 + ((-1/2)) * e16 + (-1) * e19 + ((-1/2)) * e24 + ((-1/2)) * e25 + ((1/2)) * e26 + (1) * e27 + ((-1/2)) * e32 + ((1/2)) * e33 + (1) * e34 + (-1) * e36
  have h13 : c 13 = 0 := by
    linear_combination ((-1/2)) * e10 + ((-1/2)) * e15
  have h14 : c 14 = 0 := by
    linear_combination (1) * e8 + ((-1/2)) * e11 + (1) * e14 + ((-1/2)) * e19 + ((1/2)) * e23 + ((1/2)) * e27
  have h15 : c 15 = 0 := by
    linear_combination ((-1/2)) * e12 + ((-1/2)) * e22
  have h16 : c 16 = 0 := by
    linear_combination (1) * e0 + ((3/2)) * e2 + (1) * e3 + (-1) * e4 + ((-1/2)) * e5 + ((-1/2)) * e6 + (-2) * e14 + (1) * e19 + ((-1/2)) * e20 + ((3/2)) * e24 + ((1/2)) * e25 + ((-1/2)) * e26 + (-1) * e27 + ((-1/2)) * e28 + ((3/2)) * e32 + ((-3/2)) * e33 + (-1) * e34 + (1) * e36 + (-1) * e38 + (-1) * e40
  have h17 : c 17 = 1 * c 9 := by
    linear_combination ((1/2)) * e2 + ((1/2)) * e24 + ((1/2)) * e32 + ((-1/2)) * e33
  have h18 : c 18 = 0 := by
    linear_combination (-1) * e8 + ((1/2)) * e10 + (1) * e14 + ((-1/2)) * e15
  have h19 : c 19 = 0 := by
    linear_combination ((-1/2)) * e4 + ((-1/2)) * e16
  have h20 : c 20 = 0 := by
    linear_combination ((-1/2)) * e17
  have h21 : c 21 = 2 * (c 11/2) := by
    linear_combination (-2) * e0 + (-1) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + (1) * e8 + (2) * e14 + ((-1/2)) * e16 + (-1) * e19 + (-1) * e24 + ((-1/2)) * e25 + (1) * e27 + (-1) * e32 + (1) * e33 + (-1) * e36 + (1) * e37
  have h22 : c 22 = 0 := by
    linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e6 + ((1/2)) * e10 + (-2) * e14 + ((1/2)) * e15 + ((1/2)) * e16 + ((1/2)) * e19 + ((1/2)) * e24 + ((1/2)) * e25 + ((-1/2)) * e27 + ((1/2)) * e32 + ((-1/2)) * e33 + (1) * e36 + (-1) * e39
  have h23 : c 23 = 0 := by
    linear_combination ((-1/2)) * e20 + ((-1/2)) * e28
  have h24 : c 24 = 1 * c 9 := by
    linear_combination (1) * e8
  have h25 : c 25 = 0 := by
    linear_combination (1) * e0 + (-1) * e1 + (1) * e8 + ((1/2)) * e32 + ((1/2)) * e33
  have h26 : c 26 = 0 := by
    linear_combination ((-1/2)) * e9 + ((-1/2)) * e21
  have h27 : c 27 = 0 := by
    linear_combination (-1) * e8 + ((1/2)) * e12 + (1) * e14 + ((-1/2)) * e19 + ((-1/2)) * e22 + ((1/2)) * e27
  have h28 : c 28 = 0 := by
    linear_combination ((-1/2)) * e11 + ((-1/2)) * e23
  have h29 : c 29 = 0 := by
    linear_combination (-1) * e0 + ((-3/2)) * e2 + (-1) * e3 + (1) * e4 + ((1/2)) * e5 + ((1/2)) * e6 + (1) * e8 + (-1) * e13 + (2) * e14 + (-1) * e19 + ((1/2)) * e20 + ((-3/2)) * e24 + ((-1/2)) * e25 + ((1/2)) * e26 + (1) * e27 + ((1/2)) * e28 + ((-3/2)) * e32 + ((3/2)) * e33 + (1) * e34 + (-1) * e36 + (1) * e38 + (1) * e40
  have h30 : c 30 = 1 * (c 4/2) + 1 * (c 11/2) + 1 * c 9 := by
    linear_combination ((95/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e3 + ((7/48)) * e4 + ((7/48)) * e5 + ((7/48)) * e6 + ((19/24)) * e8 + ((7/12)) * e14 + ((-7/48)) * e16 + ((-7/24)) * e19 + ((-7/48)) * e24 + ((-7/48)) * e25 + ((7/48)) * e26 + ((7/24)) * e27 + ((7/24)) * e33 + ((7/24)) * e34 + ((-7/24)) * e36 + ((-625/288)) * e41
  have h31 : c 31 = 0 := by
    linear_combination (-1) * e0 + ((1/2)) * e2 + (-1) * e14 + ((-1/2)) * e24
  have h32 : c 32 = 0 := by
    linear_combination ((-1/2)) * e6 + ((-1/2)) * e25
  have h33 : c 33 = 0 := by
    linear_combination (-1) * e0 + ((1/2)) * e5 + (-1) * e14 + ((1/2)) * e19 + ((-1/2)) * e26 + ((-1/2)) * e27
  have h34 : c 34 = 2 * (c 4/2) := by
    linear_combination (1) * e2 + (1) * e3 + ((-1/2)) * e4 + ((-1/2)) * e5 + ((-1/2)) * e6 + (1) * e7 + (-2) * e14 + ((1/2)) * e16 + (1) * e19 + ((-1/2)) * e20 + (1) * e24 + ((1/2)) * e25 + ((-1/2)) * e26 + (-1) * e27 + ((-1/2)) * e28 + (1) * e32 + (-1) * e33 + (-1) * e34 + (1) * e35 + (1) * e36 + (-1) * e38
  have h35 : c 35 = 1 * c 9 := by
    linear_combination ((1/2)) * e2 + (1) * e14 + ((1/2)) * e24 + ((1/2)) * e32 + ((-1/2)) * e33
  have h36 : c 36 = 0 := by
    linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + ((-1/2)) * e10 + (2) * e14 + ((-1/2)) * e15 + ((-1/2)) * e16 + (-1) * e19 + ((-1/2)) * e24 + ((-1/2)) * e25 + ((-1/2)) * e32 + ((1/2)) * e33 + (-1) * e36 + (1) * e39
  have h37 : c 37 = 2 * (c 11/2) := by
    linear_combination (-2) * e0 + (-1) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + (1) * e8 + ((-1/2)) * e16 + (1) * e18 + ((-1/2)) * e19 + (-1) * e24 + ((-1/2)) * e25 + ((1/2)) * e27 + (-1) * e32 + (1) * e33 + (-1) * e36 + (1) * e37
  have h38 : c 38 = 0 := by
    linear_combination (-1) * e14 + ((1/2)) * e20 + ((-1/2)) * e28
  have h39 : c 39 = 1 * c 9 := by
    linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e6 + (-2) * e14 + ((1/2)) * e16 + (1) * e19 + ((1/2)) * e24 + ((1/2)) * e25 + (-1) * e27 + ((1/2)) * e32 + ((-1/2)) * e33 + (1) * e36
  have h40 : c 40 = 0 := by
    linear_combination ((-1/2)) * e29
  have h41 : c 41 = 0 := by
    linear_combination ((-1/2)) * e30 + ((-1/2)) * e31
  have h42 : c 42 = 1 * c 9 := by
    linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e6 + (-1) * e14 + ((1/2)) * e16 + ((1/2)) * e19 + ((1/2)) * e24 + ((1/2)) * e25 + ((-1/2)) * e27 + ((1/2)) * e32 + ((-1/2)) * e33 + (1) * e36
  have h43 : c 43 = 0 := by
    linear_combination (-1) * e14 + ((1/2)) * e19 + ((-1/2)) * e27 + ((1/2)) * e30 + ((-1/2)) * e31
  have h44 : c 44 = 1 * (c 4/2) + 1 * (c 11/2) + 1 * c 9 := by
    linear_combination ((59/36)) * e0 + ((-2/3)) * e1 + ((-55/48)) * e2 + ((-5/3)) * e3 + ((31/48)) * e4 + ((31/48)) * e5 + ((31/48)) * e6 + ((19/24)) * e8 + ((31/12)) * e14 + ((-31/48)) * e16 + ((-31/24)) * e19 + ((1/2)) * e20 + ((-55/48)) * e24 + ((-31/48)) * e25 + ((31/48)) * e26 + ((31/24)) * e27 + ((1/2)) * e28 + (-1) * e32 + ((31/24)) * e33 + ((31/24)) * e34 + ((-31/24)) * e36 + (1) * e38 + ((-625/288)) * e41
  intro G
  rw [hc]
  change polynomial c ![G 0 0,G 0 1,G 0 2,G 1 0,G 1 1,G 1 2,G 2 0,G 2 1,G 2 2] = _
  rw [polynomial_vec]
  simp only [h0, h1, h2, h3, h5, h6, h7, h8, h10, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h25, h26, h27, h28, h29, h30, h31, h32, h33, h34, h35, h36, h37, h38, h39, h40, h41, h42, h43, h44]
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]
  ring

end
end S11D3Invariants
