import S11D4Invariants.Orientation

/-! Four explicit quadratic forms and full-group sufficiency. -/
namespace S11D4Invariants
noncomputable section

def traceSquare : Quad :=
  monomial 0 0 + monomial 5 5 + monomial 10 10 + monomial 15 15 + 2 • monomial 0 5 + 2 • monomial 0 10 + 2 • monomial 0 15 + 2 • monomial 5 10 + 2 • monomial 5 15 + 2 • monomial 10 15

def traceOfSquare : Quad :=
  monomial 0 0 + monomial 5 5 + monomial 10 10 + monomial 15 15 + 2 • monomial 1 4 + 2 • monomial 2 8 + 2 • monomial 3 12 + 2 • monomial 6 9 + 2 • monomial 7 13 + 2 • monomial 11 14

def frobeniusSquare : Quad :=
  monomial 0 0 + monomial 1 1 + monomial 2 2 + monomial 3 3 + monomial 4 4 + monomial 5 5 + monomial 6 6 + monomial 7 7 + monomial 8 8 + monomial 9 9 + monomial 10 10 + monomial 11 11 + monomial 12 12 + monomial 13 13 + monomial 14 14 + monomial 15 15

def orientationForm : Quad :=
  monomial 1 11 - monomial 1 14 - monomial 4 11 + monomial 4 14 -
  monomial 2 7 + monomial 2 13 + monomial 8 7 - monomial 8 13 +
  monomial 3 6 - monomial 3 9 - monomial 12 6 + monomial 12 9

theorem traceSquare_apply (G : Mat) : traceSquare G = G.trace ^ 2 := by
  simp [traceSquare, monomial, coordinateMap, coordinates, Matrix.trace, Fin.sum_univ_four]
  ring

theorem traceOfSquare_apply (G : Mat) : traceOfSquare G = (G*G).trace := by
  simp [traceOfSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_four]
  ring

theorem frobeniusSquare_apply (G : Mat) : frobeniusSquare G = (G*G.transpose).trace := by
  simp [frobeniusSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, Fin.sum_univ_four]
  ring

theorem orientationForm_apply (G : Mat) : orientationForm G = orientation G := by
  simp [orientationForm, monomial, coordinateMap, coordinates, orientation]
  ring

def invariantForm (v : Fin 4 → ℝ) : Quad :=
  v 0 • traceSquare + v 1 • traceOfSquare + v 2 • frobeniusSquare + v 3 • orientationForm

theorem invariantForm_apply (v : Fin 4 → ℝ) (G : Mat) :
    invariantForm v G = v 0 * G.trace ^ 2 + v 1 * (G*G).trace +
      v 2 * (G*G.transpose).trace + v 3 * orientation G := by
  simp [invariantForm, traceSquare_apply, traceOfSquare_apply, frobeniusSquare_apply,
    orientationForm_apply]

theorem invariantForm_SO (v : Fin 4 → ℝ) : SOInvariant (invariantForm v) := by
  intro R hR G
  simp only [invariantForm_apply, trace_conjugate hR.1, conjugate_transpose,
    conjugate_mul hR.1, orientation_conjugate, hR.2, one_mul]

theorem invariantForm_reflection (v : Fin 4 → ℝ) (G : Mat) :
    invariantForm v (conjugate reflection G) =
      v 0 * G.trace ^ 2 + v 1 * (G*G).trace +
      v 2 * (G*G.transpose).trace - v 3 * orientation G := by
  simp only [invariantForm_apply, trace_conjugate reflection_orthogonal,
    conjugate_transpose, conjugate_mul reflection_orthogonal,
    orientation_conjugate, reflection_det]
  ring

end
end S11D4Invariants
