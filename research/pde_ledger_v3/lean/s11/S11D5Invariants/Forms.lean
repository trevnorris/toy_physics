import S11D5Invariants.Rotation

/-! Three fully O(5)-invariant quadratic forms, with explicit normalization. -/
namespace S11D5Invariants
noncomputable section

def traceSquare : Quad :=
  monomial 0 0 + monomial 6 6 + monomial 12 12 + monomial 18 18 + monomial 24 24 + 2 • monomial 0 6 + 2 • monomial 0 12 + 2 • monomial 0 18 + 2 • monomial 0 24 + 2 • monomial 6 12 + 2 • monomial 6 18 + 2 • monomial 6 24 + 2 • monomial 12 18 + 2 • monomial 12 24 + 2 • monomial 18 24

def traceOfSquare : Quad :=
  monomial 0 0 + monomial 6 6 + monomial 12 12 + monomial 18 18 + monomial 24 24 + 2 • monomial 1 5 + 2 • monomial 2 10 + 2 • monomial 3 15 + 2 • monomial 4 20 + 2 • monomial 7 11 + 2 • monomial 8 16 + 2 • monomial 9 21 + 2 • monomial 13 17 + 2 • monomial 14 22 + 2 • monomial 19 23

def frobeniusSquare : Quad :=
  monomial 0 0 + monomial 1 1 + monomial 2 2 + monomial 3 3 + monomial 4 4 + monomial 5 5 + monomial 6 6 + monomial 7 7 + monomial 8 8 + monomial 9 9 + monomial 10 10 + monomial 11 11 + monomial 12 12 + monomial 13 13 + monomial 14 14 + monomial 15 15 + monomial 16 16 + monomial 17 17 + monomial 18 18 + monomial 19 19 + monomial 20 20 + monomial 21 21 + monomial 22 22 + monomial 23 23 + monomial 24 24

theorem traceSquare_apply (G : Mat) : traceSquare G = G.trace ^ 2 := by
  simp [traceSquare, monomial, coordinateMap, coordinates, Matrix.trace, sum_five]
  ring

theorem traceOfSquare_apply (G : Mat) : traceOfSquare G = (G*G).trace := by
  simp [traceOfSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, sum_five]
  ring

theorem frobeniusSquare_apply (G : Mat) : frobeniusSquare G = (G*G.transpose).trace := by
  simp [frobeniusSquare, monomial, coordinateMap, coordinates, Matrix.trace,
    Matrix.mul_apply, sum_five]
  ring

def invariantForm (v : Fin 3 → ℝ) : Quad :=
  v 0 • traceSquare + v 1 • traceOfSquare + v 2 • frobeniusSquare

theorem invariantForm_apply (v : Fin 3 → ℝ) (G : Mat) :
    invariantForm v G = v 0 * G.trace ^ 2 + v 1 * (G*G).trace +
      v 2 * (G*G.transpose).trace := by
  simp [invariantForm, traceSquare_apply, traceOfSquare_apply, frobeniusSquare_apply]

theorem invariantForm_O (v : Fin 3 → ℝ) : OInvariant (invariantForm v) := by
  intro R hR G
  simp only [invariantForm_apply, trace_conjugate hR, conjugate_transpose,
    conjugate_mul hR]

theorem invariantForm_SO (v : Fin 3 → ℝ) : SOInvariant (invariantForm v) :=
  fun R hR => invariantForm_O v R hR.1

end
end S11D5Invariants
