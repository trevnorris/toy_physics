import Mathlib.LinearAlgebra.QuadraticForm.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic

/-! Exhaustive quadratic forms on all real 3x3 matrices. Generated finite
coordinate algebra, not CAS-output transcription; see S11_lean_d3_generate.py. -/
namespace S11D3Invariants
noncomputable section
set_option maxRecDepth 2048
set_option maxHeartbeats 1600000

abbrev Mat := Matrix (Fin 3) (Fin 3) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 9 → ℝ
abbrev Coefficients := Fin 45 → ℝ

def coordinates (G : Mat) : Coordinates := ![G 0 0,G 0 1,G 0 2,G 1 0,G 1 1,G 1 2,G 2 0,G 2 1,G 2 2]
def decode (v : Coordinates) : Mat := ![![v 0,v 1,v 2],![v 3,v 4,v 5],![v 6,v 7,v 8]]

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode]

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 9) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates]
  map_smul' r G := by fin_cases i <;> simp [coordinates]

def frame : Fin 9 → Mat := ![![![1,0,0],![0,0,0],![0,0,0]],![![0,1,0],![0,0,0],![0,0,0]],![![0,0,1],![0,0,0],![0,0,0]],![![0,0,0],![1,0,0],![0,0,0]],![![0,0,0],![0,1,0],![0,0,0]],![![0,0,0],![0,0,1],![0,0,0]],![![0,0,0],![0,0,0],![1,0,0]],![![0,0,0],![0,0,0],![0,1,0]],![![0,0,0],![0,0,0],![0,0,1]]]

theorem frame_expansion (G : Mat) :
    G = coordinates G 0 • frame 0 + coordinates G 1 • frame 1 + coordinates G 2 • frame 2 + coordinates G 3 • frame 3 + coordinates G 4 • frame 4 + coordinates G 5 • frame 5 + coordinates G 6 • frame 6 + coordinates G 7 • frame 7 + coordinates G 8 • frame 8 := by
  ext i j
  fin_cases i <;> fin_cases j
  · change G 0 0 = G 0 0 * 1 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 0 1 = G 0 0 * 0 + G 0 1 * 1 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 0 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 1 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 1 0 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 1 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 1 1 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 1 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 1 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 1 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 2 0 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 1 + G 2 1 * 0 + G 2 2 * 0
    ring
  · change G 2 1 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 1 + G 2 2 * 0
    ring
  · change G 2 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 1
    ring

def polynomial (c : Coefficients) (v : Coordinates) : ℝ :=
    c 0 * v 0 * v 0 +
    c 1 * v 0 * v 1 +
    c 2 * v 0 * v 2 +
    c 3 * v 0 * v 3 +
    c 4 * v 0 * v 4 +
    c 5 * v 0 * v 5 +
    c 6 * v 0 * v 6 +
    c 7 * v 0 * v 7 +
    c 8 * v 0 * v 8 +
    c 9 * v 1 * v 1 +
    c 10 * v 1 * v 2 +
    c 11 * v 1 * v 3 +
    c 12 * v 1 * v 4 +
    c 13 * v 1 * v 5 +
    c 14 * v 1 * v 6 +
    c 15 * v 1 * v 7 +
    c 16 * v 1 * v 8 +
    c 17 * v 2 * v 2 +
    c 18 * v 2 * v 3 +
    c 19 * v 2 * v 4 +
    c 20 * v 2 * v 5 +
    c 21 * v 2 * v 6 +
    c 22 * v 2 * v 7 +
    c 23 * v 2 * v 8 +
    c 24 * v 3 * v 3 +
    c 25 * v 3 * v 4 +
    c 26 * v 3 * v 5 +
    c 27 * v 3 * v 6 +
    c 28 * v 3 * v 7 +
    c 29 * v 3 * v 8 +
    c 30 * v 4 * v 4 +
    c 31 * v 4 * v 5 +
    c 32 * v 4 * v 6 +
    c 33 * v 4 * v 7 +
    c 34 * v 4 * v 8 +
    c 35 * v 5 * v 5 +
    c 36 * v 5 * v 6 +
    c 37 * v 5 * v 7 +
    c 38 * v 5 * v 8 +
    c 39 * v 6 * v 6 +
    c 40 * v 6 * v 7 +
    c 41 * v 6 * v 8 +
    c 42 * v 7 * v 7 +
    c 43 * v 7 * v 8 +
    c 44 * v 8 * v 8

def monomial (i j : Fin 9) : Quad :=
  QuadraticMap.linMulLin (coordinateMap i) (coordinateMap j)

theorem quadratic_representation (Q : Quad) :
    ∃ c : Coefficients, ∀ G, Q G = polynomial c (coordinates G) := by
  let b := Q.associated
  let c : Coefficients := ![b (frame 0) (frame 0),b (frame 0) (frame 1) + b (frame 1) (frame 0),b (frame 0) (frame 2) + b (frame 2) (frame 0),b (frame 0) (frame 3) + b (frame 3) (frame 0),b (frame 0) (frame 4) + b (frame 4) (frame 0),b (frame 0) (frame 5) + b (frame 5) (frame 0),b (frame 0) (frame 6) + b (frame 6) (frame 0),b (frame 0) (frame 7) + b (frame 7) (frame 0),b (frame 0) (frame 8) + b (frame 8) (frame 0),b (frame 1) (frame 1),b (frame 1) (frame 2) + b (frame 2) (frame 1),b (frame 1) (frame 3) + b (frame 3) (frame 1),b (frame 1) (frame 4) + b (frame 4) (frame 1),b (frame 1) (frame 5) + b (frame 5) (frame 1),b (frame 1) (frame 6) + b (frame 6) (frame 1),b (frame 1) (frame 7) + b (frame 7) (frame 1),b (frame 1) (frame 8) + b (frame 8) (frame 1),b (frame 2) (frame 2),b (frame 2) (frame 3) + b (frame 3) (frame 2),b (frame 2) (frame 4) + b (frame 4) (frame 2),b (frame 2) (frame 5) + b (frame 5) (frame 2),b (frame 2) (frame 6) + b (frame 6) (frame 2),b (frame 2) (frame 7) + b (frame 7) (frame 2),b (frame 2) (frame 8) + b (frame 8) (frame 2),b (frame 3) (frame 3),b (frame 3) (frame 4) + b (frame 4) (frame 3),b (frame 3) (frame 5) + b (frame 5) (frame 3),b (frame 3) (frame 6) + b (frame 6) (frame 3),b (frame 3) (frame 7) + b (frame 7) (frame 3),b (frame 3) (frame 8) + b (frame 8) (frame 3),b (frame 4) (frame 4),b (frame 4) (frame 5) + b (frame 5) (frame 4),b (frame 4) (frame 6) + b (frame 6) (frame 4),b (frame 4) (frame 7) + b (frame 7) (frame 4),b (frame 4) (frame 8) + b (frame 8) (frame 4),b (frame 5) (frame 5),b (frame 5) (frame 6) + b (frame 6) (frame 5),b (frame 5) (frame 7) + b (frame 7) (frame 5),b (frame 5) (frame 8) + b (frame 8) (frame 5),b (frame 6) (frame 6),b (frame 6) (frame 7) + b (frame 7) (frame 6),b (frame 6) (frame 8) + b (frame 8) (frame 6),b (frame 7) (frame 7),b (frame 7) (frame 8) + b (frame 8) (frame 7),b (frame 8) (frame 8)]
  refine ⟨c, fun G => ?_⟩
  have hb : Q G = b G G := (Q.associated_eq_self_apply ℝ G).symm
  rw [hb, congrArg (fun Z => b Z Z) (frame_expansion G)]
  simp only [polynomial, c, Matrix.cons_val]
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
  ring

end
end S11D3Invariants
