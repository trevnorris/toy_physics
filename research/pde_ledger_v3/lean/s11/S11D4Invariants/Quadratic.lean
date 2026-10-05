import Mathlib.LinearAlgebra.QuadraticForm.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic

/-! Exhaustive quadratic forms on all real 4x4 matrices. Generated finite
coordinate algebra, not native output transcription; see S11_lean_d4_generate.py. -/
namespace S11D4Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 1600000

abbrev Mat := Matrix (Fin 4) (Fin 4) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 16 → ℝ
abbrev Coefficients := Fin 136 → ℝ

def coordinates (G : Mat) : Coordinates := ![G 0 0,G 0 1,G 0 2,G 0 3,G 1 0,G 1 1,G 1 2,G 1 3,G 2 0,G 2 1,G 2 2,G 2 3,G 3 0,G 3 1,G 3 2,G 3 3]
def decode (v : Coordinates) : Mat := ![![v 0,v 1,v 2,v 3],![v 4,v 5,v 6,v 7],![v 8,v 9,v 10,v 11],![v 12,v 13,v 14,v 15]]

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode]

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 16) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates]
  map_smul' r G := by fin_cases i <;> simp [coordinates]

def frame : Fin 16 → Mat := ![![![1,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]],![![0,1,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]],![![0,0,1,0],![0,0,0,0],![0,0,0,0],![0,0,0,0]],![![0,0,0,1],![0,0,0,0],![0,0,0,0],![0,0,0,0]],![![0,0,0,0],![1,0,0,0],![0,0,0,0],![0,0,0,0]],![![0,0,0,0],![0,1,0,0],![0,0,0,0],![0,0,0,0]],![![0,0,0,0],![0,0,1,0],![0,0,0,0],![0,0,0,0]],![![0,0,0,0],![0,0,0,1],![0,0,0,0],![0,0,0,0]],![![0,0,0,0],![0,0,0,0],![1,0,0,0],![0,0,0,0]],![![0,0,0,0],![0,0,0,0],![0,1,0,0],![0,0,0,0]],![![0,0,0,0],![0,0,0,0],![0,0,1,0],![0,0,0,0]],![![0,0,0,0],![0,0,0,0],![0,0,0,1],![0,0,0,0]],![![0,0,0,0],![0,0,0,0],![0,0,0,0],![1,0,0,0]],![![0,0,0,0],![0,0,0,0],![0,0,0,0],![0,1,0,0]],![![0,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,1,0]],![![0,0,0,0],![0,0,0,0],![0,0,0,0],![0,0,0,1]]]

theorem frame_expansion (G : Mat) :
    G = coordinates G 0 • frame 0 + coordinates G 1 • frame 1 + coordinates G 2 • frame 2 + coordinates G 3 • frame 3 + coordinates G 4 • frame 4 + coordinates G 5 • frame 5 + coordinates G 6 • frame 6 + coordinates G 7 • frame 7 + coordinates G 8 • frame 8 + coordinates G 9 • frame 9 + coordinates G 10 • frame 10 + coordinates G 11 • frame 11 + coordinates G 12 • frame 12 + coordinates G 13 • frame 13 + coordinates G 14 • frame 14 + coordinates G 15 • frame 15 := by
  ext i j
  fin_cases i <;> fin_cases j
  · change G 0 0 = G 0 0 * 1 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 0 1 = G 0 0 * 0 + G 0 1 * 1 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 0 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 1 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 0 3 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 1 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 1 0 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 1 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 1 1 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 1 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 1 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 1 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 1 3 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 1 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 2 0 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 1 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 2 1 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 1 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 2 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 1 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 2 3 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 1 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 3 0 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 1 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 3 1 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 1 + G 3 2 * 0 + G 3 3 * 0
    ring
  · change G 3 2 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 1 + G 3 3 * 0
    ring
  · change G 3 3 = G 0 0 * 0 + G 0 1 * 0 + G 0 2 * 0 + G 0 3 * 0 + G 1 0 * 0 + G 1 1 * 0 + G 1 2 * 0 + G 1 3 * 0 + G 2 0 * 0 + G 2 1 * 0 + G 2 2 * 0 + G 2 3 * 0 + G 3 0 * 0 + G 3 1 * 0 + G 3 2 * 0 + G 3 3 * 1
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
    c 9 * v 0 * v 9 +
    c 10 * v 0 * v 10 +
    c 11 * v 0 * v 11 +
    c 12 * v 0 * v 12 +
    c 13 * v 0 * v 13 +
    c 14 * v 0 * v 14 +
    c 15 * v 0 * v 15 +
    c 16 * v 1 * v 1 +
    c 17 * v 1 * v 2 +
    c 18 * v 1 * v 3 +
    c 19 * v 1 * v 4 +
    c 20 * v 1 * v 5 +
    c 21 * v 1 * v 6 +
    c 22 * v 1 * v 7 +
    c 23 * v 1 * v 8 +
    c 24 * v 1 * v 9 +
    c 25 * v 1 * v 10 +
    c 26 * v 1 * v 11 +
    c 27 * v 1 * v 12 +
    c 28 * v 1 * v 13 +
    c 29 * v 1 * v 14 +
    c 30 * v 1 * v 15 +
    c 31 * v 2 * v 2 +
    c 32 * v 2 * v 3 +
    c 33 * v 2 * v 4 +
    c 34 * v 2 * v 5 +
    c 35 * v 2 * v 6 +
    c 36 * v 2 * v 7 +
    c 37 * v 2 * v 8 +
    c 38 * v 2 * v 9 +
    c 39 * v 2 * v 10 +
    c 40 * v 2 * v 11 +
    c 41 * v 2 * v 12 +
    c 42 * v 2 * v 13 +
    c 43 * v 2 * v 14 +
    c 44 * v 2 * v 15 +
    c 45 * v 3 * v 3 +
    c 46 * v 3 * v 4 +
    c 47 * v 3 * v 5 +
    c 48 * v 3 * v 6 +
    c 49 * v 3 * v 7 +
    c 50 * v 3 * v 8 +
    c 51 * v 3 * v 9 +
    c 52 * v 3 * v 10 +
    c 53 * v 3 * v 11 +
    c 54 * v 3 * v 12 +
    c 55 * v 3 * v 13 +
    c 56 * v 3 * v 14 +
    c 57 * v 3 * v 15 +
    c 58 * v 4 * v 4 +
    c 59 * v 4 * v 5 +
    c 60 * v 4 * v 6 +
    c 61 * v 4 * v 7 +
    c 62 * v 4 * v 8 +
    c 63 * v 4 * v 9 +
    c 64 * v 4 * v 10 +
    c 65 * v 4 * v 11 +
    c 66 * v 4 * v 12 +
    c 67 * v 4 * v 13 +
    c 68 * v 4 * v 14 +
    c 69 * v 4 * v 15 +
    c 70 * v 5 * v 5 +
    c 71 * v 5 * v 6 +
    c 72 * v 5 * v 7 +
    c 73 * v 5 * v 8 +
    c 74 * v 5 * v 9 +
    c 75 * v 5 * v 10 +
    c 76 * v 5 * v 11 +
    c 77 * v 5 * v 12 +
    c 78 * v 5 * v 13 +
    c 79 * v 5 * v 14 +
    c 80 * v 5 * v 15 +
    c 81 * v 6 * v 6 +
    c 82 * v 6 * v 7 +
    c 83 * v 6 * v 8 +
    c 84 * v 6 * v 9 +
    c 85 * v 6 * v 10 +
    c 86 * v 6 * v 11 +
    c 87 * v 6 * v 12 +
    c 88 * v 6 * v 13 +
    c 89 * v 6 * v 14 +
    c 90 * v 6 * v 15 +
    c 91 * v 7 * v 7 +
    c 92 * v 7 * v 8 +
    c 93 * v 7 * v 9 +
    c 94 * v 7 * v 10 +
    c 95 * v 7 * v 11 +
    c 96 * v 7 * v 12 +
    c 97 * v 7 * v 13 +
    c 98 * v 7 * v 14 +
    c 99 * v 7 * v 15 +
    c 100 * v 8 * v 8 +
    c 101 * v 8 * v 9 +
    c 102 * v 8 * v 10 +
    c 103 * v 8 * v 11 +
    c 104 * v 8 * v 12 +
    c 105 * v 8 * v 13 +
    c 106 * v 8 * v 14 +
    c 107 * v 8 * v 15 +
    c 108 * v 9 * v 9 +
    c 109 * v 9 * v 10 +
    c 110 * v 9 * v 11 +
    c 111 * v 9 * v 12 +
    c 112 * v 9 * v 13 +
    c 113 * v 9 * v 14 +
    c 114 * v 9 * v 15 +
    c 115 * v 10 * v 10 +
    c 116 * v 10 * v 11 +
    c 117 * v 10 * v 12 +
    c 118 * v 10 * v 13 +
    c 119 * v 10 * v 14 +
    c 120 * v 10 * v 15 +
    c 121 * v 11 * v 11 +
    c 122 * v 11 * v 12 +
    c 123 * v 11 * v 13 +
    c 124 * v 11 * v 14 +
    c 125 * v 11 * v 15 +
    c 126 * v 12 * v 12 +
    c 127 * v 12 * v 13 +
    c 128 * v 12 * v 14 +
    c 129 * v 12 * v 15 +
    c 130 * v 13 * v 13 +
    c 131 * v 13 * v 14 +
    c 132 * v 13 * v 15 +
    c 133 * v 14 * v 14 +
    c 134 * v 14 * v 15 +
    c 135 * v 15 * v 15

def monomial (i j : Fin 16) : Quad :=
  QuadraticMap.linMulLin (coordinateMap i) (coordinateMap j)

theorem quadratic_representation (Q : Quad) :
    ∃ c : Coefficients, ∀ G, Q G = polynomial c (coordinates G) := by
  let b := Q.associated
  let c : Coefficients := ![b (frame 0) (frame 0),b (frame 0) (frame 1) + b (frame 1) (frame 0),b (frame 0) (frame 2) + b (frame 2) (frame 0),b (frame 0) (frame 3) + b (frame 3) (frame 0),b (frame 0) (frame 4) + b (frame 4) (frame 0),b (frame 0) (frame 5) + b (frame 5) (frame 0),b (frame 0) (frame 6) + b (frame 6) (frame 0),b (frame 0) (frame 7) + b (frame 7) (frame 0),b (frame 0) (frame 8) + b (frame 8) (frame 0),b (frame 0) (frame 9) + b (frame 9) (frame 0),b (frame 0) (frame 10) + b (frame 10) (frame 0),b (frame 0) (frame 11) + b (frame 11) (frame 0),b (frame 0) (frame 12) + b (frame 12) (frame 0),b (frame 0) (frame 13) + b (frame 13) (frame 0),b (frame 0) (frame 14) + b (frame 14) (frame 0),b (frame 0) (frame 15) + b (frame 15) (frame 0),b (frame 1) (frame 1),b (frame 1) (frame 2) + b (frame 2) (frame 1),b (frame 1) (frame 3) + b (frame 3) (frame 1),b (frame 1) (frame 4) + b (frame 4) (frame 1),b (frame 1) (frame 5) + b (frame 5) (frame 1),b (frame 1) (frame 6) + b (frame 6) (frame 1),b (frame 1) (frame 7) + b (frame 7) (frame 1),b (frame 1) (frame 8) + b (frame 8) (frame 1),b (frame 1) (frame 9) + b (frame 9) (frame 1),b (frame 1) (frame 10) + b (frame 10) (frame 1),b (frame 1) (frame 11) + b (frame 11) (frame 1),b (frame 1) (frame 12) + b (frame 12) (frame 1),b (frame 1) (frame 13) + b (frame 13) (frame 1),b (frame 1) (frame 14) + b (frame 14) (frame 1),b (frame 1) (frame 15) + b (frame 15) (frame 1),b (frame 2) (frame 2),b (frame 2) (frame 3) + b (frame 3) (frame 2),b (frame 2) (frame 4) + b (frame 4) (frame 2),b (frame 2) (frame 5) + b (frame 5) (frame 2),b (frame 2) (frame 6) + b (frame 6) (frame 2),b (frame 2) (frame 7) + b (frame 7) (frame 2),b (frame 2) (frame 8) + b (frame 8) (frame 2),b (frame 2) (frame 9) + b (frame 9) (frame 2),b (frame 2) (frame 10) + b (frame 10) (frame 2),b (frame 2) (frame 11) + b (frame 11) (frame 2),b (frame 2) (frame 12) + b (frame 12) (frame 2),b (frame 2) (frame 13) + b (frame 13) (frame 2),b (frame 2) (frame 14) + b (frame 14) (frame 2),b (frame 2) (frame 15) + b (frame 15) (frame 2),b (frame 3) (frame 3),b (frame 3) (frame 4) + b (frame 4) (frame 3),b (frame 3) (frame 5) + b (frame 5) (frame 3),b (frame 3) (frame 6) + b (frame 6) (frame 3),b (frame 3) (frame 7) + b (frame 7) (frame 3),b (frame 3) (frame 8) + b (frame 8) (frame 3),b (frame 3) (frame 9) + b (frame 9) (frame 3),b (frame 3) (frame 10) + b (frame 10) (frame 3),b (frame 3) (frame 11) + b (frame 11) (frame 3),b (frame 3) (frame 12) + b (frame 12) (frame 3),b (frame 3) (frame 13) + b (frame 13) (frame 3),b (frame 3) (frame 14) + b (frame 14) (frame 3),b (frame 3) (frame 15) + b (frame 15) (frame 3),b (frame 4) (frame 4),b (frame 4) (frame 5) + b (frame 5) (frame 4),b (frame 4) (frame 6) + b (frame 6) (frame 4),b (frame 4) (frame 7) + b (frame 7) (frame 4),b (frame 4) (frame 8) + b (frame 8) (frame 4),b (frame 4) (frame 9) + b (frame 9) (frame 4),b (frame 4) (frame 10) + b (frame 10) (frame 4),b (frame 4) (frame 11) + b (frame 11) (frame 4),b (frame 4) (frame 12) + b (frame 12) (frame 4),b (frame 4) (frame 13) + b (frame 13) (frame 4),b (frame 4) (frame 14) + b (frame 14) (frame 4),b (frame 4) (frame 15) + b (frame 15) (frame 4),b (frame 5) (frame 5),b (frame 5) (frame 6) + b (frame 6) (frame 5),b (frame 5) (frame 7) + b (frame 7) (frame 5),b (frame 5) (frame 8) + b (frame 8) (frame 5),b (frame 5) (frame 9) + b (frame 9) (frame 5),b (frame 5) (frame 10) + b (frame 10) (frame 5),b (frame 5) (frame 11) + b (frame 11) (frame 5),b (frame 5) (frame 12) + b (frame 12) (frame 5),b (frame 5) (frame 13) + b (frame 13) (frame 5),b (frame 5) (frame 14) + b (frame 14) (frame 5),b (frame 5) (frame 15) + b (frame 15) (frame 5),b (frame 6) (frame 6),b (frame 6) (frame 7) + b (frame 7) (frame 6),b (frame 6) (frame 8) + b (frame 8) (frame 6),b (frame 6) (frame 9) + b (frame 9) (frame 6),b (frame 6) (frame 10) + b (frame 10) (frame 6),b (frame 6) (frame 11) + b (frame 11) (frame 6),b (frame 6) (frame 12) + b (frame 12) (frame 6),b (frame 6) (frame 13) + b (frame 13) (frame 6),b (frame 6) (frame 14) + b (frame 14) (frame 6),b (frame 6) (frame 15) + b (frame 15) (frame 6),b (frame 7) (frame 7),b (frame 7) (frame 8) + b (frame 8) (frame 7),b (frame 7) (frame 9) + b (frame 9) (frame 7),b (frame 7) (frame 10) + b (frame 10) (frame 7),b (frame 7) (frame 11) + b (frame 11) (frame 7),b (frame 7) (frame 12) + b (frame 12) (frame 7),b (frame 7) (frame 13) + b (frame 13) (frame 7),b (frame 7) (frame 14) + b (frame 14) (frame 7),b (frame 7) (frame 15) + b (frame 15) (frame 7),b (frame 8) (frame 8),b (frame 8) (frame 9) + b (frame 9) (frame 8),b (frame 8) (frame 10) + b (frame 10) (frame 8),b (frame 8) (frame 11) + b (frame 11) (frame 8),b (frame 8) (frame 12) + b (frame 12) (frame 8),b (frame 8) (frame 13) + b (frame 13) (frame 8),b (frame 8) (frame 14) + b (frame 14) (frame 8),b (frame 8) (frame 15) + b (frame 15) (frame 8),b (frame 9) (frame 9),b (frame 9) (frame 10) + b (frame 10) (frame 9),b (frame 9) (frame 11) + b (frame 11) (frame 9),b (frame 9) (frame 12) + b (frame 12) (frame 9),b (frame 9) (frame 13) + b (frame 13) (frame 9),b (frame 9) (frame 14) + b (frame 14) (frame 9),b (frame 9) (frame 15) + b (frame 15) (frame 9),b (frame 10) (frame 10),b (frame 10) (frame 11) + b (frame 11) (frame 10),b (frame 10) (frame 12) + b (frame 12) (frame 10),b (frame 10) (frame 13) + b (frame 13) (frame 10),b (frame 10) (frame 14) + b (frame 14) (frame 10),b (frame 10) (frame 15) + b (frame 15) (frame 10),b (frame 11) (frame 11),b (frame 11) (frame 12) + b (frame 12) (frame 11),b (frame 11) (frame 13) + b (frame 13) (frame 11),b (frame 11) (frame 14) + b (frame 14) (frame 11),b (frame 11) (frame 15) + b (frame 15) (frame 11),b (frame 12) (frame 12),b (frame 12) (frame 13) + b (frame 13) (frame 12),b (frame 12) (frame 14) + b (frame 14) (frame 12),b (frame 12) (frame 15) + b (frame 15) (frame 12),b (frame 13) (frame 13),b (frame 13) (frame 14) + b (frame 14) (frame 13),b (frame 13) (frame 15) + b (frame 15) (frame 13),b (frame 14) (frame 14),b (frame 14) (frame 15) + b (frame 15) (frame 14),b (frame 15) (frame 15)]
  refine ⟨c, fun G => ?_⟩
  have hb : Q G = b G G := (Q.associated_eq_self_apply ℝ G).symm
  rw [hb, congrArg (fun Z => b Z Z) (frame_expansion G)]
  simp only [polynomial, c, Matrix.cons_val]
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
  ring

end
end S11D4Invariants
