import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock11
import S11D5Invariants.ConstraintBlock12
import S11D5Invariants.ConstraintBlock14
import S11D5Invariants.ConstraintBlock15
import S11D5Invariants.ConstraintBlock16

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock0 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 0 = 1 * (c 6/2) + 1 * (c 29/2) + 1 * c 25) ∧
    (c 1 = 0) ∧
    (c 2 = 0) ∧
    (c 3 = 0) ∧
    (c 4 = 0) ∧
    (c 5 = 0) ∧
    (c 7 = 0) ∧
    (c 8 = 0) ∧
    (c 9 = 0) ∧
    (c 10 = 0) ∧
    (c 11 = 0) ∧
    (c 12 = 2 * (c 6/2)) ∧
    (c 13 = 0) ∧
    (c 14 = 0) ∧
    (c 15 = 0) ∧
    (c 16 = 0) ∧
    (c 17 = 0) ∧
    (c 18 = 2 * (c 6/2)) ∧
    (c 19 = 0) ∧
    (c 20 = 0) := by
  obtain ⟨e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,_,_,_,e14,e15,_,_,_,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,e29,e30,_,_,_,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,_,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,e56,e57,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,_,_,_,_,_,e69,_,e71,e72,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,e92,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,e114,_,_,_,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,e123,e124,_,e126,e127,e128,e129,e130,e131,e132,e133,e134,_,_,_,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e220,e221,e222,e223,e224,e225,e226,e227,_,_,_,_,e232,_,e234,e235,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,e242,e243,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e295,e296,_,e298,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock15 hQ c hc
  obtain ⟨_,e321⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((59/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e5 + ((7/48)) * e6 + ((7/48)) * e9 + ((7/48)) * e10 + ((19/24)) * e24 + ((7/12)) * e46 + ((-7/48)) * e50 + ((-7/24)) * e55 + ((-7/48)) * e118 + ((-7/48)) * e121 + ((7/48)) * e122 + ((7/24)) * e127 + ((7/24)) * e221 + ((7/24)) * e222 + ((-7/24)) * e224 + ((-625/288)) * e321
  · linear_combination ((-1/2)) * e220 + ((-1/2)) * e221
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e118
  · linear_combination ((-1/2)) * e3 + ((-1/2)) * e119
  · linear_combination ((-1/2)) * e4 + ((-1/2)) * e120
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + (-1) * e24 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + (-1) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e222 + (1) * e224
  · linear_combination (1) * e0 + ((-1/2)) * e6 + (-1) * e46 + ((1/2)) * e50
  · linear_combination (1) * e0 + ((-1/2)) * e7 + (-1) * e69 + ((1/2)) * e72
  · linear_combination (1) * e0 + ((-1/2)) * e8 + (-1) * e90 + ((1/2)) * e92
  · linear_combination ((-1/2)) * e9 + ((-1/2)) * e122
  · linear_combination (1) * e0 + ((-1/2)) * e10 + (-1) * e46 + ((1/2)) * e55 + ((1/2)) * e121 + ((-1/2)) * e127
  · linear_combination (1) * e0 + (1) * e2 + (1) * e5 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e56 + (1) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + (-1) * e127 + ((-1/2)) * e128 + (1) * e220 + (-1) * e221 + (-1) * e222 + (1) * e223 + (1) * e224 + (-1) * e232
  · linear_combination (1) * e0 + ((1/2)) * e2 + ((-1/2)) * e7 + ((1/2)) * e29 + ((-1/2)) * e57 + (-1) * e69 + ((1/2)) * e71 + ((1/2)) * e72 + ((1/2)) * e118 + ((-1/2)) * e129 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e225 + (-1) * e234
  · linear_combination (1) * e0 + ((1/2)) * e2 + ((-1/2)) * e8 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e92 + ((1/2)) * e118 + ((-1/2)) * e130 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e226 + (-1) * e235
  · linear_combination ((-1/2)) * e14 + ((-1/2)) * e124
  · linear_combination (1) * e0 + ((-1/2)) * e15 + (-1) * e46 + ((1/2)) * e60 + ((1/2)) * e123 + ((-1/2)) * e131
  · linear_combination (1) * e0 + ((1/2)) * e2 + ((-1/2)) * e15 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((1/2)) * e123 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e227 + (-1) * e240
  · linear_combination ((1/2)) * e6 + ((1/2)) * e50 + ((-1/2)) * e62 + ((-1/2)) * e133 + (1) * e223 + (1) * e236 + (-1) * e242 + (1) * e295 + (-1) * e298
  · linear_combination (1) * e0 + (1) * e2 + ((-1/2)) * e8 + ((1/2)) * e30 + ((-1/2)) * e58 + ((-1/2)) * e63 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e92 + (1) * e118 + ((-1/2)) * e130 + ((-1/2)) * e134 + (1) * e220 + (-1) * e221 + (1) * e226 + (-1) * e235 + (-1) * e243 + (1) * e296 + ((-1/2)) * e300 + ((1/2)) * e301
  · linear_combination ((-1/2)) * e19 + ((-1/2)) * e126

end
end S11D5Invariants
