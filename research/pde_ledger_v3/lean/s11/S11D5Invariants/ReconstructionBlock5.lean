import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock7
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

theorem reconstructedBlock5 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 103 = 0) ∧
    (c 104 = 0) ∧
    (c 105 = 0) ∧
    (c 106 = 0) ∧
    (c 107 = 0) ∧
    (c 108 = 0) ∧
    (c 109 = 0) ∧
    (c 110 = 2 * (c 29/2)) ∧
    (c 111 = 0) ∧
    (c 112 = 0) ∧
    (c 113 = 0) ∧
    (c 114 = 0) ∧
    (c 115 = 1 * c 25) ∧
    (c 116 = 0) ∧
    (c 117 = 0) ∧
    (c 118 = 0) ∧
    (c 119 = 0) ∧
    (c 120 = 0) ∧
    (c 121 = 0) ∧
    (c 122 = 0) := by
  obtain ⟨e0,e1,e2,_,e4,e5,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,e25,e26,e27,_,_,e30,e31,e32,e33,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,e56,_,_,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨_,e61,_,_,_,_,e66,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,e97,e98,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,e101,e102,e103,_,_,e106,e107,e108,e109,e110,e111,e112,e113,e114,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,e126,e127,e128,_,_,_,e132,_,_,_,e136,_,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,e149,e150,e151,e152,e153,e154,e155,e156,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨e220,e221,e222,_,e224,_,_,_,_,_,_,e231,e232,_,_,_,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,_,_,_,e245,_,_,_,_,_,_,_,_,e254,_,e256,e257,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨_,_,_,_,_,e305,_,e307,_,_,_,_,_,_,_,e315,e316,_,_,_⟩ := necessaryBlock15 hQ c hc
  obtain ⟨e320,_⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((-1/2)) * e97 + ((-1/2)) * e149
  · linear_combination ((-1/2)) * e98 + ((-1/2)) * e150
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + ((1/2)) * e30 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + (1) * e90 + ((-1/2)) * e91 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e254 + (-1) * e294 + (1) * e307
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e37 + ((1/2)) * e61 + ((-1/2)) * e101 + ((-1/2)) * e114 + ((-1/2)) * e118 + ((1/2)) * e132 + ((-1/2)) * e151 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e240 + (-1) * e256
  · linear_combination ((-1/2)) * e101 + ((-1/2)) * e151
  · linear_combination ((-1/2)) * e102 + ((-1/2)) * e152
  · linear_combination ((-1/2)) * e103 + ((-1/2)) * e153
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((-1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e19 + (1) * e24 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e118 + ((-1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e126 + (1) * e127 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e231 + (-1) * e293 + (-1) * e294 + (1) * e305 + (-1) * e315 + (-1) * e316 + (1) * e320
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e42 + ((1/2)) * e66 + ((-1/2)) * e106 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((1/2)) * e136 + ((-1/2)) * e154 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e245 + (-1) * e257
  · linear_combination ((-1/2)) * e106 + ((-1/2)) * e154
  · linear_combination ((-1/2)) * e107 + ((-1/2)) * e155
  · linear_combination ((-1/2)) * e108 + ((-1/2)) * e156
  · linear_combination (1) * e24
  · linear_combination (1) * e0 + (-1) * e1 + (1) * e24 + ((1/2)) * e220 + ((1/2)) * e221
  · linear_combination ((-1/2)) * e25 + ((-1/2)) * e109
  · linear_combination ((-1/2)) * e26 + ((-1/2)) * e110
  · linear_combination ((-1/2)) * e27 + ((-1/2)) * e111
  · linear_combination (-1) * e24 + ((1/2)) * e32 + (1) * e46 + ((-1/2)) * e55 + ((-1/2)) * e112 + ((1/2)) * e127
  · linear_combination ((-1/2)) * e31 + ((-1/2)) * e113
  · linear_combination (-1) * e0 + ((-3/2)) * e2 + (-1) * e5 + (1) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + (1) * e24 + (-1) * e33 + (2) * e46 + (-1) * e55 + ((1/2)) * e56 + ((-3/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + (1) * e127 + ((1/2)) * e128 + ((-3/2)) * e220 + ((3/2)) * e221 + (1) * e222 + (-1) * e224 + (1) * e232 + (1) * e236

end
end S11D5Invariants
