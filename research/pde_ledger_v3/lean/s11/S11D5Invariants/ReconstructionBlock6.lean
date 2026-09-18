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

theorem reconstructedBlock6 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 123 = 0) ∧
    (c 124 = 0) ∧
    (c 125 = 0) ∧
    (c 126 = 0) ∧
    (c 127 = 0) ∧
    (c 128 = 0) ∧
    (c 129 = 0) ∧
    (c 130 = 0) ∧
    (c 131 = 0) ∧
    (c 132 = 0) ∧
    (c 133 = 0) ∧
    (c 134 = 0) ∧
    (c 135 = 1 * (c 6/2) + 1 * (c 29/2) + 1 * c 25) ∧
    (c 136 = 0) ∧
    (c 137 = 0) ∧
    (c 138 = 0) ∧
    (c 139 = 0) ∧
    (c 140 = 0) ∧
    (c 141 = 2 * (c 6/2)) ∧
    (c 142 = 0) := by
  obtain ⟨e0,e1,e2,e3,e4,e5,e6,e7,_,e9,e10,e11,e12,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,e29,e30,_,_,_,e34,e35,e36,e37,e38,e39⟩ := necessaryBlock1 hQ c hc
  obtain ⟨e40,e41,e42,e43,e44,e45,e46,_,_,_,e50,_,_,_,_,e55,e56,e57,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,_,e65,_,e67,e68,e69,_,e71,e72,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,e114,e115,e116,e117,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,_,e127,e128,e129,e130,e131,e132,e133,e134,e135,_,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e220,e221,e222,e223,e224,e225,_,_,_,_,_,_,e232,_,e234,e235,_,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,e242,e243,_,_,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,e299⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,_,_,_,_,_,_,_,e315,_,_,e318,_⟩ := necessaryBlock15 hQ c hc
  obtain ⟨_,e321⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((1/2)) * e3 + ((-1/2)) * e4 + (1) * e24 + ((-1/2)) * e29 + ((1/2)) * e30 + (-1) * e34 + ((1/2)) * e57 + ((-1/2)) * e58 + (1) * e69 + ((-1/2)) * e71 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e119 + ((-1/2)) * e120 + ((1/2)) * e129 + ((-1/2)) * e130 + (1) * e234 + (-1) * e235 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318
  · linear_combination (1) * e24 + (-1) * e35 + ((1/2)) * e300 + ((1/2)) * e301
  · linear_combination (-1) * e24 + ((1/2)) * e37 + (1) * e46 + ((-1/2)) * e60 + ((-1/2)) * e114 + ((1/2)) * e131
  · linear_combination ((-1/2)) * e36 + ((-1/2)) * e115
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + (1) * e24 + ((-1/2)) * e30 + ((1/2)) * e37 + (-1) * e38 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((1/2)) * e58 + ((1/2)) * e60 + ((-1/2)) * e61 + (1) * e90 + ((-1/2)) * e91 + ((1/2)) * e114 + ((-1/2)) * e118 + ((1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((1/2)) * e130 + ((-1/2)) * e131 + ((-1/2)) * e132 + (-1) * e224 + (1) * e235 + (-1) * e240 + (1) * e293 + (-1) * e294 + (1) * e299 + ((-1/2)) * e300 + ((-1/2)) * e301 + (1) * e315 + (-1) * e318
  · linear_combination ((-1/2)) * e2 + (1) * e24 + (-1) * e39 + ((1/2)) * e62 + ((-1/2)) * e118 + ((1/2)) * e133 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e242
  · linear_combination ((-1/2)) * e2 + (1) * e24 + (-1) * e40 + ((1/2)) * e63 + ((-1/2)) * e118 + ((1/2)) * e134 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e243
  · linear_combination (-1) * e24 + ((1/2)) * e42 + (1) * e46 + ((-1/2)) * e65 + ((-1/2)) * e116 + ((1/2)) * e135
  · linear_combination ((-1/2)) * e41 + ((-1/2)) * e117
  · linear_combination (1) * e24 + (-1) * e43 + ((1/2)) * e302 + ((1/2)) * e303
  · linear_combination ((-1/2)) * e2 + (1) * e24 + (-1) * e44 + ((1/2)) * e67 + ((-1/2)) * e118 + ((1/2)) * e137 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e247
  · linear_combination ((-1/2)) * e2 + (1) * e24 + (-1) * e45 + ((1/2)) * e68 + ((-1/2)) * e118 + ((1/2)) * e138 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e248
  · linear_combination ((95/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e5 + ((7/48)) * e6 + ((7/48)) * e9 + ((7/48)) * e10 + ((19/24)) * e24 + ((7/12)) * e46 + ((-7/48)) * e50 + ((-7/24)) * e55 + ((-7/48)) * e118 + ((-7/48)) * e121 + ((7/48)) * e122 + ((7/24)) * e127 + ((7/24)) * e221 + ((7/24)) * e222 + ((-7/24)) * e224 + ((-625/288)) * e321
  · linear_combination (-1) * e0 + ((1/2)) * e2 + (-1) * e46 + ((-1/2)) * e118
  · linear_combination (-1) * e0 + ((1/2)) * e3 + (-1) * e69 + ((-1/2)) * e119
  · linear_combination (-1) * e0 + ((1/2)) * e4 + (-1) * e90 + ((-1/2)) * e120
  · linear_combination ((-1/2)) * e10 + ((-1/2)) * e121
  · linear_combination (-1) * e0 + ((1/2)) * e9 + (-1) * e46 + ((1/2)) * e55 + ((-1/2)) * e122 + ((-1/2)) * e127
  · linear_combination (1) * e2 + (1) * e5 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + (1) * e11 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e56 + (1) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + (-1) * e127 + ((-1/2)) * e128 + (1) * e220 + (-1) * e221 + (-1) * e222 + (1) * e223 + (1) * e224 + (-1) * e232
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e7 + (1) * e12 + ((1/2)) * e29 + ((-1/2)) * e57 + (-1) * e69 + ((1/2)) * e71 + ((1/2)) * e72 + ((1/2)) * e118 + ((-1/2)) * e129 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e225 + (-1) * e234

end
end S11D5Invariants
