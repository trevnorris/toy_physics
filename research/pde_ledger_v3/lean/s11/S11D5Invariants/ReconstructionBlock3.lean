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

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock3 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 63 = 0) ∧
    (c 64 = 0) ∧
    (c 65 = 0) ∧
    (c 66 = 0) ∧
    (c 67 = 0) ∧
    (c 68 = 0) ∧
    (c 69 = 0) ∧
    (c 70 = 0) ∧
    (c 71 = 0) ∧
    (c 72 = 1 * c 25) ∧
    (c 73 = 0) ∧
    (c 74 = 0) ∧
    (c 75 = 0) ∧
    (c 76 = 0) ∧
    (c 77 = 0) ∧
    (c 78 = 0) ∧
    (c 79 = 0) ∧
    (c 80 = 0) ∧
    (c 81 = 0) ∧
    (c 82 = 0) := by
  obtain ⟨e0,_,e2,e3,e4,_,e6,e7,_,e9,e10,_,_,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,e27,_,e29,e30,_,_,_,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,e41,e42,_,_,_,e46,_,_,_,e50,_,e52,_,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,_,e65,e66,e67,e68,e69,_,e71,e72,e73,_,_,_,e77,e78,e79⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,e111,_,_,e114,_,e116,e117,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,_,e127,_,_,e130,e131,e132,e133,e134,e135,e136,e137,e138,e139⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e140,e141,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,_,e230,_,_,_,_,e235,_,e237,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,e241,_,_,e244,e245,e246,_,_,e249,e250,e251,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,e299⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,e304,_,_,_,_,_,_,_,_,_,_,e315,_,_,e318,_⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination (-2) * e0 + ((1/2)) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + ((-1/2)) * e30 + (1) * e37 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((1/2)) * e58 + ((1/2)) * e60 + (-1) * e61 + (1) * e90 + ((-1/2)) * e91 + (1) * e114 + ((1/2)) * e118 + ((1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((1/2)) * e130 + ((-1/2)) * e131 + (-1) * e132 + (1) * e220 + (-1) * e221 + (-1) * e224 + (1) * e235 + (-2) * e240 + (-1) * e241 + (1) * e293 + (-1) * e294 + (1) * e299 + ((-1/2)) * e300 + ((-1/2)) * e301 + (1) * e315 + (-1) * e318
  · linear_combination ((-1/2)) * e61 + ((-1/2)) * e132
  · linear_combination ((-1/2)) * e62 + ((-1/2)) * e133
  · linear_combination ((-1/2)) * e63 + ((-1/2)) * e134
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e41 + (1) * e46 + ((-1/2)) * e65 + ((1/2)) * e117 + ((-1/2)) * e118 + ((1/2)) * e135 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e244
  · linear_combination (1) * e2 + ((1/2)) * e42 + ((-1/2)) * e66 + ((1/2)) * e116 + (1) * e118 + ((-1/2)) * e136 + (1) * e220 + (-1) * e221 + (-1) * e245 + (-1) * e246 + ((1/2)) * e302 + ((1/2)) * e303
  · linear_combination ((-1/2)) * e66 + ((-1/2)) * e136
  · linear_combination ((-1/2)) * e67 + ((-1/2)) * e137
  · linear_combination ((-1/2)) * e68 + ((-1/2)) * e138
  · linear_combination ((1/2)) * e3 + ((1/2)) * e119 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e293
  · linear_combination ((-1/2)) * e3 + (1) * e24 + ((-1/2)) * e27 + (1) * e90 + ((1/2)) * e111 + ((-1/2)) * e119 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e230 + (-1) * e293 + (1) * e304
  · linear_combination (-1) * e24 + ((1/2)) * e29 + (1) * e69 + ((-1/2)) * e71
  · linear_combination ((-1/2)) * e7 + ((-1/2)) * e72
  · linear_combination (-1) * e2 + ((-1/2)) * e3 + ((1/2)) * e4 + ((-1/2)) * e30 + (1) * e46 + (-1) * e52 + ((1/2)) * e58 + (-2) * e69 + (1) * e90 + ((-1/2)) * e91 + (-1) * e118 + ((-1/2)) * e119 + ((1/2)) * e120 + ((1/2)) * e130 + (-1) * e220 + (1) * e221 + (1) * e235 + (1) * e237 + ((-1/2)) * e300 + ((-1/2)) * e301 + (1) * e315 + (-1) * e318
  · linear_combination ((-1/2)) * e73
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e30 + ((1/2)) * e58 + ((-1/2)) * e79 + ((-1/2)) * e91 + ((-1/2)) * e118 + ((1/2)) * e130 + ((-1/2)) * e141 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e235 + (-1) * e251
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((1/2)) * e29 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + (1) * e69 + ((-1/2)) * e71 + ((-1/2)) * e118 + ((-1/2)) * e121 + (1) * e127 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e249
  · linear_combination (2) * e0 + (1) * e2 + ((1/2)) * e3 + ((-1/2)) * e4 + ((-1/2)) * e6 + ((-1/2)) * e10 + ((1/2)) * e30 + (-3) * e46 + ((1/2)) * e50 + (1) * e52 + ((1/2)) * e55 + ((-1/2)) * e58 + (2) * e69 + (-1) * e90 + ((1/2)) * e91 + (1) * e118 + ((1/2)) * e119 + ((-1/2)) * e120 + ((1/2)) * e121 + ((-1/2)) * e127 + ((-1/2)) * e130 + (1) * e220 + (-1) * e221 + (1) * e224 + (-1) * e235 + (-1) * e237 + (-1) * e250 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318
  · linear_combination ((-1/2)) * e77 + ((-1/2)) * e139
  · linear_combination ((-1/2)) * e78 + ((-1/2)) * e140

end
end S11D5Invariants
