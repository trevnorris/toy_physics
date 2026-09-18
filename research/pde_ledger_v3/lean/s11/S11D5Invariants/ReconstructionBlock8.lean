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

theorem reconstructedBlock8 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 163 = 0) ∧
    (c 164 = 0) ∧
    (c 165 = 0) ∧
    (c 166 = 0) ∧
    (c 167 = 0) ∧
    (c 168 = 0) ∧
    (c 169 = 0) ∧
    (c 170 = 0) ∧
    (c 171 = 0) ∧
    (c 172 = 1 * c 25) ∧
    (c 173 = 0) ∧
    (c 174 = 0) ∧
    (c 175 = 0) ∧
    (c 176 = 0) ∧
    (c 177 = 0) ∧
    (c 178 = 0) ∧
    (c 179 = 0) ∧
    (c 180 = 2 * (c 29/2)) ∧
    (c 181 = 0) ∧
    (c 182 = 0) := by
  obtain ⟨e0,_,e2,e3,e4,_,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,e27,_,e29,e30,_,_,_,_,_,e36,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,e41,e42,_,_,_,e46,_,_,_,e50,_,e52,_,_,e55,_,_,e58,e59⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,e64,e65,e66,e67,e68,e69,e70,e71,_,_,_,e75,e76,e77,e78,e79⟩ := necessaryBlock3 hQ c hc
  obtain ⟨e80,e81,e82,e83,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,e111,_,_,e114,e115,e116,e117,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,_,e127,_,_,e130,e131,e132,e133,e134,e135,e136,e137,e138,e139⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e140,e141,e142,e143,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,_,e230,e231,_,_,_,e235,_,e237,_,e239⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,_,_,e244,e245,e246,_,_,e249,e250,_,e252,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,e304,e305,_,_,_,_,_,_,_,_,_,e315,_,_,e318,_⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e36 + (-1) * e46 + (1) * e59 + ((1/2)) * e115 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e239
  · linear_combination (-1) * e46 + ((1/2)) * e61 + ((-1/2)) * e132
  · linear_combination (-1) * e46 + ((1/2)) * e62 + ((-1/2)) * e133
  · linear_combination (-1) * e46 + ((1/2)) * e63 + ((-1/2)) * e134
  · linear_combination (-1) * e2 + ((-1/2)) * e42 + ((-1/2)) * e65 + ((1/2)) * e66 + ((-1/2)) * e116 + (-1) * e118 + ((-1/2)) * e135 + ((1/2)) * e136 + (-1) * e220 + (1) * e221 + (1) * e245 + (1) * e246 + ((-1/2)) * e302 + ((-1/2)) * e303
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e41 + (-1) * e46 + (1) * e64 + ((1/2)) * e117 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e244
  · linear_combination (-1) * e46 + ((1/2)) * e66 + ((-1/2)) * e136
  · linear_combination (-1) * e46 + ((1/2)) * e67 + ((-1/2)) * e137
  · linear_combination (-1) * e46 + ((1/2)) * e68 + ((-1/2)) * e138
  · linear_combination ((1/2)) * e3 + (1) * e69 + ((1/2)) * e119 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e293
  · linear_combination ((-1/2)) * e3 + (1) * e24 + ((-1/2)) * e27 + (-1) * e69 + (1) * e70 + ((1/2)) * e111 + ((-1/2)) * e119 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e230 + (-1) * e293 + (1) * e304
  · linear_combination (-2) * e0 + (-1) * e2 + ((-1/2)) * e3 + ((1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e10 + ((-1/2)) * e30 + (2) * e46 + ((-1/2)) * e50 + (-1) * e52 + ((1/2)) * e58 + (-1) * e69 + (-1) * e76 + (1) * e90 + ((-1/2)) * e91 + (-1) * e118 + ((-1/2)) * e119 + ((1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e130 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e235 + (1) * e237 + (1) * e250 + ((-1/2)) * e300 + ((-1/2)) * e301 + (1) * e315 + (-1) * e318
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((1/2)) * e29 + (1) * e46 + ((-1/2)) * e50 + ((-1/2)) * e55 + ((-1/2)) * e71 + (1) * e75 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e127 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e249
  · linear_combination (-1) * e69 + ((1/2)) * e77 + ((-1/2)) * e139
  · linear_combination (-1) * e69 + ((1/2)) * e78 + ((-1/2)) * e140
  · linear_combination (-1) * e69 + ((1/2)) * e79 + ((-1/2)) * e141
  · linear_combination ((1/2)) * e2 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + (1) * e69 + (-1) * e81 + ((1/2)) * e82 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e142 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e240 + (1) * e252
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((-1/2)) * e3 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + (1) * e24 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((1/2)) * e60 + (-1) * e69 + (1) * e80 + ((-1/2)) * e118 + ((-1/2)) * e119 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((-1/2)) * e131 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e231 + (-1) * e293 + (-1) * e294 + (1) * e305
  · linear_combination (-1) * e69 + ((1/2)) * e82 + ((-1/2)) * e142
  · linear_combination (-1) * e69 + ((1/2)) * e83 + ((-1/2)) * e143

end
end S11D5Invariants
