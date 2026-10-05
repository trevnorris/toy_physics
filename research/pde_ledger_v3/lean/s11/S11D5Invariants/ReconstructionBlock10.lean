import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock7
import S11D5Invariants.ConstraintBlock8
import S11D5Invariants.ConstraintBlock11
import S11D5Invariants.ConstraintBlock12
import S11D5Invariants.ConstraintBlock13
import S11D5Invariants.ConstraintBlock14
import S11D5Invariants.ConstraintBlock15

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock10 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 203 = 0) ∧
    (c 204 = 0) ∧
    (c 205 = 1 * c 25) ∧
    (c 206 = 0) ∧
    (c 207 = 0) ∧
    (c 208 = 0) ∧
    (c 209 = 0) ∧
    (c 210 = 0) ∧
    (c 211 = 0) ∧
    (c 212 = 0) ∧
    (c 213 = 0) ∧
    (c 214 = 0) ∧
    (c 215 = 0) ∧
    (c 216 = 0) ∧
    (c 217 = 0) ∧
    (c 218 = 0) ∧
    (c 219 = 0) ∧
    (c 220 = 1 * c 25) ∧
    (c 221 = 0) ∧
    (c 222 = 0) := by
  obtain ⟨e0,_,e2,_,e4,_,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,_,e30,_,_,_,_,_,_,e37,e38,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,e43,_,_,e46,_,_,_,e50,_,_,_,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,_,_,_,_,e65,e66,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,e107,e108,_,_,_,_,_,e114,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,_,e127,_,_,e130,e131,_,_,_,e135,e136,_,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e155,e156,e157,e158,e159⟩ := necessaryBlock7 hQ c hc
  obtain ⟨e160,_,_,e163,e164,e165,_,_,e168,e169,e170,e171,e172,e173,e174,e175,e176,e177,e178,e179⟩ := necessaryBlock8 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,_,_,_,_,_,_,e235,_,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨_,_,_,_,_,e245,_,_,_,_,_,_,_,_,_,_,_,_,e258,e259⟩ := necessaryBlock12 hQ c hc
  obtain ⟨e260,e261,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock13 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,e299⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,_,_,_,_,_,_,_,e315,_,_,e318,_⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination (-1) * e90 + ((1/2)) * e107 + ((-1/2)) * e155
  · linear_combination (-1) * e90 + ((1/2)) * e108 + ((-1/2)) * e156
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e10 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + (-1) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224
  · linear_combination ((-1/2)) * e157
  · linear_combination ((-1/2)) * e158 + ((-1/2)) * e171
  · linear_combination ((-1/2)) * e159 + ((-1/2)) * e172
  · linear_combination ((-1/2)) * e160 + ((-1/2)) * e173
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((1/2)) * e37 + (3) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e60 + ((-1/2)) * e114 + ((-1/2)) * e118 + ((-1/2)) * e121 + (1) * e127 + ((1/2)) * e131 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e258
  · linear_combination (4) * e0 + ((3/2)) * e2 + ((-1/2)) * e4 + (-1) * e6 + ((-1/2)) * e9 + (-1) * e10 + ((1/2)) * e14 + (-2) * e24 + ((1/2)) * e30 + (1) * e38 + (-3) * e46 + (1) * e50 + (2) * e55 + ((-1/2)) * e58 + ((-1/2)) * e60 + (-1) * e90 + ((1/2)) * e91 + ((3/2)) * e118 + ((-1/2)) * e120 + (1) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-2) * e127 + ((-1/2)) * e130 + ((1/2)) * e131 + (1) * e220 + (-1) * e221 + (2) * e224 + (-1) * e235 + (-1) * e259 + (-1) * e293 + (1) * e294 + (-1) * e299 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318
  · linear_combination ((-1/2)) * e163 + ((-1/2)) * e174
  · linear_combination ((-1/2)) * e164 + ((-1/2)) * e175
  · linear_combination ((-1/2)) * e165 + ((-1/2)) * e176
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((1/2)) * e42 + (3) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e65 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((-1/2)) * e121 + (1) * e127 + ((1/2)) * e135 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e260
  · linear_combination (2) * e0 + (1) * e2 + ((-1/2)) * e6 + ((-1/2)) * e10 + (-2) * e24 + ((1/2)) * e42 + (1) * e43 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e66 + ((1/2)) * e116 + (1) * e118 + ((1/2)) * e121 + (-1) * e127 + ((-1/2)) * e136 + (1) * e220 + (-1) * e221 + (1) * e224 + (-1) * e245 + (-1) * e261 + ((-1/2)) * e302 + ((-1/2)) * e303
  · linear_combination ((-1/2)) * e168 + ((-1/2)) * e177
  · linear_combination ((-1/2)) * e169 + ((-1/2)) * e178
  · linear_combination ((-1/2)) * e170 + ((-1/2)) * e179
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e10 + (-1) * e46 + ((1/2)) * e50 + ((1/2)) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224
  · linear_combination (-1) * e46 + ((1/2)) * e55 + ((-1/2)) * e127 + ((1/2)) * e158 + ((-1/2)) * e171
  · linear_combination (-1) * e46 + ((1/2)) * e55 + ((-1/2)) * e127 + ((1/2)) * e159 + ((-1/2)) * e172

end
end S11D5Invariants
