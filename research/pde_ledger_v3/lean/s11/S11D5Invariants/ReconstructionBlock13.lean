import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock7
import S11D5Invariants.ConstraintBlock9
import S11D5Invariants.ConstraintBlock10
import S11D5Invariants.ConstraintBlock11
import S11D5Invariants.ConstraintBlock12
import S11D5Invariants.ConstraintBlock13
import S11D5Invariants.ConstraintBlock14
import S11D5Invariants.ConstraintBlock15
import S11D5Invariants.ConstraintBlock16

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock13 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 263 = 0) ∧
    (c 264 = 0) ∧
    (c 265 = 0) ∧
    (c 266 = 0) ∧
    (c 267 = 2 * (c 29/2)) ∧
    (c 268 = 0) ∧
    (c 269 = 0) ∧
    (c 270 = 1 * c 25) ∧
    (c 271 = 0) ∧
    (c 272 = 0) ∧
    (c 273 = 0) ∧
    (c 274 = 0) ∧
    (c 275 = 0) ∧
    (c 276 = 0) ∧
    (c 277 = 0) ∧
    (c 278 = 0) ∧
    (c 279 = 0) ∧
    (c 280 = 1 * c 25) ∧
    (c 281 = 0) ∧
    (c 282 = 0) := by
  obtain ⟨e0,_,e2,_,e4,_,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,_,e30,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,_,_,_,_,e65,e66,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,e102,e103,e104,_,_,e107,e108,_,_,_,_,_,_,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,e126,e127,_,_,e130,e131,_,_,_,e135,e136,_,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,e152,e153,_,e155,e156,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e190,e191,e192,e193,e194,e195,_,_,e198,e199⟩ := necessaryBlock9 hQ c hc
  obtain ⟨e200,e201,e202,e203,e204,e205,e206,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock10 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,_,_,e231,_,_,_,e235,_,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨_,_,_,_,_,e245,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨e260,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e278,e279⟩ := necessaryBlock13 hQ c hc
  obtain ⟨e280,e281,e282,e283,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨_,_,_,_,_,e305,_,_,e308,_,_,_,_,_,_,e315,e316,_,_,_⟩ := necessaryBlock15 hQ c hc
  obtain ⟨e320,_⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((1/2)) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e102 + ((1/2)) * e118 + ((-1/2)) * e130 + ((-1/2)) * e152 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e235 + (1) * e278
  · linear_combination ((1/2)) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e103 + ((1/2)) * e118 + ((-1/2)) * e130 + ((-1/2)) * e153 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e235 + (1) * e279
  · linear_combination ((-1/2)) * e190 + ((-1/2)) * e191
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e190 + ((-1/2)) * e191
  · linear_combination (-2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e19 + (1) * e24 + ((1/2)) * e30 + ((1/2)) * e42 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e58 + ((1/2)) * e65 + ((-1/2)) * e66 + (-1) * e90 + ((1/2)) * e91 + (1) * e104 + ((1/2)) * e116 + ((1/2)) * e118 + ((-1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e126 + (1) * e127 + ((-1/2)) * e130 + ((-1/2)) * e135 + ((-1/2)) * e136 + (-1) * e224 + (1) * e231 + (-1) * e235 + (-1) * e245 + (1) * e280 + (-1) * e293 + (-1) * e294 + (1) * e305 + (-1) * e315 + (-1) * e316 + (1) * e320
  · linear_combination ((1/2)) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e107 + ((1/2)) * e118 + ((-1/2)) * e130 + ((-1/2)) * e155 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e235 + (1) * e281
  · linear_combination ((1/2)) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e108 + ((1/2)) * e118 + ((-1/2)) * e130 + ((-1/2)) * e156 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e235 + (1) * e282
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e14 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-1) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224 + (1) * e294
  · linear_combination ((-1/2)) * e192
  · linear_combination ((-1/2)) * e193 + ((-1/2)) * e201
  · linear_combination ((-1/2)) * e194 + ((-1/2)) * e202
  · linear_combination ((-1/2)) * e195 + ((-1/2)) * e203
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + ((1/2)) * e42 + (3) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e65 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((1/2)) * e135 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e260 + (-1) * e294 + (1) * e308
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e42 + ((1/2)) * e66 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((1/2)) * e136 + ((-1/2)) * e198 + ((-1/2)) * e204 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e245 + (-1) * e283
  · linear_combination ((-1/2)) * e198 + ((-1/2)) * e204
  · linear_combination ((-1/2)) * e199 + ((-1/2)) * e205
  · linear_combination ((-1/2)) * e200 + ((-1/2)) * e206
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e14 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e60 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-1) * e127 + ((1/2)) * e131 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224 + (1) * e294
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e193 + ((-1/2)) * e201
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e194 + ((-1/2)) * e202

end
end S11D5Invariants
