import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
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

theorem reconstructedBlock14 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 283 = 0) ∧
    (c 284 = 0) ∧
    (c 285 = 0) ∧
    (c 286 = 0) ∧
    (c 287 = 0) ∧
    (c 288 = 0) ∧
    (c 289 = 1 * c 25) ∧
    (c 290 = 0) ∧
    (c 291 = 0) ∧
    (c 292 = 0) ∧
    (c 293 = 0) ∧
    (c 294 = 0) ∧
    (c 295 = 0) ∧
    (c 296 = 0) ∧
    (c 297 = 1 * (c 6/2) + 1 * (c 29/2) + 1 * c 25) ∧
    (c 298 = 0) ∧
    (c 299 = 0) ∧
    (c 300 = 0) ∧
    (c 301 = 0) ∧
    (c 302 = 0) := by
  obtain ⟨e0,e1,e2,_,e4,e5,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,_,e30,_,_,_,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,_,e65,e66,e67,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,e114,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,e126,e127,_,_,e130,e131,e132,e133,e134,e135,e136,e137,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,e194,e195,e196,e197,e198,e199⟩ := necessaryBlock9 hQ c hc
  obtain ⟨e200,_,e202,e203,e204,e205,e206,e207,e208,e209,e210,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock10 hQ c hc
  obtain ⟨e220,e221,e222,_,e224,_,_,_,_,_,_,_,_,_,_,e235,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,e242,e243,_,e245,_,e247,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨e260,_,_,e263,_,_,_,e267,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock13 hQ c hc
  obtain ⟨_,_,_,e283,e284,e285,e286,e287,e288,e289,_,_,_,_,e294,_,_,_,e298,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,e308,e309,e310,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock15 hQ c hc
  obtain ⟨_,e321⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e195 + ((-1/2)) * e203
  · linear_combination ((1/2)) * e2 + ((1/2)) * e42 + ((-1/2)) * e60 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((1/2)) * e131 + ((-1/2)) * e135 + ((-1/2)) * e136 + (-1) * e197 + ((1/2)) * e198 + ((1/2)) * e204 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e245 + (1) * e283
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + ((1/2)) * e42 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((1/2)) * e60 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((-1/2)) * e131 + (1) * e196 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e260 + (-1) * e294 + (1) * e308
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e198 + ((-1/2)) * e204
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e199 + ((-1/2)) * e205
  · linear_combination (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e131 + ((1/2)) * e200 + ((-1/2)) * e206
  · linear_combination (2) * e0 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e14 + ((-1/2)) * e37 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e60 + ((1/2)) * e61 + ((-1/2)) * e114 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-1) * e127 + ((1/2)) * e131 + ((1/2)) * e132 + (1) * e224 + (1) * e240 + (1) * e294
  · linear_combination ((1/2)) * e2 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e194 + ((-1/2)) * e202 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e240 + (1) * e284
  · linear_combination ((1/2)) * e2 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e195 + ((-1/2)) * e203 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e240 + (1) * e285
  · linear_combination ((-1/2)) * e207 + ((-1/2)) * e208
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e207 + ((-1/2)) * e208
  · linear_combination (-2) * e0 + ((1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + ((1/2)) * e37 + (1) * e42 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((1/2)) * e60 + ((-1/2)) * e61 + ((-1/2)) * e66 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((-1/2)) * e136 + (1) * e196 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e224 + (-1) * e240 + (-1) * e245 + (1) * e260 + (1) * e286 + (-1) * e294 + (1) * e308
  · linear_combination ((1/2)) * e2 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e199 + ((-1/2)) * e205 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e240 + (1) * e287
  · linear_combination ((1/2)) * e2 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e200 + ((-1/2)) * e206 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e240 + (1) * e288
  · linear_combination ((95/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e5 + ((-17/48)) * e6 + ((7/48)) * e9 + ((7/48)) * e10 + ((19/24)) * e24 + ((7/12)) * e46 + ((-31/48)) * e50 + ((-7/24)) * e55 + ((1/2)) * e62 + ((-7/48)) * e118 + ((-7/48)) * e121 + ((7/48)) * e122 + ((7/24)) * e127 + ((1/2)) * e133 + ((7/24)) * e221 + ((7/24)) * e222 + ((-7/24)) * e224 + (-1) * e236 + (1) * e242 + (1) * e298 + ((-625/288)) * e321
  · linear_combination (-1) * e0 + (1) * e2 + ((1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e30 + ((1/2)) * e50 + ((-1/2)) * e58 + ((-1/2)) * e62 + ((-1/2)) * e63 + (-1) * e90 + ((1/2)) * e91 + (1) * e118 + ((-1/2)) * e120 + ((-1/2)) * e130 + ((-1/2)) * e133 + ((-1/2)) * e134 + (1) * e220 + (-1) * e221 + (-1) * e235 + (1) * e236 + (-1) * e242 + (-1) * e243 + (1) * e263 + (-1) * e298 + ((-1/2)) * e300 + ((1/2)) * e301 + (1) * e309
  · linear_combination ((-1/2)) * e209 + ((-1/2)) * e210
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e209 + ((-1/2)) * e210
  · linear_combination ((1/2)) * e2 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e209 + ((-1/2)) * e210 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e245 + (1) * e289
  · linear_combination (-1) * e0 + (1) * e2 + ((1/2)) * e6 + ((1/2)) * e19 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e50 + ((-1/2)) * e62 + ((1/2)) * e65 + ((-1/2)) * e66 + ((-1/2)) * e67 + ((1/2)) * e116 + (1) * e118 + ((-1/2)) * e126 + ((-1/2)) * e133 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((-1/2)) * e137 + (1) * e220 + (-1) * e221 + (1) * e236 + (-1) * e242 + (-1) * e245 + (-1) * e247 + (1) * e267 + (-1) * e298 + ((-1/2)) * e302 + ((1/2)) * e303 + (1) * e310

end
end S11D5Invariants
