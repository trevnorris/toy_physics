import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock7
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

theorem reconstructedBlock15 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 303 = 2 * (c 6/2)) ∧
    (c 304 = 1 * c 25) ∧
    (c 305 = 0) ∧
    (c 306 = 0) ∧
    (c 307 = 0) ∧
    (c 308 = 2 * (c 29/2)) ∧
    (c 309 = 0) ∧
    (c 310 = 1 * c 25) ∧
    (c 311 = 0) ∧
    (c 312 = 0) ∧
    (c 313 = 0) ∧
    (c 314 = 0) ∧
    (c 315 = 1 * c 25) ∧
    (c 316 = 0) ∧
    (c 317 = 0) ∧
    (c 318 = 0) ∧
    (c 319 = 1 * c 25) ∧
    (c 320 = 0) ∧
    (c 321 = 0) ∧
    (c 322 = 1 * c 25) := by
  obtain ⟨e0,_,e2,_,e4,_,e6,_,_,e9,e10,_,_,_,_,_,_,_,_,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,e23,e24,_,_,_,_,_,e30,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨_,_,e62,e63,_,e65,e66,e67,e68,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,e104,_,_,_,e108,_,_,_,_,_,_,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,_,_,e126,e127,_,_,e130,_,_,e133,e134,e135,e136,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e156,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,e211,e212,e213,e214,e215,e216,e217,e218,e219⟩ := necessaryBlock10 hQ c hc
  obtain ⟨e220,e221,_,e223,e224,_,_,_,_,_,_,e231,_,_,_,e235,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨_,_,e242,e243,_,e245,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,e269,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock13 hQ c hc
  obtain ⟨e280,_,e282,_,_,_,_,_,_,_,e290,e291,e292,e293,e294,e295,_,_,e298,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,e305,_,_,_,_,_,e311,e312,e313,_,e315,e316,e317,_,e319⟩ := necessaryBlock15 hQ c hc
  obtain ⟨e320,_⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination (-1) * e0 + (1) * e6 + (1) * e23 + (1) * e50 + ((-1/2)) * e62 + ((-1/2)) * e68 + ((-1/2)) * e133 + ((-1/2)) * e138 + (1) * e223 + (2) * e236 + (-1) * e242 + (-1) * e248 + (1) * e269 + (1) * e295 + (-2) * e298 + (1) * e311 + (1) * e317 + (-1) * e319
  · linear_combination (-1) * e2 + ((1/2)) * e4 + ((-1/2)) * e30 + ((1/2)) * e58 + ((1/2)) * e63 + (1) * e90 + ((-1/2)) * e91 + (-1) * e118 + ((1/2)) * e120 + ((1/2)) * e130 + ((1/2)) * e134 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e235 + (1) * e243 + (1) * e293 + ((1/2)) * e300 + ((-1/2)) * e301 + (1) * e315
  · linear_combination ((-1/2)) * e211 + ((-1/2)) * e212
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e211 + ((-1/2)) * e212
  · linear_combination ((1/2)) * e2 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e211 + ((-1/2)) * e212 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e245 + (1) * e290
  · linear_combination (-2) * e0 + ((3/2)) * e2 + ((-1/2)) * e4 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e19 + (1) * e24 + ((1/2)) * e30 + ((1/2)) * e42 + (1) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e58 + ((-1/2)) * e63 + ((1/2)) * e65 + ((-1/2)) * e66 + ((-1/2)) * e67 + (-1) * e90 + ((1/2)) * e91 + (1) * e104 + ((1/2)) * e116 + ((3/2)) * e118 + ((-1/2)) * e120 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e126 + (1) * e127 + ((-1/2)) * e130 + ((-1/2)) * e134 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((-1/2)) * e137 + (1) * e220 + (-1) * e221 + (-1) * e224 + (1) * e231 + (-1) * e235 + (-1) * e243 + (-1) * e245 + (-1) * e247 + (1) * e280 + (-1) * e293 + (-1) * e294 + ((-1/2)) * e300 + ((1/2)) * e301 + ((-1/2)) * e302 + ((1/2)) * e303 + (1) * e305 + (1) * e312 + (-1) * e315 + (-1) * e316 + (1) * e320
  · linear_combination (1) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + ((-1/2)) * e63 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e108 + (1) * e118 + ((-1/2)) * e130 + ((-1/2)) * e134 + ((-1/2)) * e156 + (1) * e220 + (-1) * e221 + (-1) * e235 + (-1) * e243 + (1) * e282 + ((-1/2)) * e300 + ((1/2)) * e301 + (1) * e313
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e19 + (-2) * e46 + ((1/2)) * e50 + (1) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e126 + (-1) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224 + (1) * e294 + (1) * e316
  · linear_combination ((-1/2)) * e213
  · linear_combination ((-1/2)) * e214 + ((-1/2)) * e217
  · linear_combination ((-1/2)) * e215 + ((-1/2)) * e218
  · linear_combination ((-1/2)) * e216 + ((-1/2)) * e219
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e19 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e65 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e126 + (-1) * e127 + ((1/2)) * e135 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224 + (1) * e294 + (1) * e316
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e214 + ((-1/2)) * e217
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e215 + ((-1/2)) * e218
  · linear_combination (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e135 + ((1/2)) * e216 + ((-1/2)) * e219
  · linear_combination (2) * e0 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e19 + ((-1/2)) * e42 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e65 + ((1/2)) * e66 + ((-1/2)) * e116 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e126 + (-1) * e127 + ((1/2)) * e135 + ((1/2)) * e136 + (1) * e224 + (1) * e245 + (1) * e294 + (1) * e316
  · linear_combination ((1/2)) * e2 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e215 + ((-1/2)) * e218 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e245 + (1) * e291
  · linear_combination ((1/2)) * e2 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e216 + ((-1/2)) * e219 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e245 + (1) * e292
  · linear_combination (2) * e0 + ((-1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e19 + ((-1/2)) * e42 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e65 + ((1/2)) * e66 + ((1/2)) * e67 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e126 + (-1) * e127 + ((1/2)) * e135 + ((1/2)) * e136 + ((1/2)) * e137 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e224 + (1) * e245 + (1) * e247 + (1) * e294 + ((1/2)) * e302 + ((-1/2)) * e303 + (1) * e316

end
end S11D5Invariants
