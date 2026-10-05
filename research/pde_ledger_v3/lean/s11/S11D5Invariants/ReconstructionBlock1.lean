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

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock1 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 21 = 0) ∧
    (c 22 = 0) ∧
    (c 23 = 0) ∧
    (c 24 = 2 * (c 6/2)) ∧
    (c 26 = 0) ∧
    (c 27 = 0) ∧
    (c 28 = 0) ∧
    (c 30 = 0) ∧
    (c 31 = 0) ∧
    (c 32 = 0) ∧
    (c 33 = 0) ∧
    (c 34 = 0) ∧
    (c 35 = 0) ∧
    (c 36 = 0) ∧
    (c 37 = 0) ∧
    (c 38 = 0) ∧
    (c 39 = 0) ∧
    (c 40 = 0) ∧
    (c 41 = 0) ∧
    (c 42 = 0) := by
  obtain ⟨e0,_,e2,e3,e4,e5,e6,_,_,e9,e10,_,_,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨e20,_,_,_,e24,e25,e26,e27,e28,e29,e30,e31,e32,_,_,_,e36,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,e49,e50,_,_,_,_,e55,e56,e57,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,_,_,e65,e66,e67,e68,e69,_,e71,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,e109,e110,e111,e112,e113,e114,e115,e116,_,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,e125,_,e127,e128,e129,e130,e131,e132,e133,_,e135,e136,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e220,e221,e222,e223,e224,_,_,_,e228,_,_,_,e232,_,e234,e235,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,e242,_,_,e245,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,e295,_,e297,e298,e299⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,_,_,_,_,_,_,_,e315,_,e317,e318,e319⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination (1) * e0 + ((-1/2)) * e20 + (-1) * e46 + ((1/2)) * e65 + ((1/2)) * e125 + ((-1/2)) * e135
  · linear_combination (1) * e0 + ((1/2)) * e2 + ((-1/2)) * e20 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((1/2)) * e125 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e228 + (-1) * e245
  · linear_combination (1) * e0 + (1) * e2 + ((-1/2)) * e20 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((-1/2)) * e67 + ((1/2)) * e116 + (1) * e118 + ((1/2)) * e125 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((-1/2)) * e137 + (1) * e220 + (-1) * e221 + (1) * e228 + (-1) * e245 + (-1) * e247 + (1) * e297 + ((-1/2)) * e302 + ((1/2)) * e303
  · linear_combination ((1/2)) * e6 + ((1/2)) * e50 + ((-1/2)) * e68 + ((-1/2)) * e138 + (1) * e223 + (1) * e236 + (-1) * e248 + (1) * e295 + (-1) * e298 + (1) * e317 + (-1) * e319
  · linear_combination (1) * e24 + ((-1/2)) * e25 + (1) * e46 + ((1/2)) * e109
  · linear_combination (1) * e24 + ((-1/2)) * e26 + (1) * e69 + ((1/2)) * e110
  · linear_combination (1) * e24 + ((-1/2)) * e27 + (1) * e90 + ((1/2)) * e111
  · linear_combination (-1) * e0 + ((-1/2)) * e2 + (-1) * e5 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((1/2)) * e122 + (1) * e127 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e222 + (-1) * e224
  · linear_combination ((-1/2)) * e28 + ((-1/2)) * e49
  · linear_combination ((-1/2)) * e29 + ((-1/2)) * e71
  · linear_combination ((-1/2)) * e30 + ((-1/2)) * e91
  · linear_combination (1) * e24 + ((-1/2)) * e31 + (1) * e46 + ((-1/2)) * e55 + ((1/2)) * e113 + ((1/2)) * e127
  · linear_combination ((-1/2)) * e32 + ((-1/2)) * e112
  · linear_combination (1) * e0 + ((3/2)) * e2 + (1) * e5 + (-1) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + (-2) * e46 + (1) * e55 + ((-1/2)) * e56 + ((3/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e122 + (-1) * e127 + ((-1/2)) * e128 + ((3/2)) * e220 + ((-3/2)) * e221 + (-1) * e222 + (1) * e224 + (-1) * e232 + (-1) * e236
  · linear_combination ((-1/2)) * e3 + ((1/2)) * e4 + ((1/2)) * e29 + ((-1/2)) * e30 + ((-1/2)) * e57 + ((1/2)) * e58 + (-1) * e69 + ((1/2)) * e71 + (1) * e90 + ((-1/2)) * e91 + ((-1/2)) * e119 + ((1/2)) * e120 + ((-1/2)) * e129 + ((1/2)) * e130 + (-1) * e234 + (1) * e235 + ((-1/2)) * e300 + ((-1/2)) * e301 + (1) * e315 + (-1) * e318
  · linear_combination ((-1/2)) * e300 + ((-1/2)) * e301
  · linear_combination (1) * e24 + ((-1/2)) * e36 + (1) * e46 + ((-1/2)) * e60 + ((1/2)) * e115 + ((1/2)) * e131
  · linear_combination ((-1/2)) * e37 + ((-1/2)) * e114
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e14 + ((1/2)) * e30 + ((-1/2)) * e37 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e58 + ((-1/2)) * e60 + ((1/2)) * e61 + (-1) * e90 + ((1/2)) * e91 + ((-1/2)) * e114 + ((1/2)) * e118 + ((-1/2)) * e120 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-1) * e127 + ((-1/2)) * e130 + ((1/2)) * e131 + ((1/2)) * e132 + (1) * e224 + (-1) * e235 + (1) * e240 + (-1) * e293 + (1) * e294 + (-1) * e299 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e62 + ((1/2)) * e118 + ((-1/2)) * e133 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e242

end
end S11D5Invariants
