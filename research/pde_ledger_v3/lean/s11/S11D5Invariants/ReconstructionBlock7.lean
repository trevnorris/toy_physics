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

theorem reconstructedBlock7 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 143 = 0) ∧
    (c 144 = 0) ∧
    (c 145 = 0) ∧
    (c 146 = 0) ∧
    (c 147 = 2 * (c 6/2)) ∧
    (c 148 = 0) ∧
    (c 149 = 0) ∧
    (c 150 = 0) ∧
    (c 151 = 0) ∧
    (c 152 = 0) ∧
    (c 153 = 2 * (c 6/2)) ∧
    (c 154 = 1 * c 25) ∧
    (c 155 = 0) ∧
    (c 156 = 0) ∧
    (c 157 = 0) ∧
    (c 158 = 2 * (c 29/2)) ∧
    (c 159 = 0) ∧
    (c 160 = 0) ∧
    (c 161 = 0) ∧
    (c 162 = 0) := by
  obtain ⟨e0,_,e2,_,e4,_,e6,_,e8,e9,e10,_,_,e13,e14,e15,e16,e17,e18,e19⟩ := necessaryBlock0 hQ c hc
  obtain ⟨e20,e21,e22,e23,e24,_,e26,e27,e28,_,e30,_,_,_,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,e47,e48,e49,e50,_,_,_,e54,e55,e56,e57,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,e61,e62,e63,_,e65,e66,e67,e68,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,e92,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e110,e111,_,_,e114,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,e123,e124,e125,e126,e127,e128,e129,e130,e131,e132,e133,e134,e135,e136,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e220,e221,_,e223,e224,_,e226,e227,e228,e229,e230,e231,_,e233,_,e235,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,e241,e242,e243,_,e245,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,e295,e296,e297,e298,e299⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,_,_,_,_,_,_,_,e315,_,e317,e318,e319⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e8 + (1) * e13 + ((1/2)) * e30 + ((-1/2)) * e58 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e92 + ((1/2)) * e118 + ((-1/2)) * e130 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e226 + (-1) * e235
  · linear_combination ((-1/2)) * e15 + ((-1/2)) * e123
  · linear_combination (-1) * e0 + ((1/2)) * e14 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e124 + ((-1/2)) * e131
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e15 + (1) * e16 + ((1/2)) * e37 + (-1) * e46 + ((1/2)) * e60 + ((-1/2)) * e61 + ((1/2)) * e114 + ((1/2)) * e118 + ((1/2)) * e123 + ((-1/2)) * e131 + ((-1/2)) * e132 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e227 + (-1) * e240
  · linear_combination (-1) * e0 + ((1/2)) * e6 + (1) * e17 + ((1/2)) * e50 + ((-1/2)) * e62 + ((-1/2)) * e133 + (1) * e223 + (1) * e236 + (-1) * e242 + (1) * e295 + (-1) * e298
  · linear_combination (1) * e2 + ((-1/2)) * e8 + (1) * e18 + ((1/2)) * e30 + ((-1/2)) * e58 + ((-1/2)) * e63 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e92 + (1) * e118 + ((-1/2)) * e130 + ((-1/2)) * e134 + (1) * e220 + (-1) * e221 + (1) * e226 + (-1) * e235 + (-1) * e243 + (1) * e296 + ((-1/2)) * e300 + ((1/2)) * e301
  · linear_combination ((-1/2)) * e20 + ((-1/2)) * e125
  · linear_combination (-1) * e0 + ((1/2)) * e19 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e126 + ((-1/2)) * e135
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e20 + (1) * e21 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((1/2)) * e116 + ((1/2)) * e118 + ((1/2)) * e125 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e228 + (-1) * e245
  · linear_combination (1) * e2 + ((-1/2)) * e20 + (1) * e22 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((-1/2)) * e67 + ((1/2)) * e116 + (1) * e118 + ((1/2)) * e125 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((-1/2)) * e137 + (1) * e220 + (-1) * e221 + (1) * e228 + (-1) * e245 + (-1) * e247 + (1) * e297 + ((-1/2)) * e302 + ((1/2)) * e303
  · linear_combination (-1) * e0 + ((1/2)) * e6 + (1) * e23 + ((1/2)) * e50 + ((-1/2)) * e68 + ((-1/2)) * e138 + (1) * e223 + (1) * e236 + (-1) * e248 + (1) * e295 + (-1) * e298 + (1) * e317 + (-1) * e319
  · linear_combination ((1/2)) * e2 + (1) * e46 + ((1/2)) * e118 + ((1/2)) * e220 + ((-1/2)) * e221
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e26 + (-1) * e46 + (1) * e47 + ((1/2)) * e110 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e229
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e27 + (-1) * e46 + (1) * e48 + ((1/2)) * e111 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e230
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((-1/2)) * e28 + (2) * e46 + ((-1/2)) * e49 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e118 + ((-1/2)) * e121 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e233
  · linear_combination (-2) * e0 + (-1) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + (1) * e24 + ((-1/2)) * e50 + (1) * e54 + ((-1/2)) * e55 + (-1) * e118 + ((-1/2)) * e121 + ((1/2)) * e127 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e231
  · linear_combination (-1) * e46 + ((1/2)) * e56 + ((-1/2)) * e128
  · linear_combination (-1) * e46 + ((1/2)) * e57 + ((-1/2)) * e129
  · linear_combination (-1) * e46 + ((1/2)) * e58 + ((-1/2)) * e130
  · linear_combination (2) * e0 + ((-1/2)) * e2 + ((-1/2)) * e4 + ((-1/2)) * e6 + ((-1/2)) * e9 + ((-1/2)) * e10 + ((1/2)) * e14 + ((1/2)) * e30 + (-1) * e37 + (-1) * e46 + ((1/2)) * e50 + (1) * e55 + ((-1/2)) * e58 + (-1) * e60 + (1) * e61 + (-1) * e90 + ((1/2)) * e91 + (-1) * e114 + ((-1/2)) * e118 + ((-1/2)) * e120 + ((1/2)) * e121 + ((-1/2)) * e122 + ((1/2)) * e124 + (-1) * e127 + ((-1/2)) * e130 + (1) * e132 + (-1) * e220 + (1) * e221 + (1) * e224 + (-1) * e235 + (2) * e240 + (1) * e241 + (-1) * e293 + (1) * e294 + (-1) * e299 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318

end
end S11D5Invariants
