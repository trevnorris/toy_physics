import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock4
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock11
import S11D5Invariants.ConstraintBlock12
import S11D5Invariants.ConstraintBlock15

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock2 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 43 = 0) ∧
    (c 44 = 0) ∧
    (c 45 = 0) ∧
    (c 46 = 0) ∧
    (c 47 = 0) ∧
    (c 48 = 0) ∧
    (c 49 = 1 * c 25) ∧
    (c 50 = 0) ∧
    (c 51 = 0) ∧
    (c 52 = 0) ∧
    (c 53 = 0) ∧
    (c 54 = 0) ∧
    (c 55 = 0) ∧
    (c 56 = 0) ∧
    (c 57 = 2 * (c 29/2)) ∧
    (c 58 = 0) ∧
    (c 59 = 0) ∧
    (c 60 = 0) ∧
    (c 61 = 0) ∧
    (c 62 = 0) := by
  obtain ⟨e0,_,e2,e3,e4,_,e6,_,_,_,e10,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,e26,e27,e28,_,e30,_,_,_,_,_,e36,_,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,e41,e42,_,_,_,e46,_,_,e49,e50,e51,_,_,_,e55,e56,e57,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨e60,_,_,e63,_,e65,_,e67,e68,e69,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e90,e91,_,_,_,_,_,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,e110,e111,_,_,_,e115,e116,e117,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,_,_,_,_,_,e127,e128,e129,e130,e131,_,_,e134,e135,_,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,e229,e230,e231,_,e233,_,e235,_,e237,e238,e239⟩ := necessaryBlock11 hQ c hc
  obtain ⟨_,_,_,e243,_,_,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨e300,e301,e302,e303,_,_,_,_,_,_,_,_,_,_,_,e315,_,_,e318,_⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e63 + ((1/2)) * e118 + ((-1/2)) * e134 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e243
  · linear_combination (1) * e24 + ((-1/2)) * e41 + (1) * e46 + ((-1/2)) * e65 + ((1/2)) * e117 + ((1/2)) * e135
  · linear_combination ((-1/2)) * e42 + ((-1/2)) * e116
  · linear_combination ((-1/2)) * e302 + ((-1/2)) * e303
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e67 + ((1/2)) * e118 + ((-1/2)) * e137 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e247
  · linear_combination ((1/2)) * e2 + ((-1/2)) * e68 + ((1/2)) * e118 + ((-1/2)) * e138 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e248
  · linear_combination ((1/2)) * e2 + ((1/2)) * e118 + ((1/2)) * e220 + ((-1/2)) * e221
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e26 + (1) * e69 + ((1/2)) * e110 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e229
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e27 + (1) * e90 + ((1/2)) * e111 + ((-1/2)) * e118 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e230
  · linear_combination (-1) * e24 + ((1/2)) * e28 + (1) * e46 + ((-1/2)) * e49
  · linear_combination ((-1/2)) * e6 + ((-1/2)) * e50
  · linear_combination ((-1/2)) * e51
  · linear_combination (1) * e2 + ((1/2)) * e3 + ((-1/2)) * e4 + ((1/2)) * e30 + ((-1/2)) * e58 + (1) * e69 + (-1) * e90 + ((1/2)) * e91 + (1) * e118 + ((1/2)) * e119 + ((-1/2)) * e120 + ((-1/2)) * e130 + (1) * e220 + (-1) * e221 + (-1) * e235 + (-1) * e237 + ((1/2)) * e300 + ((1/2)) * e301 + (-1) * e315 + (1) * e318
  · linear_combination (1) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + ((1/2)) * e91 + (1) * e118 + ((-1/2)) * e130 + (1) * e220 + (-1) * e221 + (-1) * e235 + (-1) * e238 + ((1/2)) * e300 + ((1/2)) * e301
  · linear_combination (-2) * e0 + (-1) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + (1) * e24 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + (-1) * e118 + ((-1/2)) * e121 + (1) * e127 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e231
  · linear_combination (2) * e0 + ((1/2)) * e2 + ((-1/2)) * e6 + ((-1/2)) * e10 + ((1/2)) * e28 + (-2) * e46 + ((1/2)) * e49 + ((1/2)) * e50 + ((1/2)) * e55 + ((1/2)) * e118 + ((1/2)) * e121 + ((-1/2)) * e127 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e224 + (-1) * e233
  · linear_combination ((-1/2)) * e56 + ((-1/2)) * e128
  · linear_combination ((-1/2)) * e57 + ((-1/2)) * e129
  · linear_combination ((-1/2)) * e58 + ((-1/2)) * e130
  · linear_combination ((-1/2)) * e2 + (1) * e24 + ((-1/2)) * e36 + (1) * e46 + ((-1/2)) * e60 + ((1/2)) * e115 + ((-1/2)) * e118 + ((1/2)) * e131 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e239

end
end S11D5Invariants
