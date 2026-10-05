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

theorem reconstructedBlock4 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 83 = 0) ∧
    (c 84 = 2 * (c 29/2)) ∧
    (c 85 = 0) ∧
    (c 86 = 0) ∧
    (c 87 = 0) ∧
    (c 88 = 0) ∧
    (c 89 = 0) ∧
    (c 90 = 0) ∧
    (c 91 = 0) ∧
    (c 92 = 0) ∧
    (c 93 = 0) ∧
    (c 94 = 1 * c 25) ∧
    (c 95 = 0) ∧
    (c 96 = 0) ∧
    (c 97 = 0) ∧
    (c 98 = 0) ∧
    (c 99 = 0) ∧
    (c 100 = 0) ∧
    (c 101 = 0) ∧
    (c 102 = 0) := by
  obtain ⟨e0,_,e2,e3,e4,_,e6,_,e8,e9,e10,_,_,_,e14,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,_,e30,_,_,_,_,_,_,e37,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,e41,e42,_,_,_,e46,_,_,_,e50,_,_,e53,_,e55,_,_,e58,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨_,e61,_,_,_,e65,e66,_,_,e69,_,_,_,_,e74,_,_,_,_,e79⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,e82,e83,e84,_,_,e87,e88,e89,e90,e91,e92,e93,_,_,e96,_,_,_⟩ := necessaryBlock4 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,e114,_,e116,e117,e118,e119⟩ := necessaryBlock5 hQ c hc
  obtain ⟨e120,e121,e122,_,e124,_,_,e127,_,_,e130,_,e132,_,_,e135,e136,_,_,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,e141,e142,e143,e144,e145,e146,e147,e148,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock7 hQ c hc
  obtain ⟨e220,e221,_,_,e224,_,_,_,_,_,_,e231,_,_,_,e235,_,_,e238,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨e240,_,_,_,e244,e245,_,_,_,_,_,e251,e252,e253,e254,e255,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,e293,e294,_,_,_,_,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨e300,e301,_,_,_,e305,e306,_,_,_,_,_,_,_,_,e315,_,_,_,_⟩ := necessaryBlock15 hQ c hc
  refine ⟨?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_,?_⟩
  · linear_combination ((-1/2)) * e79 + ((-1/2)) * e141
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((-1/2)) * e3 + ((1/2)) * e6 + ((1/2)) * e9 + ((1/2)) * e10 + ((-1/2)) * e14 + (1) * e24 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + ((-1/2)) * e118 + ((-1/2)) * e119 + ((-1/2)) * e121 + ((1/2)) * e122 + ((-1/2)) * e124 + (1) * e127 + (-1) * e220 + (1) * e221 + (-1) * e224 + (1) * e231 + (-1) * e293 + (-1) * e294 + (1) * e305
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e37 + ((1/2)) * e61 + ((-1/2)) * e82 + ((-1/2)) * e114 + ((-1/2)) * e118 + ((1/2)) * e132 + ((-1/2)) * e142 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e240 + (-1) * e252
  · linear_combination ((-1/2)) * e82 + ((-1/2)) * e142
  · linear_combination ((-1/2)) * e83 + ((-1/2)) * e143
  · linear_combination ((-1/2)) * e84 + ((-1/2)) * e144
  · linear_combination ((-1/2)) * e3 + (1) * e24 + ((-1/2)) * e41 + (1) * e46 + ((-1/2)) * e65 + ((1/2)) * e117 + ((-1/2)) * e119 + ((1/2)) * e135 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e244 + (-1) * e293 + (1) * e306
  · linear_combination ((-1/2)) * e2 + ((-1/2)) * e42 + ((1/2)) * e66 + ((-1/2)) * e87 + ((-1/2)) * e116 + ((-1/2)) * e118 + ((1/2)) * e136 + ((-1/2)) * e145 + ((-1/2)) * e220 + ((1/2)) * e221 + (1) * e245 + (-1) * e253
  · linear_combination ((-1/2)) * e87 + ((-1/2)) * e145
  · linear_combination ((-1/2)) * e88 + ((-1/2)) * e146
  · linear_combination ((-1/2)) * e89 + ((-1/2)) * e147
  · linear_combination ((1/2)) * e4 + ((1/2)) * e120 + ((1/2)) * e220 + ((-1/2)) * e221 + (1) * e293 + (1) * e315
  · linear_combination (-1) * e24 + ((1/2)) * e30 + (1) * e90 + ((-1/2)) * e91
  · linear_combination ((-1/2)) * e8 + ((-1/2)) * e92
  · linear_combination (-1) * e2 + ((-1/2)) * e30 + (1) * e46 + (-1) * e53 + ((1/2)) * e58 + (-1) * e90 + ((-1/2)) * e91 + (-1) * e118 + ((1/2)) * e130 + (-1) * e220 + (1) * e221 + (1) * e235 + (1) * e238 + ((-1/2)) * e300 + ((-1/2)) * e301
  · linear_combination ((1/2)) * e2 + ((1/2)) * e30 + ((-1/2)) * e58 + (1) * e69 + (-1) * e74 + ((1/2)) * e79 + (-1) * e90 + ((1/2)) * e91 + ((1/2)) * e118 + ((-1/2)) * e130 + ((1/2)) * e141 + ((1/2)) * e220 + ((-1/2)) * e221 + (-1) * e235 + (1) * e251
  · linear_combination ((-1/2)) * e93
  · linear_combination (-2) * e0 + ((-1/2)) * e2 + ((1/2)) * e6 + ((1/2)) * e10 + ((1/2)) * e30 + (2) * e46 + ((-1/2)) * e50 + (-1) * e55 + (1) * e90 + ((-1/2)) * e91 + ((-1/2)) * e118 + ((-1/2)) * e121 + (1) * e127 + ((-1/2)) * e220 + ((1/2)) * e221 + (-1) * e224 + (1) * e254
  · linear_combination (2) * e0 + (1) * e2 + ((-1/2)) * e6 + ((-1/2)) * e10 + ((1/2)) * e30 + (-3) * e46 + ((1/2)) * e50 + (1) * e53 + ((1/2)) * e55 + ((-1/2)) * e58 + (1) * e90 + ((1/2)) * e91 + (1) * e118 + ((1/2)) * e121 + ((-1/2)) * e127 + ((-1/2)) * e130 + (1) * e220 + (-1) * e221 + (1) * e224 + (-1) * e235 + (-1) * e238 + (-1) * e255 + ((1/2)) * e300 + ((1/2)) * e301
  · linear_combination ((-1/2)) * e96 + ((-1/2)) * e148

end
end S11D5Invariants
