import S11D5Invariants.ConstraintBlock0
import S11D5Invariants.ConstraintBlock1
import S11D5Invariants.ConstraintBlock2
import S11D5Invariants.ConstraintBlock3
import S11D5Invariants.ConstraintBlock5
import S11D5Invariants.ConstraintBlock6
import S11D5Invariants.ConstraintBlock10
import S11D5Invariants.ConstraintBlock11
import S11D5Invariants.ConstraintBlock12
import S11D5Invariants.ConstraintBlock14
import S11D5Invariants.ConstraintBlock15
import S11D5Invariants.ConstraintBlock16

/-! Bounded coefficient reconstruction using unchanged exact linear combinations. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem reconstructedBlock16 {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    (c 323 = 0) ∧
    (c 324 = 1 * (c 6/2) + 1 * (c 29/2) + 1 * c 25) := by
  obtain ⟨e0,e1,e2,_,_,e5,e6,_,_,e9,e10,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock0 hQ c hc
  obtain ⟨_,_,_,_,e24,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock1 hQ c hc
  obtain ⟨_,_,e42,_,_,_,e46,_,_,_,e50,_,_,_,_,e55,_,_,_,_⟩ := necessaryBlock2 hQ c hc
  obtain ⟨_,_,_,_,_,e65,e66,e67,e68,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock3 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e116,_,e118,_⟩ := necessaryBlock5 hQ c hc
  obtain ⟨_,e121,e122,_,_,_,_,e127,_,_,_,_,_,_,_,e135,e136,e137,e138,_⟩ := necessaryBlock6 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,e216,_,_,e219⟩ := necessaryBlock10 hQ c hc
  obtain ⟨e220,e221,e222,_,e224,_,_,_,_,_,_,_,_,_,_,_,e236,_,_,_⟩ := necessaryBlock11 hQ c hc
  obtain ⟨_,_,_,_,_,e245,_,e247,e248,_,_,_,_,_,_,_,_,_,_,_⟩ := necessaryBlock12 hQ c hc
  obtain ⟨_,_,_,_,_,_,_,_,_,_,_,_,e292,_,_,_,_,_,e298,_⟩ := necessaryBlock14 hQ c hc
  obtain ⟨_,_,e302,e303,_,_,_,_,_,_,_,_,_,_,e314,_,_,_,_,e319⟩ := necessaryBlock15 hQ c hc
  obtain ⟨_,e321⟩ := necessaryBlock16 hQ c hc
  refine ⟨?_,?_⟩
  · linear_combination (1) * e2 + ((1/2)) * e42 + (-1) * e46 + ((1/2)) * e65 + ((-1/2)) * e66 + ((-1/2)) * e67 + ((1/2)) * e116 + (1) * e118 + ((-1/2)) * e135 + ((-1/2)) * e136 + ((-1/2)) * e137 + ((1/2)) * e216 + ((-1/2)) * e219 + (1) * e220 + (-1) * e221 + (-1) * e245 + (-1) * e247 + (1) * e292 + ((-1/2)) * e302 + ((1/2)) * e303 + (1) * e314
  · linear_combination ((95/36)) * e0 + ((-2/3)) * e1 + ((-7/48)) * e2 + ((-2/3)) * e5 + ((-17/48)) * e6 + ((7/48)) * e9 + ((7/48)) * e10 + ((19/24)) * e24 + ((7/12)) * e46 + ((-31/48)) * e50 + ((-7/24)) * e55 + ((1/2)) * e68 + ((-7/48)) * e118 + ((-7/48)) * e121 + ((7/48)) * e122 + ((7/24)) * e127 + ((1/2)) * e138 + ((7/24)) * e221 + ((7/24)) * e222 + ((-7/24)) * e224 + (-1) * e236 + (1) * e248 + (1) * e298 + (1) * e319 + ((-625/288)) * e321

end
end S11D5Invariants
