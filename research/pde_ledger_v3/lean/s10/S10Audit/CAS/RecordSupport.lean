import S10Audit.CAS.CoincidenceBindings

namespace S10Audit.CAS
open S10Pilot S10Anisotropic Polynomial
noncomputable section
set_option backward.isDefEq.respectTransparency false

/-- What the CAS reported, including an explicitly unresolved sign. -/
inductive ReportedSign where
  | zero | positive | negative | undecided
  deriving DecidableEq

def ReportedSign.Holds (s : ReportedSign) (x : ℝ) : Prop :=
  match s with
  | .zero => x = 0
  | .positive => 0 < x
  | .negative => x < 0
  | .undecided => True

def classifiedSign (i : Fin 3) : ReportedSign :=
  if i = 0 then .zero else .positive

theorem reference_root_sign (rho mu sigma : ℝ) (k : Vec 3) (i : Fin 3)
    (hr : 0 < rho) (hm : 0 < mu) (hs : 0 < sigma) (hk : k ≠ 0) :
    (classifiedSign i).Holds (referenceRoot rho mu sigma k i) := by
  fin_cases i
  · rfl
  · change 0 < (mu / rho) * normSq k
    exact mul_pos (div_pos hm hr) (dot_self_pos hk)
  · change 0 < (mu / rho) * extraValue 0 sigma k
    exact mul_pos (div_pos hm hr) (extraValue_pos hs hk)

theorem reference_root_isRoot (rho mu sigma : ℝ) (k : Vec 3) (i : Fin 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k i) k) = 0 := by
  rw [rootPolynomial_complete rho mu sigma _ k hr hs, rootPolynomial_roots rho mu sigma k hr hs]
  fin_cases i <;> simp [referenceRoot]

end
end S10Audit.CAS
