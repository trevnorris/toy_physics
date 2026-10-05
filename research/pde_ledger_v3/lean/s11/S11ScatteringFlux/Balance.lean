import S11ScatteringFlux.Current

/-! F3: orientation, explicit zero-denominator exclusion and conditional balance.
No physical current conservation or positivity is assumed silently. -/
namespace S11ScatteringFlux
noncomputable section

inductive End where | left | right
  deriving DecidableEq

def outward : End → ℝ | .left => -1 | .right => 1
def incident (e : End) (j : ℝ) : ℝ := -outward e * j
def outgoing (left right : ℝ) : ℝ := outward .left * left + outward .right * right
def fraction (num den : ℝ) : Option ℝ := if den = 0 then none else some (num / den)

theorem end_coverage (e : End) : e = .left ∨ e = .right := by cases e <;> simp
theorem flux_sign_coverage (j : ℝ) : j < 0 ∨ j = 0 ∨ 0 < j := lt_trichotomy j 0
theorem incident_left (j : ℝ) : incident .left j = j := by simp [incident, outward]
theorem incident_right (j : ℝ) : incident .right j = -j := by simp [incident, outward]
theorem outgoing_eq (l r : ℝ) : outgoing l r = -l + r := by simp [outgoing, outward]
theorem fraction_undefined_iff (num den : ℝ) : fraction num den = none ↔ den = 0 := by
  simp [fraction]
theorem fraction_defined (num den : ℝ) (h : den ≠ 0) :
    fraction num den = some (num / den) := by simp [fraction, h]

theorem fraction_nonneg {num den : ℝ} (hn : 0 ≤ num) (hd : 0 < den) :
    0 ≤ num / den := div_nonneg hn hd.le
theorem fraction_le_one_iff {num den : ℝ} (hd : 0 < den) :
    num / den ≤ 1 ↔ num ≤ den := by simpa using (div_le_one hd)

/-- The supplied defect includes any unaccounted physical current or source term.
The premise is an application obligation, not a theorem about an arbitrary S matrix. -/
theorem conditional_balance {converted survived defect incoming : ℝ}
    (h : converted + survived + defect = incoming) (hd : incoming ≠ 0) :
    converted / incoming + survived / incoming = 1 - defect / incoming := by
  calc
    converted / incoming + survived / incoming = (converted + survived) / incoming :=
      (add_div _ _ _).symm
    _ = (incoming - defect) / incoming := by congr 1; linarith
    _ = 1 - defect / incoming := by rw [sub_div, div_self hd]

theorem conservation_requires_zero_defect {converted survived defect incoming : ℝ}
    (h : converted + survived + defect = incoming) (hd : incoming ≠ 0) :
    converted / incoming + survived / incoming = 1 ↔ defect = 0 := by
  rw [conditional_balance h hd]
  constructor
  · intro hh
    have hz : defect / incoming = 0 := by linarith
    exact (div_eq_zero_iff.mp hz).resolve_right hd
  · intro hh
    simp [hh]

end
end S11ScatteringFlux
