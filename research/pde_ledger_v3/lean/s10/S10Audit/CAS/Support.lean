import S10Audit.NormalizedTrees

/-! Reference semantics for the first transcript bridge: one-axis anisotropic
inertia in D=3. Squared frequency is an independent real variable. -/

namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section

inductive Symbol where
  | rho | mu | sigma | z | k (i : Fin 3)
  deriving DecidableEq

def values (rho mu sigma z : ℝ) (k : Vec 3) : Symbol → ℝ
  | .rho => rho
  | .mu => mu
  | .sigma => sigma
  | .z => z
  | .k i => k i

def units : Symbol → Dim
  | .rho => dimensions (-3) 0 1
  | .mu => dimensions (-1) (-2) 1
  | .sigma => dimensions 0 0 0
  | .z => dimensions 0 (-2) 0
  | .k _ => dimensions (-1) 0 0

theorem rho_units : units .rho = inferredRho 3 lengthDim := by
  simpa [units] using (inferred_rho_LTM 3).symm

theorem mu_units : units .mu = inferredMu 3 lengthDim := by
  have h := (inferred_mu_LTM 3).symm
  norm_num at h
  exact h

theorem dimensions_zero : dimensions 0 0 0 = 1 := by
  ext i
  fin_cases i <;> rfl

theorem dimensions_mul (a b c d e f : ℚ) :
    dimensions a b c * dimensions d e f = dimensions (a + d) (b + e) (c + f) := by
  ext i
  apply Dimension.Exponent.ringEquivRat.injective
  fin_cases i <;> simp [dimensions, Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat]

theorem dimensions_div (a b c d e f : ℚ) :
    dimensions a b c / dimensions d e f = dimensions (a - d) (b - e) (c - f) := by
  ext i
  apply Dimension.Exponent.ringEquivRat.injective
  fin_cases i <;> simp [dimensions, Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat]

theorem dimensions_pow (a b c : ℚ) (n : ℕ) :
    dimensions a b c ^ n = dimensions (n * a) (n * b) (n * c) := by
  have hn : ((n : Dimension.Exponent) : ℚ) = (n : ℚ) :=
    map_natCast Dimension.Exponent.ringEquivRat n
  ext i
  apply Dimension.Exponent.ringEquivRat.injective
  fin_cases i <;>
    simp [dimensions, Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat, nsmul_eq_mul, hn]

theorem castDim {u : Symbol → Dim} {e : Expr Symbol} {d d' : Dim}
    (h : Expr.HasDim u e d) (hd : d = d') : Expr.HasDim u e d' := hd ▸ h

/-- This polynomial matrix has no coefficient or wavevector denominators. -/
def referenceMatrix (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => rho * z * (if i = 0 then sigma else 1) * (if i = j then 1 else 0) -
    mu * (normSq k * (if i = j then 1 else 0) - k i * k j)

theorem referenceMatrix_action (rho mu sigma omega : ℝ) (k : Vec 3) :
    referenceMatrix rho mu sigma (omega ^ 2) k =
      actionMatrix .anisotropic 0 rho mu sigma 1 omega k := by
  ext i j
  simp only [referenceMatrix, actionMatrix, packageOperator, S10Anisotropic.modalOperator,
    S10Pilot.modalOperator, Pi.add_apply, Pi.smul_apply, smul_eq_mul, dot_unit_right]
  by_cases hij : i = j <;> by_cases hi : i = 0 <;> by_cases hj : j = 0 <;>
    simp_all [unit, eq_comm]
  all_goals ring

theorem referenceMatrix_normalized (rho mu sigma z : ℝ) (k a : Vec 3) (hm : mu ≠ 0) :
    (referenceMatrix rho mu sigma z k).mulVec a =
      mu • normalizedOperator 0 sigma (rho * z / mu) k a := by
  ext i
  fin_cases i <;>
    simp [referenceMatrix, Matrix.mulVec, dotProduct, normalizedOperator,
      unit, normSq, dot, Fin.sum_univ_succ] <;> field_simp <;> ring

def referenceRoot (rho mu sigma : ℝ) (k : Vec 3) : Fin 3 → ℝ
  | 0 => 0
  | 1 => coneValue rho mu k
  | _ => extraConeValue 0 sigma rho mu k

def referenceBasis (sigma : ℝ) (k : Vec 3) : Fin 3 → Vec 3
  | 0 => (k 2)⁻¹ • k
  | 1 => chartVector k 1 2
  | _ => (-(sigma * k 0 * k 2))⁻¹ • extraVector 0 sigma k

/-- Basis denominators are chart restrictions, separate from physical strata. -/
def BasisDomain (sigma : ℝ) (k : Vec 3) : Fin 3 → Prop
  | 0 => k 2 ≠ 0
  | 1 => k 1 ≠ 0
  | _ => sigma ≠ 0 ∧ k 0 ≠ 0 ∧ k 2 ≠ 0

theorem referenceBasis_last (sigma : ℝ) (k : Vec 3) (r : Fin 3) (h : BasisDomain sigma k r) :
    referenceBasis sigma k r 2 = 1 := by
  fin_cases r
  · change (k 2)⁻¹ * k 2 = 1
    exact inv_mul_cancel₀ h
  · simp [referenceBasis, chartVector, unit]
  · rcases h with ⟨hs, h0, h2⟩
    simp [referenceBasis, extraVector, unit]
    field_simp

open Lean Elab Tactic in
/-- A proof-producing tactic for the generated rational identities. Every
operation produces an ordinary kernel-checked proof term. -/
elab "cas_equal" : tactic => do
  if (← getGoals).isEmpty then return
  evalTactic (← `(tactic| norm_num [referenceMatrix, referenceRoot, referenceBasis,
    Matrix.det_fin_three, Matrix.mulVec, dotProduct, chartVector, extraVector,
    extraNumerator, extraConeValue, extraValue, perpSq, coneValue, unit, normSq,
    dot, Fin.sum_univ_three, Fin.ext_iff, Fin.coe_ofNat_eq_mod, -Fin.val_eq_zero_iff]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| field_simp))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| ring_nf))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| norm_num))

open Lean Elab Tactic in
elab "cas_eval" h:Lean.Parser.Tactic.rwRule : tactic => do
  evalTactic (← `(tactic| rw [$h]))
  if !(← getGoals).isEmpty then
    evalTactic (← `(tactic| cas_equal))

end
end S10Audit.CAS
