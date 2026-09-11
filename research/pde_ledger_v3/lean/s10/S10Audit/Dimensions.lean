import Physlib.Units.Dimension
import Mathlib.Tactic

/-! Q6 expression-tree dimensions in PhysLean's dimensional algebra.
Sums are nonempty. Numeric literals infer dimensionless units by convention;
this does not give zero a unique physical dimension. No simplification of an expression precedes its dimension check. -/

namespace S10Audit
noncomputable section

instance : DimensionBasis (Fin 3) := DimensionBasis.pi _
abbrev Dim := Dimension (Fin 3)

def dimensions (l t m : ℚ) : Dim := Dimension.ofFunction
  ![Dimension.Exponent.ofRat l, Dimension.Exponent.ofRat t, Dimension.Exponent.ofRat m]

def lengthDim : Dim := dimensions 1 0 0
def timeDim : Dim := dimensions 0 1 0
def energyDim : Dim := dimensions 2 (-2) 1
def densityDim (D : ℕ) : Dim := energyDim / lengthDim ^ D

inductive Expr (α : Type) where
  | atom (a : α)
  | scalar (value : ℝ)
  | add (a b : Expr α)
  | sub (a b : Expr α)
  | mul (a b : Expr α)
  | div (a b : Expr α)
  | pow (a : Expr α) (n : ℕ)
  | sum (n : ℕ) (terms : Fin (n + 1) → Expr α)

namespace Expr
variable {α : Type}

def eval (v : α → ℝ) : Expr α → ℝ
  | .atom a => v a
  | .scalar c => c
  | .add a b => eval v a + eval v b
  | .sub a b => eval v a - eval v b
  | .mul a b => eval v a * eval v b
  | .div a b => eval v a / eval v b
  | .pow a n => eval v a ^ n
  | .sum _ f => ∑ i, eval v (f i)

def infer (u : α → Dim) : Expr α → Option Dim
  | .atom a => some (u a)
  | .scalar _ => some 1
  | .add a b | .sub a b => do
      let x ← infer u a
      let y ← infer u b
      if x = y then some x else none
  | .mul a b => do return (← infer u a) * (← infer u b)
  | .div a b => do return (← infer u a) / (← infer u b)
  | .pow a n => do return (← infer u a) ^ n
  | .sum _ f => do
      let x ← infer u (f 0)
      if ∀ i, infer u (f i) = some x then some x else none

inductive HasDim (u : α → Dim) : Expr α → Dim → Prop where
  | atom (a) : HasDim u (.atom a) (u a)
  | scalar (c) : HasDim u (.scalar c) 1
  | add {a b d} : HasDim u a d → HasDim u b d → HasDim u (.add a b) d
  | sub {a b d} : HasDim u a d → HasDim u b d → HasDim u (.sub a b) d
  | mul {a b d e} : HasDim u a d → HasDim u b e → HasDim u (.mul a b) (d * e)
  | div {a b d e} : HasDim u a d → HasDim u b e → HasDim u (.div a b) (d / e)
  | pow {a d} (n) : HasDim u a d → HasDim u (.pow a n) (d ^ n)
  | sum {n f d} : (∀ i, HasDim u (f i) d) → HasDim u (.sum n f) d

theorem HasDim.infer {u : α → Dim} {e : Expr α} {d : Dim} (h : HasDim u e d) :
    infer u e = some d := by
  induction h with
  | atom => rfl
  | scalar => rfl
  | add _ _ ha hb => simp [Expr.infer, ha, hb]
  | sub _ _ ha hb => simp [Expr.infer, ha, hb]
  | mul _ _ ha hb => simp [Expr.infer, ha, hb]
  | div _ _ ha hb => simp [Expr.infer, ha, hb]
  | pow n _ ha => simp [Expr.infer, ha]
  | sum _ ha => simp [Expr.infer, ha]

/-- Dimensional typing implies covariance under every multiplicative unit change.
The value of every atom, including a bare field, receives its assigned unit. -/
theorem HasDim.rescale {u : α → Dim} {e : Expr α} {d : Dim} (h : HasDim u e d)
    (χ : Dim →* ℝˣ) (v : α → ℝ) :
    eval (fun a => (χ (u a) : ℝ) * v a) e = (χ d : ℝ) * eval v e := by
  induction h with
  | atom => rfl
  | scalar => simp [eval]
  | add _ _ ha hb => simp only [eval, ha, hb, mul_add]
  | sub _ _ ha hb => simp only [eval, ha, hb, mul_sub]
  | mul _ _ ha hb => simp only [eval, ha, hb, map_mul, Units.val_mul]; ring
  | div _ _ ha hb =>
      simp only [eval, ha, hb, map_div, Units.val_div_eq_div_val]
      exact mul_div_mul_comm _ _ _ _
  | pow n _ ha => simp only [eval, ha, map_pow, Units.val_pow_eq_pow_val, mul_pow]
  | sum _ ha => simp only [eval, ha, Finset.mul_sum]

@[simp] theorem infer_sum_iff {u : α → Dim} {n : ℕ} {f : Fin (n + 1) → Expr α} {d : Dim} :
    infer u (.sum n f) = some d ↔ ∀ i, infer u (f i) = some d := by
  change ((infer u (f 0)).bind fun x =>
    if ∀ i, infer u (f i) = some x then some x else none) = some d ↔ _
  cases h : infer u (f 0) with
  | none => simp only [Option.bind_none, reduceCtorEq, false_iff]; intro hi; simpa [h] using hi 0
  | some x =>
      simp only [Option.bind_some]
      split_ifs with hx
      · simp only [Option.some.injEq]
        constructor
        · rintro rfl; exact hx
        · intro hi; exact Option.some.inj (h.symm.trans (hi 0))
      · simp only [false_iff]
        intro hi
        have hxd : x = d := Option.some.inj (h.symm.trans (hi 0))
        exact hx (hxd.symm ▸ hi)

theorem infer_sound {u : α → Dim} {e : Expr α} {d : Dim}
    (h : infer u e = some d) : HasDim u e d := by
  induction e generalizing d with
  | atom a =>
      have hd : u a = d := Option.some.inj h
      exact hd ▸ HasDim.atom a
  | scalar c =>
      have hd : (1 : Dim) = d := Option.some.inj h
      exact hd ▸ HasDim.scalar c
  | add a b ha hb =>
      change ((infer u a).bind fun x => (infer u b).bind fun y =>
        if x = y then some x else none) = some d at h
      obtain ⟨x, hx, h⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨y, hy, h⟩ := Option.bind_eq_some_iff.mp h
      split_ifs at h with hxy
      · have hd := Option.some.inj h
        exact hd ▸ HasDim.add (ha hx) (hxy.symm ▸ hb hy)
  | sub a b ha hb =>
      change ((infer u a).bind fun x => (infer u b).bind fun y =>
        if x = y then some x else none) = some d at h
      obtain ⟨x, hx, h⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨y, hy, h⟩ := Option.bind_eq_some_iff.mp h
      split_ifs at h with hxy
      · have hd := Option.some.inj h
        exact hd ▸ HasDim.sub (ha hx) (hxy.symm ▸ hb hy)
  | mul a b ha hb =>
      change ((infer u a).bind fun x => (infer u b).bind fun y => some (x * y)) = some d at h
      obtain ⟨x, hx, h⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨y, hy, h⟩ := Option.bind_eq_some_iff.mp h
      exact Option.some.inj h ▸ HasDim.mul (ha hx) (hb hy)
  | div a b ha hb =>
      change ((infer u a).bind fun x => (infer u b).bind fun y => some (x / y)) = some d at h
      obtain ⟨x, hx, h⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨y, hy, h⟩ := Option.bind_eq_some_iff.mp h
      exact Option.some.inj h ▸ HasDim.div (ha hx) (hb hy)
  | pow a m ha =>
      change ((infer u a).bind fun x => some (x ^ m)) = some d at h
      obtain ⟨x, hx, h⟩ := Option.bind_eq_some_iff.mp h
      exact Option.some.inj h ▸ HasDim.pow m (ha hx)
  | sum m f hf => exact HasDim.sum fun i => hf i (infer_sum_iff.mp h i)

theorem infer_iff_hasDim {u : α → Dim} {e : Expr α} {d : Dim} :
    infer u e = some d ↔ HasDim u e d := ⟨infer_sound, HasDim.infer⟩

end Expr
end
end S10Audit
