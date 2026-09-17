import Mathlib.LinearAlgebra.Dimension.Finrank
import Mathlib.Tactic

/-! NP1: typed modal algebra. Analytic inverse existence is not assumed proved. -/
namespace S11NonlinearPole
noncomputable section

variable {𝕜 X Y K : Type*} [Field 𝕜]
  [AddCommGroup X] [Module 𝕜 X] [AddCommGroup Y] [Module 𝕜 Y]
  [AddCommGroup K] [Module 𝕜 K]

/-- The whole modal pairing is invertible. Full kernel coverage is an additional
hypothesis of `residue_of_inverse_coefficients`, not built into this algebra. -/
structure PairingData (𝕜 X Y K : Type*) [Field 𝕜]
    [AddCommGroup X] [Module 𝕜 X] [AddCommGroup Y] [Module 𝕜 Y]
    [AddCommGroup K] [Module 𝕜 K] where
  derivative : X →ₗ[𝕜] Y
  right : K →ₗ[𝕜] X
  left : Y →ₗ[𝕜] K
  pairing : K ≃ₗ[𝕜] K
  pairing_eq : ∀ k, left (derivative (right k)) = pairing k

namespace PairingData

variable (d : PairingData 𝕜 X Y K)

def residueCandidate : Y →ₗ[𝕜] X :=
  d.right.comp (d.pairing.symm.toLinearMap.comp d.left)

def fieldProjection : X →ₗ[𝕜] X := d.residueCandidate.comp d.derivative

def sourceProjection : Y →ₗ[𝕜] Y := d.derivative.comp d.residueCandidate

theorem right_injective : Function.Injective d.right := by
  intro x y h
  apply d.pairing.injective
  simpa only [d.pairing_eq] using congrArg (fun v => d.left (d.derivative v)) h

theorem derivative_right_injective : Function.Injective (d.derivative.comp d.right) := by
  intro x y h
  apply d.pairing.injective
  simpa only [LinearMap.comp_apply, d.pairing_eq] using congrArg d.left h

theorem residue_sandwich :
    d.residueCandidate.comp (d.derivative.comp d.residueCandidate) = d.residueCandidate := by
  ext y
  change d.right (d.pairing.symm (d.left (d.derivative
    (d.right (d.pairing.symm (d.left y)))))) = d.right (d.pairing.symm (d.left y))
  rw [d.pairing_eq, d.pairing.symm_apply_apply]

theorem fieldProjection_fixes (k : K) : d.fieldProjection (d.right k) = d.right k := by
  change d.right (d.pairing.symm (d.left (d.derivative (d.right k)))) = d.right k
  rw [d.pairing_eq, d.pairing.symm_apply_apply]

theorem sourceProjection_fixes (k : K) :
    d.sourceProjection (d.derivative (d.right k)) = d.derivative (d.right k) := by
  change d.derivative (d.right (d.pairing.symm (d.left (d.derivative (d.right k))))) = _
  rw [d.pairing_eq, d.pairing.symm_apply_apply]

theorem fieldProjection_idempotent :
    d.fieldProjection.comp d.fieldProjection = d.fieldProjection := by
  ext x
  exact LinearMap.congr_fun d.residue_sandwich (d.derivative x)

theorem sourceProjection_idempotent :
    d.sourceProjection.comp d.sourceProjection = d.sourceProjection := by
  ext y
  exact congrArg d.derivative (LinearMap.congr_fun d.residue_sandwich y)

theorem range_fieldProjection : LinearMap.range d.fieldProjection = LinearMap.range d.right := by
  ext x
  constructor
  · rintro ⟨y, rfl⟩
    exact ⟨d.pairing.symm (d.left (d.derivative y)), rfl⟩
  · rintro ⟨k, rfl⟩
    exact ⟨d.right k, d.fieldProjection_fixes k⟩

theorem range_sourceProjection :
    LinearMap.range d.sourceProjection = LinearMap.range (d.derivative.comp d.right) := by
  ext y
  constructor
  · rintro ⟨x, rfl⟩
    exact ⟨d.pairing.symm (d.left x), rfl⟩
  · rintro ⟨k, rfl⟩
    exact ⟨d.derivative (d.right k), d.sourceProjection_fixes k⟩

theorem finrank_fieldProjection [FiniteDimensional 𝕜 K] :
    Module.finrank 𝕜 (LinearMap.range d.fieldProjection) = Module.finrank 𝕜 K := by
  rw [d.range_fieldProjection]
  exact LinearMap.finrank_range_of_inj d.right_injective

theorem finrank_sourceProjection [FiniteDimensional 𝕜 K] :
    Module.finrank 𝕜 (LinearMap.range d.sourceProjection) = Module.finrank 𝕜 K := by
  rw [d.range_sourceProjection]
  exact LinearMap.finrank_range_of_inj d.derivative_right_injective

theorem residue_unique (R : Y →ₗ[𝕜] X)
    (hRange : LinearMap.range R ≤ LinearMap.range d.right)
    (hNormalization : ∀ y, d.left (d.derivative (R y)) = d.left y) :
    R = d.residueCandidate := by
  ext y
  obtain ⟨k, hk⟩ := hRange ⟨y, rfl⟩
  have hd : d.pairing k = d.left y := by
    rw [← d.pairing_eq, hk]
    exact hNormalization y
  change R y = d.right (d.pairing.symm (d.left y))
  rw [← hk, ← hd, d.pairing.symm_apply_apply]

/-- These are the leading and constant equations from an actual simple inverse
expansion. This theorem identifies its coefficient; it does not construct the
analytic expansion or assert a semisimplicity criterion for arbitrary pencils. -/
theorem residue_of_inverse_coefficients (L₀ : X →ₗ[𝕜] Y) (R H : Y →ₗ[𝕜] X)
    (hFull : LinearMap.range d.right = LinearMap.ker L₀)
    (hLeft : ∀ x, d.left (L₀ x) = 0)
    (hLeading : ∀ y, L₀ (R y) = 0)
    (hConstant : ∀ y, L₀ (H y) + d.derivative (R y) = y) :
    R = d.residueCandidate := by
  apply d.residue_unique R
  · intro x hx
    rw [hFull]
    obtain ⟨y, rfl⟩ := hx
    exact hLeading y
  · intro y
    simpa only [map_add, hLeft, zero_add] using congrArg d.left (hConstant y)

end PairingData
end
end S11NonlinearPole
