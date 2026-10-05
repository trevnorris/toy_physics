import S10Anisotropic.Variation

/-! Compatibility of the anisotropic relative variational principle with finite total action. -/

namespace S10Anisotropic
open S10Pilot
noncomputable section
open MeasureTheory

variable {D : ℕ}

/-- Used below only with an explicit integrability hypothesis on the density. -/
def action (e : Fin D) (sigma : ℝ) (rho mu : ℝ) (u : Point D → Vec D) : ℝ :=
  ∫ x : Point D, lagrangian e sigma rho mu (fieldJet u x)

theorem perturbed_action_integrable (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ)
    (hbase : Integrable (fun x => lagrangian e sigma rho mu (fieldJet u x))) :
    Integrable (fun x => lagrangian e sigma rho mu (fieldJet (u + s • h) x)) := by
  convert! (relative_density_integrable e sigma hu hh rho mu s).add hbase using 1
  funext x
  simp only [Pi.add_apply, sub_add_cancel]

theorem relativeAction_eq_action_sub (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu s : ℝ)
    (hbase : Integrable (fun x => lagrangian e sigma rho mu (fieldJet u x))) :
    relativeAction e sigma rho mu u h s = action e sigma rho mu (u + s • h) - action e sigma rho mu u := by
  exact integral_sub (perturbed_action_integrable e sigma hu hh rho mu s hbase) hbase

/-- For finite total action, the derivative of S[u+s h] itself is the integral of EL·h. -/
theorem finiteAction_hasDerivAt (e : Fin D) (sigma : ℝ) {u h : Point D → Vec D}
    (hu : SmoothField u) (hh : TestField h) (rho mu : ℝ)
    (hbase : Integrable (fun x => lagrangian e sigma rho mu (fieldJet u x))) :
    HasDerivAt (fun s : ℝ => action e sigma rho mu (u + s • h))
      (∫ x, ∑ i : Fin D, eulerLagrange e sigma rho mu u x i * h x i) 0 := by
  have hd := (relativeAction_hasDerivAt e sigma hu hh rho mu).add_const (action e sigma rho mu u)
  convert! hd using 1
  · funext s
    rw [relativeAction_eq_action_sub e sigma hu hh rho mu s hbase]
    ring
  · exact (integrated_variation_by_parts e sigma hu hh rho mu).symm

end
end S10Anisotropic
