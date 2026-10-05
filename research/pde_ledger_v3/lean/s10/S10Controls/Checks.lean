import S10Controls.Certificate

/-! Concrete density distinctions and the dimension-one coincidence. -/

namespace S10Controls
open S10Pilot
noncomputable section

/-- Both controls reduce to the same scalar-gradient stiffness in one dimension. -/
theorem dimension_one_stiffness_eq (J : Jet 1) :
    stiffness .fullGradient J = stiffness .divergenceOnly J := by
  simp [stiffness, fullGradientStiffness, divergenceOnlyStiffness, divergence]

/-- One common derivative jet distinguishes all three supplied densities in D=2. -/
theorem concrete_distinct_densities :
    let J : Jet 2 := ![![0, 0], ![1, 0], ![0, 1]]
    lagrangian .fullGradient 1 1 J = -1 ∧
    lagrangian .divergenceOnly 1 1 J = -2 ∧
    S10Pilot.lagrangian 1 1 J = 0 := by
  norm_num [lagrangian, stiffness, fullGradientStiffness, divergenceOnlyStiffness, divergence,
    S10Pilot.lagrangian, S10Pilot.stiffness, antisym, normSq, dot, Fin.sum_univ_succ]

end
end S10Controls
