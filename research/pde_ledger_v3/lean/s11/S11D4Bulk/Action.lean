import S11D4Bulk.Even
import S11D4Odd.Boundary

/-! D4C.1: the full classified density and derivative-defined momentum.
The even and odd pieces use the same gradient, normalization and spacetime. -/
namespace S11D4Bulk
noncomputable section
open S10Pilot
open scoped ContDiff

def lagrangian (v : Vec 4) (J : Jet 4) : ℝ :=
  Even.lagrangian v J + S11D4Odd.lagrangian (v 3) J

theorem density_identity (v : Vec 4) (J : Jet 4) :
    lagrangian v J = -(1 / 2 : ℝ) *
      S11D4Invariants.invariantForm v (Even.spatialGradient J) := by
  rw [S11D4Invariants.invariantForm_apply]
  simp only [lagrangian, Even.lagrangian, Even.divergence, Even.transposePair,
    Even.gradientSquare, Even.spatialGradient, S11D4Odd.lagrangian,
    S11D4Odd.orientationJet,
    Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Fin.sum_univ_four]
  ring!

theorem all_invariant_densities (Q : S11D4Invariants.Quad)
    (hQ : S11D4Invariants.SOInvariant Q) :
    ∃! v : Vec 4, Q = S11D4Invariants.invariantForm v :=
  S11D4Invariants.SO_unique Q hQ

def variationDensity (v : Vec 4) (J H : Jet 4) : ℝ :=
  Even.variationDensity v J H + S11D4Odd.variationDensity (v 3) J H

theorem lagrangian_increment (v : Vec 4) (s : ℝ) (J H : Jet 4) :
    lagrangian v (J + s • H) = lagrangian v J +
      s * variationDensity v J H + s ^ 2 * lagrangian v H := by
  simp only [lagrangian, variationDensity, Even.lagrangian_increment,
    S11D4Odd.lagrangian_increment]
  ring

theorem lagrangian_variation (v : Vec 4) (J H : Jet 4) :
    HasDerivAt (fun s : ℝ => lagrangian v (J + s • H))
      (variationDensity v J H) 0 :=
  (Even.lagrangian_variation v J H).add (S11D4Odd.lagrangian_variation (v 3) J H)

def momentum (v : Vec 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) : ℝ :=
  deriv (fun s : ℝ => lagrangian v (J + s • basisJet j i)) 0

theorem momentum_eq (v : Vec 4) (J : Jet 4) (j : Fin 5) (i : Fin 4) :
    momentum v J j i = Even.momentum v J j i + S11D4Odd.momentum (v 3) J j i := by
  unfold momentum Even.momentum S11D4Odd.momentum
  rw [(lagrangian_variation v J (basisJet j i)).deriv,
    (Even.lagrangian_variation v J (basisJet j i)).deriv,
    (S11D4Odd.lagrangian_variation (v 3) J (basisJet j i)).deriv]
  rfl

theorem momentum_time (v : Vec 4) (J : Jet 4) (i : Fin 4) : momentum v J 0 i = 0 := by
  simp [momentum_eq, Even.momentum_eq, S11D4Odd.momentum_eq]

def eulerLagrange (v : Vec 4) (u : Point 4 → Vec 4) (x : Point 4) : Vec 4 :=
  fun i => -∑ j : Fin 5, coordDeriv j (fun y => momentum v (fieldJet u y) j i) x

theorem smooth_momentum {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (j : Fin 5) (i : Fin 4) :
    ContDiff ℝ ∞ (fun x => momentum v (fieldJet u x) j i) := by
  simp only [momentum_eq]
  exact (Even.smooth_momentum hu v j i).add (S11D4Odd.smooth_momentum hu (v 3) j i)

theorem smooth_eulerLagrange {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (i : Fin 4) : ContDiff ℝ ∞ (fun x => eulerLagrange v u x i) := by
  have hd : ∀ j, ContDiff ℝ ∞ (coordDeriv j (fun y => momentum v (fieldJet u y) j i)) :=
    fun j => smooth_coordDeriv (smooth_momentum hu v j i) j
  unfold eulerLagrange
  fun_prop

theorem eulerLagrange_even {u : Point 4 → Vec 4} (hu : SmoothField u)
    (v : Vec 4) (x : Point 4) : eulerLagrange v u x = Even.eulerLagrange v u x := by
  have hs : eulerLagrange v u x =
      Even.eulerLagrange v u x + S11D4Odd.eulerLagrange (v 3) u x := by
    ext i
    simp only [Pi.add_apply, eulerLagrange, Even.eulerLagrange, S11D4Odd.eulerLagrange, momentum_eq]
    simp_rw [S11D4Odd.partial_add (Even.smooth_momentum hu v _ i)
      (S11D4Odd.smooth_momentum hu (v 3) _ i)]
    simp [Finset.sum_add_distrib, add_comm]
  rw [hs, S11D4Odd.eulerLagrange_zero hu, add_zero]

end
end S11D4Bulk
