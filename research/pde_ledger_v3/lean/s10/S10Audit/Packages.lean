import S10Controls.Scalar
import S10Anisotropic.Scaling
import S10Pilot.Specialization

/-! A shared selector references the already verified supplied densities.
The stiffness operand is the exact definition used by the selected action. -/

namespace S10Audit
open S10Pilot
noncomputable section

inductive Package where
  | main | fullGradient | divergenceOnly | signFlip | anisotropic | coefficientScale
  deriving DecidableEq

def packageStiffness (p : Package) {D : ℕ} (J : Jet D) : ℝ := match p with
  | .fullGradient => S10Controls.stiffness .fullGradient J
  | .divergenceOnly => S10Controls.stiffness .divergenceOnly J
  | _ => S10Pilot.stiffness J

def packageAction (p : Package) {D : ℕ} (e : Fin D)
    (rho mu sigma scale : ℝ) (J : Jet D) : ℝ := match p with
  | .main => S10Pilot.lagrangian rho mu J
  | .fullGradient => S10Controls.lagrangian .fullGradient rho mu J
  | .divergenceOnly => S10Controls.lagrangian .divergenceOnly rho mu J
  | .signFlip => S10ScalarControls.lagrangian (-1) rho mu J
  | .anisotropic => S10Anisotropic.lagrangian e sigma rho mu J
  | .coefficientScale => S10ScalarControls.lagrangian scale rho mu J

def stiffnessCoefficient (p : Package) (mu scale : ℝ) : ℝ := match p with
  | .signFlip => -mu
  | .coefficientScale => scale * mu
  | _ => mu

def packageKinetic (p : Package) {D : ℕ} (e : Fin D) (sigma : ℝ) (J : Jet D) : ℝ :=
  if p = .anisotropic then S10Anisotropic.kinetic e sigma (J 0) else normSq (J 0)

theorem packageAction_uses_stiffness (p : Package) {D : ℕ} (e : Fin D)
    (rho mu sigma scale : ℝ) (J : Jet D) :
    packageAction p e rho mu sigma scale J =
      rho / 2 * packageKinetic p e sigma J -
        stiffnessCoefficient p mu scale / 2 * packageStiffness p J := by
  cases p <;> simp [packageAction, packageKinetic, stiffnessCoefficient, packageStiffness,
    S10Pilot.lagrangian, S10Controls.lagrangian, S10Anisotropic.lagrangian,
    S10ScalarControls.lagrangian]

end
end S10Audit
