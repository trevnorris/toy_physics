# Moving-Throat PDE — Stage 20: Exact Continuum-Kernel Extraction of `A`, `DeltaK_ax`, `beta_0`, and `M_mix`

## Purpose

Stage 19 reduced the selected-branch theorem gate to the exact admissibility geometry
driven by

- the normalization ratio `R_target`,
- the mixed baseline `M_mix`,
- and the bare anisotropy ratio `delta`.

But those quantities were still being treated as microscopic inputs.

The right next move is therefore to derive them from the first explicit **linearized finite-throat continuum operator** that is still consistent with the reduced hierarchy already used in Stages 10–19.

That is what this stage does.

The main result is that the first continuum operator already produces all of the key reduced quantities in closed form. After exact mode projection and mass normalization, one finds

`A = [ K_U (K_eta + 6 T_Omega) - c_(eta U)^2 ] / (mu_eta K_U),`

`DeltaK_ax = pi^2 T_w / (mu_eta L^2),`

`beta_0 = (mu_W / mu_eta) ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`         / [ ( K_U K_W^(eff) - c_(UW)^2 sigma )^2 ],`

`alpha_mix = ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`            / [ mu_eta K_U ( K_U K_W^(eff) - c_(UW)^2 sigma ) ],`

`M_mix = 8 alpha_mix / (pi^2 A),`

with

`K_W^(eff) := K_W + pi^2 T_W / (4 L^2),`
`sigma = 88 / (9 pi^2).`

So the Stage-17/19 branch variables are no longer abstract. They are exact low-mode functionals of one explicit continuum kernel.

---

## 1. Minimal linearized continuum operator

Work on the finite throat interval

`s in [0,L]`.

Keep the same wall/support ontology already used in the reduced stages:

- a quadrupole wall field `eta(s,t)` with N/N boundaries,
- a brane-like internal field `U(s,t)` on the same interval,
- a support/BdG field `phi(s,t)` on the D/N half-wave branch,
- a mixed `W(s,t)` field representing the `A_w/F_(mu w)/J^w` lane on the same D/N branch.

The minimal quadratic Lagrangian density is taken to be

`L_eta = (mu_eta/2) dot(eta)^2 - (T_w/2) (eta')^2 - (K_eta + 6 T_Omega) eta^2 / 2,`

`L_U   = (mu_U/2) dot(U)^2 - K_U U^2 / 2,`

`L_phi = (mu_phi/2) dot(phi)^2 - (T_phi/2) (phi')^2 - K_phi phi^2 / 2,`

`L_W   = (mu_W/2) dot(W)^2 - (T_W/2) (W')^2 - K_W W^2 / 2,`

with local bilinear couplings

`L_(eta U)   = - c_(eta U) int_0^L eta U ds,`

`L_(eta phi) = - c_(eta phi) int_0^L eta phi ds,`

`L_(eta W)   = - c_(eta W) int_0^L eta W ds,`

`L_(UW)      = + c_(UW) int_0^L U W ds.`

Two comments matter.

First, the wall operator is exactly the same quadrupole wall operator already carried from the earlier stages.

Second, the internal `U` field is kept in the same **brane-like flat-doublet limit** already used implicitly in Stages 12–16: at this minimal order there is no axial-gradient penalty in `U`, so the first two N/N channels stay degenerate.
That is the smallest continuum choice that reproduces the earlier reduced hierarchy.

---

## 2. Exact mode basis and overlaps

Use the exact first two N/N modes

`u_0(s) = 1 / sqrt(L),`

`u_1(s) = sqrt(2/L) cos(pi s / L),`

and the exact lowest D/N half-wave

`f_0(s) = sqrt(2/L) sin(pi s / (2L)).`

The exact D/N overlap vector is

`v_i = int_0^L u_i(s) f_0(s) ds,`

so

`kappa_0 = 2 sqrt(2) / pi,`

`kappa_1 = - 4 / (3 pi),`

`v = (kappa_0, kappa_1)^T,`

`sigma = v.v = 88 / (9 pi^2).`

This is the same overlap data already driving the Stage-10/15 profile and source maps, but now it is being used directly as the projection datum of a continuum operator.

---

## 3. Mass-normalized modal coordinates

Expand

`eta(s,t) = q_0(t) u_0(s) + q_1(t) u_1(s),`

`U(s,t)   = u_0^(int)(t) u_0(s) + u_1^(int)(t) u_1(s),`

`phi(s,t) = phi_0(t) f_0(s),`

`W(s,t)   = W_0(t) f_0(s).`

Now pass to mass-normalized modal coordinates

`Q_i   = sqrt(mu_eta) q_i,`

`U_i   = sqrt(mu_U) u_i^(int),`

`Phi   = sqrt(mu_phi) phi_0,`

`Wbar  = sqrt(mu_W) W_0.`

Then the reduced bare kernels are

`K_0 = (K_eta + 6 T_Omega) / mu_eta,`

`DeltaK_ax = pi^2 T_w / (mu_eta L^2),`

`varpi^2 = ( K_phi + pi^2 T_phi / (4 L^2) ) / mu_phi,`

`Omega_U^2 = K_U / mu_U,`

`Omega_W^2 = ( K_W + pi^2 T_W / (4 L^2) ) / mu_W`
`          = K_W^(eff) / mu_W.`

The mass-normalized couplings are

`g_U = c_(eta U) / sqrt(mu_eta mu_U),`

`g_B = c_(eta phi) / sqrt(mu_eta mu_phi),`

`g_W = c_(eta W) / sqrt(mu_eta mu_W),`

`g_R = c_(UW) / sqrt(mu_U mu_W).`

So the projected couplings become

`L_(eta U)   = - g_U Q.U,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(eta W)   = - g_W (v.Q) Wbar,`

`L_(UW)      = + g_R (v.U) Wbar.`

That is exactly the reduced pattern used abstractly in Stage 12. The difference is that every symbol in it has now been tied back to the continuum operator.

---

## 4. Exact Schur complement of the internal block

In frequency space the internal kernels are

`A_U(omega)   = Omega_U^2 - omega^2,`

`A_phi(omega) = varpi^2 - omega^2,`

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega).`

Eliminating `(U, Phi, Wbar)` gives the exact wall self-energy

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T,`

with

`Xi(omega) = g_U^2 / A_U(omega),`

`alpha(omega) = g_B^2 / A_phi(omega)`
`             + ( A_U(omega) g_W + g_R g_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ],`

`Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma.`

So the Stage-11 rank-1 loading law is not phenomenological.
It is the exact first Schur complement of the continuum operator.

On the conservative static branch `Pi_out(0)=0`, define

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma,`

`Chi = Omega_U^2 g_W + g_R g_U.`

Then the key reduced quantities are

`A = K_0 - g_U^2 / Omega_U^2,`

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 ),`

`beta_0 = Chi^2 / Delta_0^2,`

`M_mix = 8 alpha_mix / (pi^2 A),`

`delta = DeltaK_ax / A.`

This is exactly the reduced data set that was still abstract at Stage 19.

---

## 5. Closed continuum formulas

Now substitute the mass-normalized projections back into the definitions above.

Introduce

`K_eta^(eff) := K_eta + 6 T_Omega,`

`K_W^(eff)   := K_W + pi^2 T_W / (4 L^2).`

Then

`A = [ K_U K_eta^(eff) - c_(eta U)^2 ] / ( mu_eta K_U ),`

`DeltaK_ax = pi^2 T_w / ( mu_eta L^2 ),`

`Delta_0 = [ K_U K_W^(eff) - c_(UW)^2 sigma ] / ( mu_U mu_W ),`

`Chi = [ K_U c_(eta W) + c_(UW) c_(eta U) ]`
`      / [ mu_U sqrt(mu_eta mu_W) ].`

So the static mixed loading is

`alpha_mix`
`= ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ mu_eta K_U ( K_U K_W^(eff) - c_(UW)^2 sigma ) ],`

the outgoing transfer factor is

`beta_0`
`= (mu_W / mu_eta)`
`  ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ ( K_U K_W^(eff) - c_(UW)^2 sigma )^2 ],`

and the dimensionless mixed baseline is

`M_mix`
`= 8 ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ pi^2 ( K_U K_eta^(eff) - c_(eta U)^2 )`
`          ( K_U K_W^(eff) - c_(UW)^2 sigma ) ].`

Finally, the bare anisotropy ratio is

`delta`
`= pi^2 T_w K_U`
`  / [ L^2 ( K_U K_eta^(eff) - c_(eta U)^2 ) ].`

So the branch data entering Stage 17–19 are now explicit continuum-kernel functions.

---

## 6. Exact continuum stability inequalities

The selected-branch stability gates now become exact inequalities on the continuum operator.

### 6.1 Wall/internal flat-branch stability

`A > 0`
is equivalent to

`K_U K_eta^(eff) > c_(eta U)^2.`

So the wall/brane-like internal coupling cannot over-soften the flat wall branch.

### 6.2 Internal mixed-sector positivity

`Delta_0 > 0`
is equivalent to

`K_U K_W^(eff) > c_(UW)^2 sigma.`

So the internal `U/W` block must stay below its continuum overmixing threshold.

These are the first exact continuum-kernel positivity surfaces that underwrite the Stage-16/19 reduced branch geometry.

---

## 7. Best current theorem gate after Stage 20

The open theorem gap has now narrowed again.

The selected quadrupole branch is no longer parameterized by abstract reduced inputs.
The first explicit finite-throat continuum operator already determines

- `A`,
- `DeltaK_ax`,
- `alpha_mix`,
- `beta_0`,
- `M_mix`,
- and therefore `delta`.

So the remaining task is no longer “invent microscopic inputs.”
It is:

> determine whether the completed moving-throat PDE places the actual defect on the admissible Stage-18/19 branch geometry generated by these continuum-kernel quantities.
