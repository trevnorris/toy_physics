(* ===================================================================== *)
(*  ORR–SOMMERFELD FROM 2D VORTICITY EQN + ROBUST EQUIVALENCE CHECK       *)
(* ===================================================================== *)

ClearAll["Global`*"];

(* Assumptions go right here *)
$Assumptions = Element[{k, nu, x, y, t}, Reals] && k > 0 && nu >= 0;
simp[expr_] := FullSimplify[expr, Assumptions -> $Assumptions];

phase = Exp[gamma*t + I*k*y];

(* Base shear flow: v0 = V0(x) y-hat *)
v0 = {0, V0[x], 0};

(* Streamfunction perturbation *)
psi1 = phi[x]*phase;

(* Incompressible 2D perturbation velocity: v1 = (∂y psi, -∂x psi, 0) *)
v1 = {D[psi1, y], -D[psi1, x], 0};

(* Total velocity *)
vTot = v0 + eps*v1;

ωz[v_] := D[v[[2]], x] - D[v[[1]], y];
lap2[f_] := D[f, {x, 2}] + D[f, {y, 2}];

ωTot = ωz[vTot];

(* 2D vorticity equation *)
eqω = D[ωTot, t] + vTot[[1]]*D[ωTot, x] + vTot[[2]]*D[ωTot, y] - nu*lap2[ωTot];

(* Linearize in eps and divide out Fourier phase *)
lin = Expand[Coefficient[eqω, eps, 1]];
odeDerived = simp[Expand[lin/phase]];

Print["--- DERIVED LINEAR ODE (LHS == 0) ---"];
Print[odeDerived == 0];

(* ----------------------------- *)
(* Build OSstandard in 3 steps   *)
(* (this prevents truncation)    *)
(* ----------------------------- *)
OSstandard = 0;

OSstandard += (gamma + I*k*V0[x])*(D[phi[x], {x, 2}] - k^2*phi[x]);
OSstandard += -I*k*D[V0[x], {x, 2}]*phi[x];
OSstandard += -nu*(D[phi[x], {x, 4}] - 2*k^2*D[phi[x], {x, 2}] + k^4*phi[x]);

Print["\n--- OSstandard (must contain V0'' and nu blocks) ---"];
Print[OSstandard];

OSexpandedMinus = Expand[-OSstandard];
OSexpandedPlus  = Expand[ OSstandard];

Print["\n--- CHECK A: derived + OSstandard (should be 0) ---"];
Print[simp[Expand[odeDerived + OSstandard]]];

Print["\n--- CHECK B: derived - (-OSstandard expanded) (should be 0) ---"];
Print[simp[Expand[odeDerived - OSexpandedMinus]]];

Print["\n--- COLLECTED DERIVED ---"];
Print[Collect[Expand[odeDerived], {D[phi[x], {x, 4}], D[phi[x], {x, 2}], phi[x]}]];

Print["\n--- COLLECTED OSstandard ---"];
Print[Collect[Expand[OSstandard], {D[phi[x], {x, 4}], D[phi[x], {x, 2}], phi[x]}]];

(*"
Output:

--- DERIVED LINEAR ODE (LHS == 0) ---
-((gamma + 2*k^2*nu + I*k*V0[x])*Derivative[2][phi][x]) + k*phi[x]*(k*(gamma + k^2*nu + I*k*V0[x]) + I*Derivative[2][V0][x]) + nu*Derivative[4][phi][x] == 0

--- OSstandard (must contain V0'' and nu blocks) ---
(gamma + I*k*V0[x])*(-(k^2*phi[x]) + Derivative[2][phi][x]) - I*k*phi[x]*Derivative[2][V0][x] - nu*(k^4*phi[x] - 2*k^2*Derivative[2][phi][x] + Derivative[4][phi][x])

--- CHECK A: derived + OSstandard (should be 0) ---
0

--- CHECK B: derived - (-OSstandard expanded) (should be 0) ---
0

--- COLLECTED DERIVED ---
(-gamma - 2*k^2*nu - I*k*V0[x])*Derivative[2][phi][x] + phi[x]*(gamma*k^2 + k^4*nu + I*k^3*V0[x] + I*k*Derivative[2][V0][x]) + nu*Derivative[4][phi][x]

--- COLLECTED OSstandard ---
(gamma + 2*k^2*nu + I*k*V0[x])*Derivative[2][phi][x] + phi[x]*(-(gamma*k^2) - k^4*nu - I*k^3*V0[x] - I*k*Derivative[2][V0][x]) - nu*Derivative[4][phi][x]
"*)
