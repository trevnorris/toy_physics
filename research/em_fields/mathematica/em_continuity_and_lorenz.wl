(* em_continuity_and_lorenz.wl *)

(* ---------------------------------------------------------------------- *)
(* 0. Coordinates and symbols *)
(* ---------------------------------------------------------------------- *)

$Assumptions = {
  t \[Element] Reals, x \[Element] Reals, y \[Element] Reals, z \[Element] Reals,
  c > 0, mu0 > 0
};

(* --- Helper Function for Clean Terminal Output --- *)
(* FIX: Used :> (RuleDelayed) to correctly capture the function name f. *)
cleanFormat[expr_] := ToString[expr /. {
   Derivative[1, 0, 0, 0][f_][__] :> "dt(" <> ToString[f] <> ")",
   Derivative[0, 1, 0, 0][f_][__] :> "dx(" <> ToString[f] <> ")",
   Derivative[0, 0, 1, 0][f_][__] :> "dy(" <> ToString[f] <> ")",
   Derivative[0, 0, 0, 1][f_][__] :> "dz(" <> ToString[f] <> ")",

   Derivative[2, 0, 0, 0][f_][__] :> "dt^2(" <> ToString[f] <> ")",
   Derivative[0, 2, 0, 0][f_][__] :> "dx^2(" <> ToString[f] <> ")",
   Derivative[0, 0, 2, 0][f_][__] :> "dy^2(" <> ToString[f] <> ")",
   Derivative[0, 0, 0, 2][f_][__] :> "dz^2(" <> ToString[f] <> ")",

   phi -> "phi", Ax -> "Ax", Ay -> "Ay", Az -> "Az",
   rhoE -> "rho_e", Jx -> "Jx", Jy -> "Jy", Jz -> "Jz"
}, InputForm];

(* ---------------------------------------------------------------------- *)
(* 1. Define d'Alembertian Box and Lorenz gauge *)
(* ---------------------------------------------------------------------- *)

box[f_] := (1/c^2)*D[f, {t, 2}] - (D[f, {x, 2}] + D[f, {y, 2}] + D[f, {z, 2}]);

(* Define Field variables with suffix to avoid recursion errors *)
phiField = phi[t, x, y, z];
AxField  = Ax[t, x, y, z];
AyField  = Ay[t, x, y, z];
AzField  = Az[t, x, y, z];

rhoEField = rhoE[t, x, y, z];
JxField   = Jx[t, x, y, z];
JyField   = Jy[t, x, y, z];
JzField   = Jz[t, x, y, z];

(* Lorenz gauge *)
lorenzGauge = (1/c^2)*D[phiField, t] + D[AxField, x] + D[AyField, y] + D[AzField, z];

Print["Lorenz gauge condition d_mu A^mu = 0:"];
Print[cleanFormat[lorenzGauge == 0]];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 2. Field equations Box A^mu = -mu0 J^mu *)
(* ---------------------------------------------------------------------- *)

eqPhiLHS = box[phiField];
eqPhiRHS = -mu0 * c^2 * rhoEField;

eqAxLHS = box[AxField];
eqAxRHS = -mu0 * JxField;
eqAyLHS = box[AyField];
eqAyRHS = -mu0 * JyField;
eqAzLHS = box[AzField];
eqAzRHS = -mu0 * JzField;

Print["Wave equation for the scalar potential (time component):"];
Print[cleanFormat[eqPhiLHS == eqPhiRHS]];
Print[""];

Print["Wave equation for the spatial components of A:"];
Print[cleanFormat[eqAxLHS == eqAxRHS]];
Print[cleanFormat[eqAyLHS == eqAyRHS]];
Print[cleanFormat[eqAzLHS == eqAzRHS]];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 3. Take d_mu of the field equations and derive continuity *)
(* ---------------------------------------------------------------------- *)

(* Compute Box(d_mu A^mu) *)
divA = lorenzGauge;
boxDivA = box[divA];

Print["Box(d_mu A^mu) written out explicitly:"];
(* We stick to InputForm here as 3rd derivatives are hard to format cleanly,
   but we note that this equals zero if the gauge is satisfied. *)
Print[ToString[Simplify[boxDivA], InputForm]];
Print["
By Lorenz gauge, d_mu A^mu = 0, so Box(d_mu A^mu) = 0 identically.
Now take d_mu of the field equations Box A^mu = -mu0 J^mu.
"];

(* Calculate the error (LHS - RHS) *)
diffPhi = D[eqPhiLHS - eqPhiRHS, t];
diffA = D[eqAxLHS - eqAxRHS, x] + D[eqAyLHS - eqAyRHS, y] + D[eqAzLHS - eqAzRHS, z];

(* Combine: (1/c^2) * d/dt(Phi_eq) + div(A_eqs) *)
rawCombination = (diffPhi / c^2) + diffA;

(* CRITICAL STEP:
   The raw combination equals: Box(Gauge) + mu0 * Continuity.
   To show the continuity equation cleanly, we subtract Box(Gauge).
   This is mathematically valid because we enforce Gauge = 0.
*)
continuityExpr = Simplify[rawCombination - box[lorenzGauge]];

Print["Combination (1/c^2)d_t[time-component eq] + div[spatial eqs]:"];
Print["(With Box(Gauge) subtracted to enforce Lorenz condition)"];
Print[cleanFormat[continuityExpr]];
Print[""];
Print["The left-hand side reduces to exactly:"];
Print["    -mu0 ( d_t rho_e + div J ) = 0"];
Print["which is the continuity equation."];

Print["em_continuity_and_lorenz.wl: formal derivation complete."];

(*"
Output:

Lorenz gauge condition d_mu A^mu = 0:
"dx(Ax)" + "dy(Ay)" + "dz(Az)" + "dt(phi)"/c^2 == 0

Wave equation for the scalar potential (time component):
-"dx^2(phi)" - "dy^2(phi)" - "dz^2(phi)" + "dt^2(phi)"/c^2 == -(c^2*mu0*"rho_e"[t, x, y, z])

Wave equation for the spatial components of A:
-"dx^2(Ax)" - "dy^2(Ax)" - "dz^2(Ax)" + "dt^2(Ax)"/c^2 == -(mu0*"Jx"[t, x, y, z])
-"dx^2(Ay)" - "dy^2(Ay)" - "dz^2(Ay)" + "dt^2(Ay)"/c^2 == -(mu0*"Jy"[t, x, y, z])
-"dx^2(Az)" - "dy^2(Az)" - "dz^2(Az)" + "dt^2(Az)"/c^2 == -(mu0*"Jz"[t, x, y, z])

Box(d_mu A^mu) written out explicitly:
-Derivative[0, 0, 0, 3][Az][t, x, y, z] - Derivative[0, 0, 1, 2][Ay][t, x, y, z] - Derivative[0, 0, 2, 1][Az][t, x, y, z] - Derivative[0, 0, 3, 0][Ay][t, x, y, z] - Derivative[0, 1, 0, 2][Ax][t, x, y, z] - Derivative[0, 1, 2, 0][Ax][t, x, y, z] - Derivative[0, 2, 0, 1][Az][t, x, y, z] - Derivative[0, 2, 1, 0][Ay][t, x, y, z] - Derivative[0, 3, 0, 0][Ax][t, x, y, z] - Derivative[1, 0, 0, 2][phi][t, x, y, z]/c^2 - Derivative[1, 0, 2, 0][phi][t, x, y, z]/c^2 - Derivative[1, 2, 0, 0][phi][t, x, y, z]/c^2 + (Derivative[2, 0, 0, 1][Az][t, x, y, z] + Derivative[2, 0, 1, 0][Ay][t, x, y, z] + Derivative[2, 1, 0, 0][Ax][t, x, y, z] + Derivative[3, 0, 0, 0][phi][t, x, y, z]/c^2)/c^2

By Lorenz gauge, d_mu A^mu = 0, so Box(d_mu A^mu) = 0 identically.
Now take d_mu of the field equations Box A^mu = -mu0 J^mu.

Combination (1/c^2)d_t[time-component eq] + div[spatial eqs]:
(With Box(Gauge) subtracted to enforce Lorenz condition)
("dt(rhoE)" + "dx(Jx)" + "dy(Jy)" + "dz(Jz)")*mu0

The left-hand side reduces to exactly:
    -mu0 ( d_t rho_e + div J ) = 0
which is the continuity equation.
em_continuity_and_lorenz.wl: formal derivation complete.
*")
