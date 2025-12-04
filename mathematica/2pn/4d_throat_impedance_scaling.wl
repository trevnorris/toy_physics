(* Script: 4d_throat_impedance_scaling.wl *)
(* Purpose: Determine the scaling of Mass M vs Radius a due to 3D->4D flow transition *)

ClearAll["Global`*"]

(* 1. Define Geometries *)
(* Brane (3D): Surface Area of sphere at radius r *)
AreaBrane[r_] := 4 * Pi * r^2;

(* Bulk (4D): The throat cross-section is a 3D ball of radius a *)
(* The "Area" flux passes through is Volume-like in 3D *)
AreaThroat[a_] := (4/3) * Pi * a^3;

(* 2. Define Velocity Profiles based on Dimensionality *)
(* v3D scales as 1/r^2 to conserve flux in 3D *)
vBrane[r_, Flux_] := Flux / (rho0 * AreaBrane[r]);

(* v4D in the throat: Assumed uniform flow through the 3D cross-section *)
vThroat[a_, Flux_] := Flux / (rho0 * AreaThroat[a]);

(* 3. Bernoulli Matching (Pressure Drop) *)
(* Total Energy/Enthalpy is conserved along streamline *)
(* P_inf + 0.5 rho v_inf^2 = P_throat + 0.5 rho v_throat^2 *)
(* We assume P_inf is fixed background pressure. *)
(* We assume P_throat is the vacuum pressure (or fixed low pressure). *)
(* Thus, DeltaP is constant. *)

EqBernoulli = DeltaP == (1/2) * rho0 * (vThroat[a, M]^2 - 0);

(* 4. Solve for Mass Flux M as a function of a *)
(* We assume DeltaP and rho0 are constants independent of geometry *)
sol = Solve[EqBernoulli, M];

MassScaling = Simplify[M /. sol[[2]]]; (* Take positive solution *)

Print["--- Scaling Result ---"];
Print["Mass Flux M scales as: "];
Print[MassScaling];

(* Check power of a *)
ExponentOfA = Exponent[MassScaling, a];
Print["Scaling Power of a: ", ExponentOfA];

(* 5. Compare to Volume Scaling *)
Print[""];
Print["Naive Volume Mass Scaling (M ~ a^3 * L): a^3 (assuming L const)"];
Print["Or M ~ a^4 (assuming L ~ a)"];

(* 6. Interpret 4D 'Resistance' *)
(* If flow is dominated by the 'squeeze' into 4D, does it change? *)
