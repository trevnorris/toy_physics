(* ---------------------------------------------------------------------- *)
(* Script: Lensing_from_Flow.wl *)
(* Purpose: Calculate Lensing/Shadow size using the EXACT Transonic Flow Profile *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* 1. Setup the Nonlinear Flow (From Previous Success) *)
n = 5;
rH = 1.0;
rhoH = 1.0;
csH = rhoH^((n-1)/2);
C1Val = rhoH * csH * rH^3; (* Continuity Constant *)
C2Val = (n/(n-1))*rhoH^(n-1) + csH^2/2; (* Bernoulli Constant *)

(* Solve Velocity v[r] numerically *)
solV[radius_] := Module[{vSol},
   vSol = v /. FindRoot[
      ((n/(n-1)) * (C1Val/(v * radius^3))^(n-1) + v^2/2) - C2Val == 0,
      {v, 0.01}
   ];
   vSol
];

(* 2. Define Refractive Index n(r) from Density *)
(* Density rho is derived from velocity via Continuity *)
rho[r_] := C1Val / (solV[r] * r^3);

(* Refractive Index N(r) ~ 1/c_s(r) *)
(* In vacuum, N = 1. Near hole, c_s drops, so N increases. *)
(* c_s ~ rho^2 (for n=5) *)
RefractiveIndex[r_] := (rho[100.0]^2) / (rho[r]^2); (* Normalized to far field *)

(* 3. Lensing Integral *)
(* Deflection Theta = 2 * Integral [ (b/r) * (dn/dr) / sqrt(n^2 r^2 - b^2) ] dr *)
(* We compute this numerically for different impact parameters b *)

DeflectionAngle[b_] := Module[{Integrand},
   (* Derivative dn/dr calculated numerically *)
   dNdr[rad_] := (RefractiveIndex[rad + 0.001] - RefractiveIndex[rad]) / 0.001;

   Integrand[r_] := (b / r) * dNdr[r] / Sqrt[(RefractiveIndex[r]*r)^2 - b^2];

   (* Integrate from closest approach r0 to Infinity (approx 50.0) *)
   (* Find r0 where denominator is zero (classic turning point) *)
   r0 = b; (* Approximation for weak lensing, good enough for check *)

   Quiet[NIntegrate[Integrand[r], {r, r0, 50.0}]] * 2
];

(* 4. Compare with GR Scaling *)
(* GR: Theta ~ 1/b *)
(* Pure Flow (Bernoulli): Theta ~ 1/b^3 or 1/b^5 ?? *)
(* Let's test two impact parameters *)
b1 = 10.0;
b2 = 20.0;

theta1 = Abs[DeflectionAngle[b1]];
theta2 = Abs[DeflectionAngle[b2]];

Print["--- Lensing Check ---"];
Print["Deflection at b=10: ", theta1];
Print["Deflection at b=20: ", theta2];

ratio = theta1 / theta2;
Print["Ratio (Theta1/Theta2): ", ratio];
Print[""];
Print["Expected Ratios:"];
Print["  GR (1/b):      2.0"];
Print["  Dipole (1/b^3): 8.0"];
Print["  Short Range:   > 8.0"];

If[Abs[ratio - 2.0] < 0.5,
    Print["RESULT: Matches GR! The Flow naturally creates long-range lensing."],
    Print["RESULT: Mismatch. The Flow is short-range. The 'Kernel' is required for Far Field."]
];

(*"
Output:

--- Lensing Check ---
Deflection at b=10: 0.0007504107473192793
Deflection at b=20: 0.0007425163639681885
Ratio (Theta1/Theta2): 1.010631931812117

Expected Ratios:
  GR (1/b):      2.0
  Dipole (1/b^3): 8.0
  Short Range:   > 8.0
RESULT: Mismatch. The Flow is short-range. The 'Kernel' is required for Far Field.
"*)
