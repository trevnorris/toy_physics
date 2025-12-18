(* ---------------------------------------------------------------------- *)
(* Script: Lensing_from_Flow_Fixed.wl *)
(* Purpose: Calculate Lensing/Shadow size using the Corrected Hybrid Model *)
(* (Scalar Potential + Vector Flow) *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Hybrid Lensing Verification"];
Print["======================================================="];

(* 1. Setup the Hybrid Flow (Scalar+Vector) *)
n = 5;                  (* Stiff Superfluid *)
rH = 1.0;               (* Horizon Radius *)
rhoH = 1.0;             (* Density at Horizon *)
csH = rhoH^((n-1)/2);   (* Speed of Sound at Horizon *)

(* FIX 1: Continuity uses r^2 for 3D flux conservation *)
C1Val = rhoH * csH * rH^2; 

(* FIX 2: Define the Gravitational Potential (Scalar Sector) *)
(* We set GM such that the potential is relevant. *)
(* Tuning parameter: For n=5, a simple choice is GM=0.5 to keep scales aligned *)
GM = 0.5; 
Phi[r_] := -GM / r;

(* Bernoulli Constant at Horizon: Enthalpy + Kinetic + Potential *)
(* h = n/(n-1) * rho^(n-1) *)
C2Val = (n/(n-1))*rhoH^(n-1) + csH^2/2 + Phi[rH]; 

(* Solve Velocity v[r] numerically with Hybrid Bernoulli Eq *)
(* Eq: h(rho) + v^2/2 + Phi(r) = Const *)
(* Substitute rho = C1 / (v * r^2) *)
solV[radius_] := Module[{vSol},
   (* We look for the subsonic branch (small v) in the far field *)
   (* Use Quiet to suppress precision warnings during root finding *)
   vSol = v /. Quiet[FindRoot[
      ((n/(n-1)) * (C1Val/(v * radius^2))^(n-1) + v^2/2 + Phi[radius]) - C2Val == 0,
      {v, 0.001} (* Initial guess: very low velocity far away *)
   ]];
   vSol
];

(* 2. Define Refractive Index N(r) *)
(* FIX 3: Density uses r^2 consistency *)
rho[r_] := C1Val / (solV[r] * r^2);

(* Refractive Index N(r) ~ 1/c_s(r) *)
(* Normalized so N -> 1 at Infinity (r = 100.0 is approx infinity) *)
RefractiveIndex[r_] := (rho[100.0]^((n-1)/2)) / (rho[r]^((n-1)/2));

(* 3. Lensing Integral *)
(* Deflection Theta = 2 * Integral [ (b/r) * (dN/dr) / sqrt((N*r)^2 - b^2) ] dr *)

DeflectionAngle[b_] := Module[{Integrand, r0, rMax},
   (* Derivative dN/dr calculated numerically for stability *)
   dNdr[rad_] := (RefractiveIndex[rad + 0.001] - RefractiveIndex[rad]) / 0.001;
   
   (* Integrand Definition *)
   Integrand[r_] := (b / r) * dNdr[r] / Sqrt[(RefractiveIndex[r]*r)^2 - b^2];

   (* Find Turning Point r0 where denominator is zero: N(r0)*r0 = b *)
   (* Approx r0 ~ b for weak lensing *)
   r0 = b; 
   (* Refine r0 slightly if needed, but for b >> rH, r0=b is safe start *)
   
   (* Integrate from turning point to Infinity (approx 1000.0) *)
   rMax = 1000.0;
   
   Quiet[NIntegrate[Integrand[r], {r, r0, rMax}]] * 2
];

(* 4. Run Verification *)
Print["\n--- Computing Lensing Angles ---"];

(* Test two impact parameters in the weak field *)
b1 = 10.0;
b2 = 20.0;

theta1 = Abs[DeflectionAngle[b1]];
theta2 = Abs[DeflectionAngle[b2]];

Print["Deflection at b=10: ", theta1];
Print["Deflection at b=20: ", theta2];

ratio = theta1 / theta2;
Print["\nRatio (Theta1/Theta2): ", ratio];

Print["\n--- Verdict ---"];
Print["Expected Ratios:"];
Print["  GR (1/b scaling):      2.0"];
Print["  Short Range (1/b^4):   16.0"];

If[Abs[ratio - 2.0] < 0.2,
    Print["SUCCESS: Matches GR! The Scalar Potential dominates the far-field optics."],
    Print["FAILURE: Still mismatches GR. Check potential coupling."]
];

(*"
Output:

=======================================================
   Paper VI: Hybrid Lensing Verification
=======================================================

--- Computing Lensing Angles ---
Deflection at b=10: 0.009965822770890285
Deflection at b=20: 0.01946976392826504

Ratio (Theta1/Theta2): 0.5118615103710887

--- Verdict ---
Expected Ratios:
  GR (1/b scaling):      2.0
  Short Range (1/b^4):   16.0
FAILURE: Still mismatches GR. Check potential coupling.
"*)
