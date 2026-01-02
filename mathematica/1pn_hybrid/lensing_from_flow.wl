(* ---------------------------------------------------------------------- *)
(* Script: Lensing_from_Flow_Fixed.wl *)
(* Purpose: Calculate Lensing/Shadow size using the Corrected Hybrid Model *)
(* (Scalar Potential + Vector Flow) *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* --- Common polytrope normalization (Paper VI convention) --- *)
n = 5;

(* Choose units: background density rho0 and background wave speed c0 *)
rho0 = 1.0;
c0 = 1.0;

(* EOS: P = K rho^n with K chosen so cs(rho0)=c0 *)
K = c0^2/(n * rho0^(n - 1));

Pressure[rho_] := K * rho^n;
cs2[rho_] := D[Pressure[rho], rho];          (* = n K rho^(n-1) *)
cs[rho_] := Sqrt[cs2[rho]];

(* Enthalpy: h = ∫ dP/rho = (n K/(n-1)) rho^(n-1) (up to additive constant) *)
Enthalpy[rho_] := (n*K/(n-1)) * rho^(n - 1);

Print["======================================================="];
Print["   Paper VI: Hybrid Lensing Verification"];
Print["======================================================="];

(* 1. Setup the Hybrid Flow (Scalar+Vector) *)
rH = 1.0;               (* Horizon Radius *)
rhoH = 1.0;             (* Density at Horizon *)
csH = cs[rhoH];         (* Speed of Sound at Horizon *)
vH = csH;               (* Sonic at the horizon *)

(* FIX 1: Continuity uses r^2 for 3D flux conservation *)
C1Val = rhoH * vH * rH^2; 

(* FIX 2: Define the Gravitational Potential (Scalar Sector) *)
(* We set GM such that the potential is relevant. *)
(* Tuning parameter: For n=5, a simple choice is GM=0.5 to keep scales aligned *)
GM = 0.5; 
Phi[r_] := -GM / r;

(* Bernoulli Constant at Horizon: Enthalpy + Kinetic + Potential *)
C2Val = Enthalpy[rhoH] + vH^2/2 + Phi[rH]; 

(* --- build a stable subsonic branch v(r) by continuation --- *)
rMax = 1000.0;
dr1 = 0.05;     (* fine near horizon *)
dr2 = 0.5;      (* coarser far out *)
rGrid = Join[Range[rH, 20.0, dr1], Range[20.0 + dr2, rMax, dr2]];

vList = {};
lastV = vH;   (* start sonic at horizon *)

Do[
  sol = Quiet@FindRoot[
    Enthalpy[C1Val/(v*r^2)] + v^2/2 + Phi[r] - C2Val == 0,
    {v, lastV, 10^-6, 2.0}
  ];
  lastV = v /. sol;
  AppendTo[vList, lastV];
, {r, rGrid}];

vInterp = Interpolation[Transpose[{rGrid, vList}], InterpolationOrder -> 3];

rhoOfR[r_] := C1Val/(vInterp[r]*r^2);
csOfR[r_] := cs[rhoOfR[r]];

(* Paper VI / optics convention: N = c/c_s; here c0 is background c *)
NOfR[r_] := c0/csOfR[r];

(* normalize N -> 1 at large r *)
Ninf = NOfR[rMax];
Nnorm[r_] := NOfR[r]/Ninf;

NInterp = Interpolation[Table[{r, Nnorm[r]}, {r, rGrid}], InterpolationOrder -> 3];
dNdr[r_] := NInterp'[r];

(* 3. Lensing Integral *)
(* Deflection Theta = 2 * Integral [ (b/r) * (dN/dr) / sqrt((N*r)^2 - b^2) ] dr *)

DeflectionAngle[b_] := Module[{Integrand, r0, rMax},
   (* Integrand Definition *)
   Integrand[r_] := (b / r) * dNdr[r] / Sqrt[(NInterp[r]*r)^2 - b^2];

   (* Find Turning Point r0 where denominator is zero: N(r0)*r0 = b *)
   r0 = r /. Quiet@FindRoot[NInterp[r]*r - b == 0, {r, b, rH + 10^-3, b}];
   
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
