(* ---------------------------------------------------------------------- *)
(* Script: lensing_from_flow.wl                                            *)
(* Purpose: Evaluate a representative n=5 optical profile with the         *)
(* correct 1/r weak-field tail and a finite-size throat correction.        *)
(*                                                                        *)
(* Audit note: the hybrid paper's strong-field section is currently        *)
(* phenomenological. This script therefore uses a weak-field-matched       *)
(* optical profile, not a first-principles throat solution.                *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

lambdaStar = N[(Sqrt[2] * Pi)/BesselJZero[0, 1], 16];
sigmaFS = N[1/(1 + lambdaStar^2), 16];
mu = 1.0;
a = 1.0;

NProfile[r_] := 1 + (2 * mu)/r + sigmaFS/(r^2 + a^2);

photonSphereEquation[r_] := D[NProfile[x]^4 * x^2, x] /. x -> r;
rph = r /. FindRoot[photonSphereEquation[r] == 0, {r, 2.0}];
bph = NProfile[rph]^2 * rph;

TurningPoint[b_] := r /. FindRoot[NProfile[r]^2 * r - b == 0, {r, 1.1, b}];

DeflectionAngle[b_] := Module[{r0, integrand, rMax = 500.0},
  r0 = TurningPoint[b];
  integrand[r_] := b/(r * Sqrt[NProfile[r]^4 * r^2 - b^2]);
  2 * NIntegrate[integrand[r], {r, r0 + 10^-5, rMax}] - Pi
];

b1 = 10.0;
b2 = 15.0;
theta1 = N[DeflectionAngle[b1], 16];
theta2 = N[DeflectionAngle[b2], 16];

Print["--- Representative Lensing Profile ---"];
Print["lambda_* = ", lambdaStar];
Print["sigma_fs = ", sigmaFS];
Print["r_ph = ", N[rph, 16]];
Print["b_ph = ", N[bph, 16]];
Print["Deflection[b=10] = ", theta1];
Print["Deflection[b=15] = ", theta2];
Print["Deflection ratio = ", N[theta1/theta2, 16]];

(*"
Output:

--- Representative Lensing Profile ---
lambda_* = 1.84748657712012805107642087596939904276`16.
sigma_fs = 0.2265926068524372168300616723354750574`15.81056168507603
r_ph = 2.19764991619614
b_ph = 8.347412193956503
Deflection[b=10] = 1.9172638577403216
Deflection[b=15] = 0.7681424872102203
Deflection ratio = 2.4959742361128074
"*)
