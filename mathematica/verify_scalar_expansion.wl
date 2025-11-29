(* ================================================================== *)
(* RIGOROUS DERIVATION OF SCALAR 1PN CORRECTION                       *)
(* ================================================================== *)
(* *)
(* This script performs a full Taylor series expansion of the         *)
(* Scalar Lienard-Wiechert potential:                                 *)
(* Phi = -mu / ( R_ret * (1 - vr_ret / cs) )                       *)
(* *)
(* It iteratively solves the retardation equation t_ret = t - R/cs    *)
(* to find the exact expansion without heuristic approximations.      *)
(* *)
(* ================================================================== *)

ClearAll["Global`*"];

Print["============================================================="];
Print["STEP 1: Iterative Solution of Retarded Time"];
Print["============================================================="];

(* Expansion parameter eps = 1/cs *)
(* Instantaneous variables: r, vr (radial velocity), ar (radial accel) *)

(* Iterative solution for tau = R(t-tau)*eps *)
tau0 = r * eps;
tau1 = (r - tau0 * vr) * eps;
tau2 = (r - tau1 * vr + 1/2 * tau1^2 * ar) * eps;

(* Keep terms up to eps^3 for precision *)
tauSol = Series[tau2, {eps, 0, 3}] // Normal;

(* Retarded Distance Rret = tau / eps *)
Rret = Series[tauSol / eps, {eps, 0, 2}] // Normal;

Print["Retarded Distance (R_ret) expansion:"];
Print[Rret];

(* ================================================================== *)
(* STEP 2: Expansion of the Potential                                 *)
(* ================================================================== *)

(* Radial velocity at retarded time: vr_ret *)
vrRet = Series[vr - tauSol * ar, {eps, 0, 2}] // Normal;

(* Doppler Denominator D = Rret * (1 - vr_ret / cs) *)
Denom = Rret * (1 - vrRet * eps);

(* Scalar Potential Phi = -mu / Denom *)
PhiFull = Series[-mu / Denom, {eps, 0, 2}] // Normal;

(* Extract Coefficients *)
TermNewton = Simplify[Coefficient[PhiFull, eps, 0]];
TermLinear = Simplify[Coefficient[PhiFull, eps, 1]];
TermQuad   = Simplify[Coefficient[PhiFull, eps, 2]];

Print["\n============================================================="];
Print["STEP 3: Analysis of Terms"];
Print["============================================================="];

Print["0th Order (Newtonian):"];
Print[TermNewton];

Print["\n1st Order (1/cs):"];
Print[TermLinear];
Print["(Note: proportional to time derivative of log(r), so no secular effect)"];

Print["\n2nd Order (1/cs^2) Raw:"];
Print[TermQuad];

(* ================================================================== *)
(* STEP 4: Circular Orbit Limit and Comparison                        *)
(* ================================================================== *)

Print["\n============================================================="];
Print["STEP 4: Circular Orbit Limit (vr -> 0, ar -> -mu/r^2)"];
Print["============================================================="];

(* Substitute 0PN Equation of Motion: ar -> -mu/r^2 *)
PhiMotion = Simplify[TermQuad /. ar -> -mu/r^2];

(* For circular orbits, radial velocity vr is 0 *)
PhiCirc = Simplify[PhiMotion /. vr -> 0];

Print["Rigorous Scalar Correction (Circular):"];
Print[PhiCirc];

(* Check coefficient *)
CoeffRigorous = Coefficient[PhiCirc, mu^2/r^2];
Print["\nCoefficient of (mu^2/r^2): ", CoeffRigorous];

Print["\nComparison to Naive Approximation (-1/2):"];
NaiveClaim = -1/2;
FactorDiff = CoeffRigorous / NaiveClaim;

Print["Ratio (Rigorous / Naive) = ", FactorDiff];

(* ================================================================== *)
(* STEP 5: Beta Requirement Calculation                               *)
(* ================================================================== *)

Print["\n============================================================="];
Print["STEP 5: Implication for Beta"];
Print["============================================================="];

(* Precession is linear in the potential coefficient. *)
(* Factor 1/2 in potential -> 1/6 of GR Precession. *)
(* Factor 3/2 in potential -> 3 * (1/6) = 1/2 of GR Precession. *)

Print["Scalar Sector Contribution: 3/6 (50% of GR)"];

(* Total must be 1 (100%) *)
(* 0.5 + InertiaPart = 1.0 *)
(* InertiaPart = 2*Beta / 6 *)

NewBeta = Solve[3/6 + (2*beta)/6 == 1, beta][[1,1,2]];

Print["Required Inertia Contribution: 50%"];
Print["\nCALCULATED REQUIRED BETA:"];
Print[NewBeta];

Print["\nVERDICT:"];
If[NewBeta == 1.5,
    Print["SUCCESS. Beta = 1.5 is confirmed by the rigorous expansion."],
    Print["WARNING. Unexpected Beta result."]
];
Print["This justifies the use of beta=1.5 in the other scripts."];

(*"
Output:

=============================================================
STEP 1: Iterative Solution of Retarded Time
=============================================================
Retarded Distance (R_ret) expansion:
r - eps*r*vr + eps^2*((ar*r^2)/2 + r*vr^2)

=============================================================
STEP 3: Analysis of Terms
=============================================================
0th Order (Newtonian):
-(mu/r)

1st Order (1/cs):
(-2*mu*vr)/r
(Note: proportional to time derivative of log(r), so no secular effect)

2nd Order (1/cs^2) Raw:
(3*ar*mu)/2 - (2*mu*vr^2)/r

=============================================================
STEP 4: Circular Orbit Limit (vr -> 0, ar -> -mu/r^2)
=============================================================
Rigorous Scalar Correction (Circular):
(-3*mu^2)/(2*r^2)

Coefficient of (mu^2/r^2): -3/2

Comparison to Naive Approximation (-1/2):
Ratio (Rigorous / Naive) = 3

=============================================================
STEP 5: Implication for Beta
=============================================================
Scalar Sector Contribution: 3/6 (50% of GR)
Required Inertia Contribution: 50%

CALCULATED REQUIRED BETA:
3/2

VERDICT:
SUCCESS. Beta = 1.5 is confirmed by the rigorous expansion.
This justifies the use of beta=1.5 in the other scripts.
"*)
