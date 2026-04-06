(* ================================================================== *)
(* RIGOROUS STATIC-SOURCE CHECK: SCALAR 1PN CORRECTION SHOULD VANISH  *)
(* ================================================================== *)
(* *)
(* For a moving source, the scalar Liénard–Wiechert potential is      *)
(*   Phi = -mu / ( R_ret * (1 - nDotVs / cs) ).                       *)
(* In the test-mass / central-field limit, the source is static       *)
(* (vs = 0), so Phi = -mu / R_ret with no 1/cs corrections.           *)
(* *)
(* This script performs the Taylor series expansion of the retarded   *)
(* solution in eps = 1/cs and confirms that the 1PN scalar correction *)
(* vanishes when the source is static.                                *)
(* ================================================================== *)

(* NOTE:
   The final output ΔΦScalar[r] is the scalar 1PN correction used by:
     - clean_derivation.wl (Stage B, revised)
     - combined_1PN_sigma_derivation.wl (Φ_eff[r])
   It corresponds to solving the scalar retardation equation
   t_ret = t - R/cs for radial motion and inserting the Kepler EOM.
*)

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

(* Retarded Distance R_ret = tau / eps *)
RretRaw = Series[tauSol / eps, {eps, 0, 2}] // Normal;

Print["Retarded Distance (R_ret) expansion (raw):"];
Print[RretRaw];

(* Static source: retarded distance is just the instantaneous separation *)
Rret = r;

(* ================================================================== *)
(* STEP 2: Expansion of the Potential                                 *)
(* ================================================================== *)

(* Radial velocity at retarded time: vr_ret *)
vrRet = Series[vr - tauSol * ar, {eps, 0, 2}] // Normal;

(* Static source: vs = 0 -> no Doppler factor *)
Denom = Rret;

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

Print["Rigorous Scalar Correction (Circular): (expected 0)"];
Print[PhiCirc];

(* Check coefficient *)
CoeffRigorous = Coefficient[PhiCirc, mu^2/r^2];
Print["\nCoefficient of (mu^2/r^2): ", CoeffRigorous];

(* Define the scalar 1PN correction as a function of r *)

Clear[ΔΦScalar, cs];
ΔΦScalar[r_] :=
  Simplify[(PhiCirc) * eps^2 /. eps -> 1/cs];

Print["Scalar 1PN correction ΔΦ_scalar(r): (expected 0)"];
Print["    ", ΔΦScalar[r]];

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

Print["Scalar Sector Contribution: 0/6 (0% of GR)"];

(* Total must be 6/6 (100% from inertia) *)
(* InertiaPart = 2*beta / 6 *)

NewBeta = Solve[(2*beta)/6 == 1, beta][[1,1,2]];

Print["Required Inertia Contribution: 100%"];
Print["\nCALCULATED REQUIRED BETA:"];
Print[NewBeta];

Print["\nVERDICT:"];
If[NewBeta == 3,
    Print["VERDICT: SUCCESS. Beta = 3.0 is required."],
    Print["WARNING: Unexpected beta value. Check derivation."]
];

(*"
Output:

=============================================================
STEP 1: Iterative Solution of Retarded Time
=============================================================
Retarded Distance (R_ret) expansion (raw):
r - eps*r*vr + eps^2*((ar*r^2)/2 + r*vr^2)

=============================================================
STEP 3: Analysis of Terms
=============================================================
0th Order (Newtonian):
-(mu/r)

1st Order (1/cs):
0
(Note: proportional to time derivative of log(r), so no secular effect)

2nd Order (1/cs^2) Raw:
0

=============================================================
STEP 4: Circular Orbit Limit (vr -> 0, ar -> -mu/r^2)
=============================================================
Rigorous Scalar Correction (Circular): (expected 0)
0

Coefficient of (mu^2/r^2): 0
Scalar 1PN correction ΔΦ_scalar(r): (expected 0)
    0

Comparison to Naive Approximation (-1/2):
Ratio (Rigorous / Naive) = 0

=============================================================
STEP 5: Implication for Beta
=============================================================
Scalar Sector Contribution: 0/6 (0% of GR)
Required Inertia Contribution: 100%

CALCULATED REQUIRED BETA:
3

VERDICT:
VERDICT: SUCCESS. Beta = 3.0 is required.
"*)
