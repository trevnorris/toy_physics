#!/usr/bin/env wolframscript

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

banner[title_] := Module[{line = StringRepeat["=", 88]},
  Print[""];
  Print[line];
  Print[title];
  Print[line];
];

checkZero[name_, expr_] := Module[{res = FullSimplify[Expand[expr]]},
  Print[name, " = ", res];
  If[res === 0,
    passCount++,
    failCount++;
    Print["FAIL: ", name];
    Throw[$Failed]
  ];
];

eps8Coeff[expr_] := FullSimplify[SeriesCoefficient[expr, {eps, 0, 8}]];

slotCoeff[expr_, upow_, vpow_] := FullSimplify[
  Coefficient[Coefficient[eps8Coeff[expr], U, upow], v, vpow]
];

banner["EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 4PN"];

lexact = Expand @ Normal @ Series[
   -1/eps^2 Sqrt[
     ((1 - U eps^2/2)/(1 + U eps^2/2))^2 -
     (1 + U eps^2/2)^4 v^2 eps^2
   ],
   {eps, 0, 8}
];

Print["L_exact/m = ", lexact];
checkZero["c^-8 coefficient of v^10 - 7/256", slotCoeff[lexact, 0, 10] - 7/256];
checkZero["c^-8 coefficient of U v^8 - 75/128", slotCoeff[lexact, 1, 8] - 75/128];
checkZero["c^-8 coefficient of U^2 v^6 - 59/16", slotCoeff[lexact, 2, 6] - 59/16];
checkZero["c^-8 coefficient of U^3 v^4 - 203/32", slotCoeff[lexact, 3, 4] - 203/32];
checkZero["c^-8 coefficient of U^4 v^2 - 31/32", slotCoeff[lexact, 4, 2] - 31/32];
checkZero["c^-8 coefficient of U^5 - 1/16", slotCoeff[lexact, 5, 0] - 1/16];

banner["CARRIED 3PN PACKAGING EXTENDED TO QUARTIC ORDER"];

denom = 1 - 4 U eps^2 + 8 U^2 eps^4 + d3 U^3 eps^6 + d4 U^4 eps^8;
lred = Expand @ Normal @ Series[
   -1/eps^2 (1 - U eps^2) Sqrt[1 - v^2 eps^2/denom],
   {eps, 0, 8}
];

lcandidateMin = Expand[
  lred
  - U^2 eps^2/2
  + U^3 eps^4/4
  - muRho3 U^4 eps^6/2
  + muRho4 U^5 eps^8/2
  + s24 U^2 v^4 eps^6
];

carriedVals = {d3 -> -45/4, muRho3 -> 1/4, s24 -> -1/16};
residualMin = Expand[lexact - (lcandidateMin /. carriedVals)];

Print["Residual with only d4 and mu_rho4 still open = ", residualMin];
checkZero["residual U^5 / c^8 slot - (1/16 - mu_rho4/2)", Coefficient[Coefficient[Coefficient[residualMin, eps, 8], U, 5], v, 0] - (1/16 - muRho4/2)];
checkZero["residual U^4 v^2 / c^8 slot + 205/32 - d4/2", Coefficient[Coefficient[Coefficient[residualMin, eps, 8], U, 4], v, 2] + 205/32 - d4/2];

banner["MINIMAL 4PN ONE-BODY REPAIR LEDGER"];

lcandidate = Expand[
  lred
  - U^2 eps^2/2
  + U^3 eps^4/4
  - muRho3 U^4 eps^6/2
  + muRho4 U^5 eps^8/2
  + s24 U^2 v^4 eps^6
  + s34 U^3 v^4 eps^8
  + s26 U^2 v^6 eps^8
];

residual = Expand[lexact - (lcandidate /. carriedVals)];
muSol = muRho4 /. First @ Solve[Coefficient[Coefficient[Coefficient[residual, eps, 8], U, 5], v, 0] == 0, muRho4];
d4Sol = d4 /. First @ Solve[Coefficient[Coefficient[Coefficient[residual, eps, 8], U, 4], v, 2] == 0, d4];
s34Sol = s34 /. First @ Solve[Coefficient[Coefficient[Coefficient[residual, eps, 8], U, 3], v, 4] == 0, s34];
s26Sol = s26 /. First @ Solve[Coefficient[Coefficient[Coefficient[residual, eps, 8], U, 2], v, 6] == 0, s26];

Print["mu_rho4 = ", muSol];
Print["d4 = ", d4Sol];
Print["s34 = ", s34Sol];
Print["s26 = ", s26Sol];

checkZero["residual after solving", residual /. {muRho4 -> muSol, d4 -> d4Sol, s34 -> s34Sol, s26 -> s26Sol}];

banner["FINAL 4PN ONE-BODY KICKOFF LEDGER"];
Print["mu_rho4 = ", muSol];
Print["d4 = ", d4Sol];
Print["s34 = ", s34Sol];
Print["s26 = ", s26Sol];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 4PN
========================================================================================
L_exact/m = -eps^(-2) + U - (eps^2*U^2)/2 + (eps^4*U^3)/4 - (eps^6*U^4)/8 + (eps^8*U^5)/16 + v^2/2 + (3*eps^2*U*v^2)/2 + 2*eps^4*U^2*v^2 + (13*eps^6*U^3*v^2)/8 + (31*eps^8*U^4*v^2)/32 + (eps^2*v^4)/8 + (7*eps^4*U*v^4)/8 + (47*eps^6*U^2*v^4)/16 + (203*eps^8*U^3*v^4)/32 + (eps^4*v^6)/16 + (11*eps^6*U*v^6)/16 + (59*eps^8*U^2*v^6)/16 + (5*eps^6*v^8)/128 + (75*eps^8*U*v^8)/128 + (7*eps^8*v^10)/256
c^-8 coefficient of v^10 - 7/256 = 0
c^-8 coefficient of U v^8 - 75/128 = 0
c^-8 coefficient of U^2 v^6 - 59/16 = 0
c^-8 coefficient of U^3 v^4 - 203/32 = 0
c^-8 coefficient of U^4 v^2 - 31/32 = 0
c^-8 coefficient of U^5 - 1/16 = 0

========================================================================================
CARRIED 3PN PACKAGING EXTENDED TO QUARTIC ORDER
========================================================================================
Residual with only d4 and mu_rho4 still open = (eps^8*U^5)/16 - (eps^8*muRho4*U^5)/2 - (205*eps^8*U^4*v^2)/32 + (d4*eps^8*U^4*v^2)/2 - (15*eps^8*U^3*v^4)/32 - (eps^8*U^2*v^6)/16
residual U^5 / c^8 slot - (1/16 - mu_rho4/2) = 0
residual U^4 v^2 / c^8 slot + 205/32 - d4/2 = 0

========================================================================================
MINIMAL 4PN ONE-BODY REPAIR LEDGER
========================================================================================
mu_rho4 = 1/8
d4 = 205/16
s34 = -15/32
s26 = -1/16
residual after solving = 0

========================================================================================
FINAL 4PN ONE-BODY KICKOFF LEDGER
========================================================================================
mu_rho4 = 1/8
d4 = 205/16
s34 = -15/32
s26 = -1/16
Passes: 9
Fails: 0
"*)
