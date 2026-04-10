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

checkEqual[name_, lhs_, rhs_] := checkZero[name, lhs - rhs];

eps6Coeff[expr_] := FullSimplify[SeriesCoefficient[expr, {eps, 0, 6}]];

slotCoeff[expr_, upow_, vpow_] := FullSimplify[
  Coefficient[Coefficient[eps6Coeff[expr], U, upow], v, vpow]
];

banner["EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 3PN"];

lexact = Expand @ Normal @ Series[
   -1/eps^2 Sqrt[
     ((1 - U eps^2/2)/(1 + U eps^2/2))^2 -
     (1 + U eps^2/2)^4 v^2 eps^2
   ],
   {eps, 0, 6}
];

Print["L_exact/m = ", lexact];
checkEqual["c^-6 coefficient of v^8", slotCoeff[lexact, 0, 8], 5/128];
checkEqual["c^-6 coefficient of U v^6", slotCoeff[lexact, 1, 6], 11/16];
checkEqual["c^-6 coefficient of U^2 v^4", slotCoeff[lexact, 2, 4], 47/16];
checkEqual["c^-6 coefficient of U^3 v^2", slotCoeff[lexact, 3, 2], 13/8];
checkEqual["c^-6 coefficient of U^4", slotCoeff[lexact, 4, 0], -1/8];

banner["CARRIED DENOMINATOR-STYLE SELF SECTOR EXTENDED TO CUBIC ORDER"];

denom = 1 - 4 U eps^2 + 8 U^2 eps^4 + d3 U^3 eps^6;
lred = Expand @ Normal @ Series[
   -1/eps^2 (1 - U eps^2) Sqrt[1 - v^2 eps^2/denom],
   {eps, 0, 6}
];

Print["L_red/m = ", lred];
checkEqual["red c^-6 coefficient of v^8", slotCoeff[lred, 0, 8], 5/128];
checkEqual["red c^-6 coefficient of U v^6", slotCoeff[lred, 1, 6], 11/16];
checkEqual["red c^-6 coefficient of U^2 v^4", slotCoeff[lred, 2, 4], 3];
checkEqual["red c^-6 coefficient of U^3 v^2", slotCoeff[lred, 3, 2], -d3/2 - 4];

banner["MINIMAL 3PN ONE-BODY REPAIR LEDGER"];

lcandidate = Expand[
  lred
  - U^2 eps^2/2
  + U^3 eps^4/4
  - muRho3 U^4 eps^6/2
  + s24 U^2 v^4 eps^6
];

residual = Expand[lexact - lcandidate];
Print["Exact target minus candidate = ", residual];

muSol = muRho3 /. First @ Solve[
  Coefficient[Coefficient[Coefficient[residual, eps, 6], U, 4], v, 0] == 0,
  muRho3
];
d3Sol = d3 /. First @ Solve[
  Coefficient[Coefficient[Coefficient[residual, eps, 6], U, 3], v, 2] == 0,
  d3
];
s24Sol = s24 /. First @ Solve[
  Coefficient[Coefficient[Coefficient[residual, eps, 6], U, 2], v, 4] == 0,
  s24
];

Print["mu_rho3 = ", muSol];
Print["d3 = ", d3Sol];
Print["s24 = ", s24Sol];

checkZero[
  "residual after solving",
  residual /. {muRho3 -> muSol, d3 -> d3Sol, s24 -> s24Sol}
];

banner["UNEXTENDED 2PN INVARIANT DENOMINATOR PREDICTION AT CUBIC ORDER"];

u = Symbol["u"];
g1 = 57/64;
g2 = 298821/131072;
mu = 32768/3249;
gSeries = 1 + g1 u + g2 u^2;
dCarry = Expand[(1 - 4 u) (1 + mu (gSeries - 1)^2)];
dCarrySeries = Expand @ Normal @ Series[dCarry, {u, 0, 3}];

Print["D_carry(u) = ", dCarrySeries];
checkEqual["carried cubic coefficient d3_carry", Coefficient[dCarrySeries, u, 3], 21783/2432];
checkEqual["target cubic coefficient d3_target", d3Sol, -45/4];

banner["MINIMAL CUBIC GEOMETRY-INVARIANT CORRECTION"];

dRepaired = Expand[(1 - 4 u) (1 + mu (gSeries - 1)^2 + nuRepair (g1 u)^3)];
dRepairedSeries = Expand @ Normal @ Series[dRepaired, {u, 0, 3}];
Print["D_repaired(u) = ", dRepairedSeries];

nuRepairSol = nuRepair /. First @ Solve[
  Coefficient[dRepairedSeries, u, 3] == d3Sol,
  nuRepair
];

Print["nu = ", nuRepairSol];
checkEqual[
  "repaired cubic coefficient",
  Coefficient[dRepairedSeries, u, 3] /. nuRepair -> nuRepairSol,
  -45/4
];

banner["FINAL 3PN ONE-BODY KICKOFF LEDGER"];
Print["mu_rho3 = ", muSol];
Print["d3 = ", d3Sol];
Print["s24 = ", s24Sol];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 3PN
========================================================================================
L_exact/m = -eps^(-2) + U - (eps^2*U^2)/2 + (eps^4*U^3)/4 - (eps^6*U^4)/8 + v^2/2 + (3*eps^2*U*v^2)/2 + 2*eps^4*U^2*v^2 + (13*eps^6*U^3*v^2)/8 + (eps^2*v^4)/8 + (7*eps^4*U*v^4)/8 + (47*eps^6*U^2*v^4)/16 + (eps^4*v^6)/16 + (11*eps^6*U*v^6)/16 + (5*eps^6*v^8)/128
c^-6 coefficient of v^8 = 0
c^-6 coefficient of U v^6 = 0
c^-6 coefficient of U^2 v^4 = 0
c^-6 coefficient of U^3 v^2 = 0
c^-6 coefficient of U^4 = 0

========================================================================================
CARRIED DENOMINATOR-STYLE SELF SECTOR EXTENDED TO CUBIC ORDER
========================================================================================
L_red/m = -eps^(-2) + U + v^2/2 + (3*eps^2*U*v^2)/2 + 2*eps^4*U^2*v^2 - 4*eps^6*U^3*v^2 - (d3*eps^6*U^3*v^2)/2 + (eps^2*v^4)/8 + (7*eps^4*U*v^4)/8 + 3*eps^6*U^2*v^4 + (eps^4*v^6)/16 + (11*eps^6*U*v^6)/16 + (5*eps^6*v^8)/128
red c^-6 coefficient of v^8 = 0
red c^-6 coefficient of U v^6 = 0
red c^-6 coefficient of U^2 v^4 = 0
red c^-6 coefficient of U^3 v^2 = 0

========================================================================================
MINIMAL 3PN ONE-BODY REPAIR LEDGER
========================================================================================
Exact target minus candidate = -1/8*(eps^6*U^4) + (eps^6*muRho3*U^4)/2 + (45*eps^6*U^3*v^2)/8 + (d3*eps^6*U^3*v^2)/2 - (eps^6*U^2*v^4)/16 - eps^6*s24*U^2*v^4
mu_rho3 = 1/4
d3 = -45/4
s24 = -1/16
residual after solving = 0

========================================================================================
UNEXTENDED 2PN INVARIANT DENOMINATOR PREDICTION AT CUBIC ORDER
========================================================================================
D_carry(u) = 1 - 4*u + 8*u^2 + (21783*u^3)/2432
carried cubic coefficient d3_carry = 0
target cubic coefficient d3_target = 0

========================================================================================
MINIMAL CUBIC GEOMETRY-INVARIANT CORRECTION
========================================================================================
D_repaired(u) = 1 - 4*u + 8*u^2 + (21783*u^3)/2432 + (185193*nuRepair*u^3)/262144
nu = -33548288/1172889
repaired cubic coefficient = 0

========================================================================================
FINAL 3PN ONE-BODY KICKOFF LEDGER
========================================================================================
mu_rho3 = 1/4
d3 = -45/4
s24 = -1/16
Passes: 13
Fails: 0
"*)
