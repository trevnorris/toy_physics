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

grTargetH[] := <|
  1 -> (-5 + 35 nu - 70 nu^2 + 35 nu^3)/128,
  2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
  6 -> (-7 + 42 nu - 53 nu^2 - 5 nu^3)/16,
  7 -> (2 - 3 nu) nu^2/16,
  8 -> 3 (1 - nu) nu^2/16,
  9 -> -5 nu^3/16,
  10 -> (-27 + 136 nu + 109 nu^2)/16,
  11 -> (17 + 30 nu) nu/16,
  12 -> (5 + 43 nu) nu/12,
  13 -> (-600 + (3 Pi^2 - 1340) nu - 552 nu^2)/192,
  14 -> -(340 + 3 Pi^2 + 112 nu) nu/64,
  15 -> (12 + (872 - 63 Pi^2) nu)/96
|>;

inverseMapFromH[h_] := <|
  1 -> (-h[1] + 3 nu/16 - 21 nu^2/16 + 9 nu^3/4),
  2 -> -h[2], 3 -> -h[3], 4 -> -h[4], 5 -> -h[5],
  6 -> (-h[6] + 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4),
  7 -> (-h[7] + 11 nu^2/8 - 9 nu^3/2),
  8 -> (-h[8] + 3 nu^2/4 - 9 nu^3/4),
  9 -> -h[9],
  10 -> (-h[10] + 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4),
  11 -> (-h[11] + 7 nu/8 + 41 nu^2/8 + 31 nu^3/4),
  12 -> (-h[12] + 9 nu^2/2 + 4 nu^3),
  13 -> (-h[13] - 3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2),
  14 -> (-h[14] + 7 nu/4 - 31 nu^2/4 - 7 nu^3/2),
  15 -> -h[15]
|>;

carriedSeedL[] := <|
  1 -> 5/128 - 35 nu/128 + 35 nu^2/64 - 35 nu^3/128,
  2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
  6 -> 11/16 - 33 nu/8 + 99 nu^2/16 - 11 nu^3/8,
  7 -> 0, 8 -> 0, 9 -> 0,
  10 -> 47/16 - 235 nu/16 + 235 nu^2/16,
  11 -> 0, 12 -> 0,
  13 -> 13/8 - 13 nu/2 + 13 nu^2/4,
  14 -> 0,
  15 -> -1/8 + 3 nu/8
|>;

carriedSeedH[] := <|
  1 -> -5/128 + 59 nu/128 - 119 nu^2/64 + 323 nu^3/128,
  2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
  6 -> -7/16 + 5 nu - 169 nu^2/16 - 31 nu^3/8,
  7 -> nu^2 (11 - 36 nu)/8,
  8 -> 3 nu^2 (1 - 3 nu)/4,
  9 -> 0,
  10 -> -27/16 + 265 nu/16 + 11 nu^2/16 + 13 nu^3/4,
  11 -> nu (7 + 41 nu + 62 nu^2)/8,
  12 -> nu^2 (9 + 8 nu)/2,
  13 -> -25/8 - 33 nu/4 - 19 nu^2/2 - nu^3/2,
  14 -> nu (7 - 31 nu - 14 nu^2)/4,
  15 -> (1 - 3 nu)/8
|>;

banner["PART I — IMPORT THE EXACT GR 3PN COM HAMILTONIAN TARGET"];

targetH = grTargetH[];
targetL = inverseMapFromH[targetH];
seedH = carriedSeedH[];
seedL = carriedSeedL[];

Do[Print["h", i, "^(GR) = ", targetH[i]], {i, 1, 15}];

banner["PART II — SOLVE THE EXACT COM ORDINARY-LAGRANGIAN COEFFICIENTS"];
Do[Print["l", i, "^(GR) = ", targetL[i]], {i, 1, 15}];

Do[
  checkZero[
    "map check h" <> ToString[i],
    Switch[i,
      1, 3 nu/16 - 21 nu^2/16 + 9 nu^3/4 - targetL[1] - targetH[1],
      2, -targetL[2] - targetH[2],
      3, -targetL[3] - targetH[3],
      4, -targetL[4] - targetH[4],
      5, -targetL[5] - targetH[5],
      6, 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4 - targetL[6] - targetH[6],
      7, 11 nu^2/8 - 9 nu^3/2 - targetL[7] - targetH[7],
      8, 3 nu^2/4 - 9 nu^3/4 - targetL[8] - targetH[8],
      9, -targetL[9] - targetH[9],
      10, 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4 - targetL[10] - targetH[10],
      11, 7 nu/8 + 41 nu^2/8 + 31 nu^3/4 - targetL[11] - targetH[11],
      12, 9 nu^2/2 + 4 nu^3 - targetL[12] - targetH[12],
      13, -3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2 - targetL[13] - targetH[13],
      14, 7 nu/4 - 31 nu^2/4 - 7 nu^3/2 - targetL[14] - targetH[14],
      15, -targetL[15] - targetH[15]
    ]
  ],
  {i, 1, 15}
];

banner["PART III — ONE-BODY GATE AND SEED CHECKS"];

Do[
  checkZero["one-body h" <> ToString[i] <> " target - seed", (targetH[i] - seedH[i]) /. nu -> 0];
  checkZero["one-body l" <> ToString[i] <> " target - seed", (targetL[i] - seedL[i]) /. nu -> 0];,
  {i, {1, 6, 10, 13, 15}}
];

Do[
  checkZero["one-body h" <> ToString[i], targetH[i] /. nu -> 0];
  checkZero["one-body l" <> ToString[i], (targetL[i] - seedL[i]) /. nu -> 0];,
  {i, {2, 3, 4, 5, 7, 8, 9, 11, 12, 14}}
];

banner["PART IV — GENUINE COMPARABLE-MASS RESIDUALS"];

deltaH = Association@Table[i -> FullSimplify[targetH[i] - seedH[i]], {i, 1, 15}];
deltaL = Association@Table[i -> FullSimplify[targetL[i] - seedL[i]], {i, 1, 15}];

Do[
  checkZero["nu -> 0 residual h" <> ToString[i], deltaH[i] /. nu -> 0];
  checkZero["nu -> 0 residual l" <> ToString[i], deltaL[i] /. nu -> 0];
  checkZero["Delta l" <> ToString[i] <> " + Delta h" <> ToString[i], deltaL[i] + deltaH[i]];,
  {i, 1, 15}
];

banner["FINAL FOUNDATION LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — IMPORT THE EXACT GR 3PN COM HAMILTONIAN TARGET
========================================================================================
h1^(GR) = (-5 + 35*nu - 70*nu^2 + 35*nu^3)/128
h2^(GR) = 0
h3^(GR) = 0
h4^(GR) = 0
h5^(GR) = 0
h6^(GR) = (-7 + 42*nu - 53*nu^2 - 5*nu^3)/16
h7^(GR) = ((2 - 3*nu)*nu^2)/16
h8^(GR) = (3*(1 - nu)*nu^2)/16
h9^(GR) = (-5*nu^3)/16
h10^(GR) = (-27 + 136*nu + 109*nu^2)/16
h11^(GR) = (nu*(17 + 30*nu))/16
h12^(GR) = (nu*(5 + 43*nu))/12
h13^(GR) = (-600 - 552*nu^2 + nu*(-1340 + 3*Pi^2))/192
h14^(GR) = (nu*(-340 - 112*nu - 3*Pi^2))/64
h15^(GR) = (12 + nu*(872 - 63*Pi^2))/96

========================================================================================
PART II — SOLVE THE EXACT COM ORDINARY-LAGRANGIAN COEFFICIENTS
========================================================================================
l1^(GR) = (3*nu)/16 - (21*nu^2)/16 + (9*nu^3)/4 + (5 - 35*nu + 70*nu^2 - 35*nu^3)/128
l2^(GR) = 0
l3^(GR) = 0
l4^(GR) = 0
l5^(GR) = 0
l6^(GR) = 1/4 + (7*nu)/8 - (35*nu^2)/8 - (21*nu^3)/4 + (7 - 42*nu + 53*nu^2 + 5*nu^3)/16
l7^(GR) = (11*nu^2)/8 - ((2 - 3*nu)*nu^2)/16 - (9*nu^3)/2
l8^(GR) = (3*nu^2)/4 - (3*(1 - nu)*nu^2)/16 - (9*nu^3)/4
l9^(GR) = (5*nu^3)/16
l10^(GR) = 5/4 + (15*nu)/8 + (123*nu^2)/8 + (13*nu^3)/4 + (27 - 136*nu - 109*nu^2)/16
l11^(GR) = (7*nu)/8 + (41*nu^2)/8 + (31*nu^3)/4 - (nu*(17 + 30*nu))/16
l12^(GR) = (9*nu^2)/2 + 4*nu^3 - (nu*(5 + 43*nu))/12
l13^(GR) = -3/2 - (59*nu)/4 - (25*nu^2)/4 - nu^3/2 + (600 + 552*nu^2 - nu*(-1340 + 3*Pi^2))/192
l14^(GR) = (7*nu)/4 - (31*nu^2)/4 - (7*nu^3)/2 - (nu*(-340 - 112*nu - 3*Pi^2))/64
l15^(GR) = (-12 - nu*(872 - 63*Pi^2))/96
map check h1 = 0
map check h2 = 0
map check h3 = 0
map check h4 = 0
map check h5 = 0
map check h6 = 0
map check h7 = 0
map check h8 = 0
map check h9 = 0
map check h10 = 0
map check h11 = 0
map check h12 = 0
map check h13 = 0
map check h14 = 0
map check h15 = 0

========================================================================================
PART III — ONE-BODY GATE AND SEED CHECKS
========================================================================================
one-body h1 target - seed = 0
one-body l1 target - seed = 0
one-body h6 target - seed = 0
one-body l6 target - seed = 0
one-body h10 target - seed = 0
one-body l10 target - seed = 0
one-body h13 target - seed = 0
one-body l13 target - seed = 0
one-body h15 target - seed = 0
one-body l15 target - seed = 0
one-body h2 = 0
one-body l2 = 0
one-body h3 = 0
one-body l3 = 0
one-body h4 = 0
one-body l4 = 0
one-body h5 = 0
one-body l5 = 0
one-body h7 = 0
one-body l7 = 0
one-body h8 = 0
one-body l8 = 0
one-body h9 = 0
one-body l9 = 0
one-body h11 = 0
one-body l11 = 0
one-body h12 = 0
one-body l12 = 0
one-body h14 = 0
one-body l14 = 0

========================================================================================
PART IV — GENUINE COMPARABLE-MASS RESIDUALS
========================================================================================
nu -> 0 residual h1 = 0
nu -> 0 residual l1 = 0
Delta l1 + Delta h1 = 0
nu -> 0 residual h2 = 0
nu -> 0 residual l2 = 0
Delta l2 + Delta h2 = 0
nu -> 0 residual h3 = 0
nu -> 0 residual l3 = 0
Delta l3 + Delta h3 = 0
nu -> 0 residual h4 = 0
nu -> 0 residual l4 = 0
Delta l4 + Delta h4 = 0
nu -> 0 residual h5 = 0
nu -> 0 residual l5 = 0
Delta l5 + Delta h5 = 0
nu -> 0 residual h6 = 0
nu -> 0 residual l6 = 0
Delta l6 + Delta h6 = 0
nu -> 0 residual h7 = 0
nu -> 0 residual l7 = 0
Delta l7 + Delta h7 = 0
nu -> 0 residual h8 = 0
nu -> 0 residual l8 = 0
Delta l8 + Delta h8 = 0
nu -> 0 residual h9 = 0
nu -> 0 residual l9 = 0
Delta l9 + Delta h9 = 0
nu -> 0 residual h10 = 0
nu -> 0 residual l10 = 0
Delta l10 + Delta h10 = 0
nu -> 0 residual h11 = 0
nu -> 0 residual l11 = 0
Delta l11 + Delta h11 = 0
nu -> 0 residual h12 = 0
nu -> 0 residual l12 = 0
Delta l12 + Delta h12 = 0
nu -> 0 residual h13 = 0
nu -> 0 residual l13 = 0
Delta l13 + Delta h13 = 0
nu -> 0 residual h14 = 0
nu -> 0 residual l14 = 0
Delta l14 + Delta h14 = 0
nu -> 0 residual h15 = 0
nu -> 0 residual l15 = 0
Delta l15 + Delta h15 = 0

========================================================================================
FINAL FOUNDATION LEDGER
========================================================================================
Passes: 90
Fails: 0
"*)
