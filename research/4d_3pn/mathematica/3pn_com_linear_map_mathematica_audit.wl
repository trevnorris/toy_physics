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

banner["PART I — EXACT COM 3PN LINEAR MAP"];

hExpectedList = {
  -l1c + 3 nu/16 - 21 nu^2/16 + 9 nu^3/4,
  -l2c,
  -l3c,
  -l4c,
  -l5c,
  -l6c + 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4,
  -l7c + 11 nu^2/8 - 9 nu^3/2,
  -l8c + 3 nu^2/4 - 9 nu^3/4,
  -l9c,
  -l10c + 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4,
  -l11c + 7 nu/8 + 41 nu^2/8 + 31 nu^3/4,
  -l12c + 9 nu^2/2 + 4 nu^3,
  -l13c - 3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2,
  -l14c + 7 nu/4 - 31 nu^2/4 - 7 nu^3/2,
  -l15c
};

Do[Print["h", i, " = ", hExpectedList[[i]]], {i, 1, 15}];

banner["PART II — INVERSE MAP l_i(h_j)"];

inverseMapList = {
  -h1 + 3 nu/16 - 21 nu^2/16 + 9 nu^3/4,
  -h2,
  -h3,
  -h4,
  -h5,
  -h6 + 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4,
  -h7 + 11 nu^2/8 - 9 nu^3/2,
  -h8 + 3 nu^2/4 - 9 nu^3/4,
  -h9,
  -h10 + 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4,
  -h11 + 7 nu/8 + 41 nu^2/8 + 31 nu^3/4,
  -h12 + 9 nu^2/2 + 4 nu^3,
  -h13 - 3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2,
  -h14 + 7 nu/4 - 31 nu^2/4 - 7 nu^3/2,
  -h15
};

Do[Print["l", i, " = ", inverseMapList[[i]]], {i, 1, 15}];

Do[
  checkZero[
    "inverse check h" <> ToString[i],
    (hExpectedList[[i]] /. {
      l1c -> inverseMapList[[1]], l2c -> inverseMapList[[2]], l3c -> inverseMapList[[3]],
      l4c -> inverseMapList[[4]], l5c -> inverseMapList[[5]], l6c -> inverseMapList[[6]],
      l7c -> inverseMapList[[7]], l8c -> inverseMapList[[8]], l9c -> inverseMapList[[9]],
      l10c -> inverseMapList[[10]], l11c -> inverseMapList[[11]], l12c -> inverseMapList[[12]],
      l13c -> inverseMapList[[13]], l14c -> inverseMapList[[14]], l15c -> inverseMapList[[15]]
    }) - Symbol["h" <> ToString[i]]
  ],
  {i, 1, 15}
];

banner["PART III — COM IMAGE OF THE CARRIED 3PN SELF/STATIC SEED"];

seedMap = {
  l1c -> 5/128 - 35 nu/128 + 35 nu^2/64 - 35 nu^3/128,
  l2c -> 0, l3c -> 0, l4c -> 0, l5c -> 0,
  l6c -> 11/16 - 33 nu/8 + 99 nu^2/16 - 11 nu^3/8,
  l7c -> 0, l8c -> 0, l9c -> 0,
  l10c -> 47/16 - 235 nu/16 + 235 nu^2/16,
  l11c -> 0, l12c -> 0,
  l13c -> 13/8 - 13 nu/2 + 13 nu^2/4,
  l14c -> 0,
  l15c -> -1/8 + 3 nu/8
};

hSeedExpected = {
  -5/128 + 59 nu/128 - 119 nu^2/64 + 323 nu^3/128,
  0,
  0,
  0,
  0,
  -7/16 + 5 nu - 169 nu^2/16 - 31 nu^3/8,
  nu^2 (11 - 36 nu)/8,
  3 nu^2 (1 - 3 nu)/4,
  0,
  -27/16 + 265 nu/16 + 11 nu^2/16 + 13 nu^3/4,
  nu (7 + 41 nu + 62 nu^2)/8,
  nu^2 (9 + 8 nu)/2,
  -25/8 - 33 nu/4 - 19 nu^2/2 - nu^3/2,
  nu (7 - 31 nu - 14 nu^2)/4,
  (1 - 3 nu)/8
};

Do[
  checkZero["seed image h" <> ToString[i], (hExpectedList[[i]] /. seedMap) - hSeedExpected[[i]]],
  {i, 1, 15}
];

banner["FINAL FOUNDATION LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — EXACT COM 3PN LINEAR MAP
========================================================================================
h1 = -l1c + (3*nu)/16 - (21*nu^2)/16 + (9*nu^3)/4
h2 = -l2c
h3 = -l3c
h4 = -l4c
h5 = -l5c
h6 = 1/4 - l6c + (7*nu)/8 - (35*nu^2)/8 - (21*nu^3)/4
h7 = -l7c + (11*nu^2)/8 - (9*nu^3)/2
h8 = -l8c + (3*nu^2)/4 - (9*nu^3)/4
h9 = -l9c
h10 = 5/4 - l10c + (15*nu)/8 + (123*nu^2)/8 + (13*nu^3)/4
h11 = -l11c + (7*nu)/8 + (41*nu^2)/8 + (31*nu^3)/4
h12 = -l12c + (9*nu^2)/2 + 4*nu^3
h13 = -3/2 - l13c - (59*nu)/4 - (25*nu^2)/4 - nu^3/2
h14 = -l14c + (7*nu)/4 - (31*nu^2)/4 - (7*nu^3)/2
h15 = -l15c

========================================================================================
PART II — INVERSE MAP l_i(h_j)
========================================================================================
l1 = -h1 + (3*nu)/16 - (21*nu^2)/16 + (9*nu^3)/4
l2 = -h2
l3 = -h3
l4 = -h4
l5 = -h5
l6 = 1/4 - h6 + (7*nu)/8 - (35*nu^2)/8 - (21*nu^3)/4
l7 = -h7 + (11*nu^2)/8 - (9*nu^3)/2
l8 = -h8 + (3*nu^2)/4 - (9*nu^3)/4
l9 = -h9
l10 = 5/4 - h10 + (15*nu)/8 + (123*nu^2)/8 + (13*nu^3)/4
l11 = -h11 + (7*nu)/8 + (41*nu^2)/8 + (31*nu^3)/4
l12 = -h12 + (9*nu^2)/2 + 4*nu^3
l13 = -3/2 - h13 - (59*nu)/4 - (25*nu^2)/4 - nu^3/2
l14 = -h14 + (7*nu)/4 - (31*nu^2)/4 - (7*nu^3)/2
l15 = -h15
inverse check h1 = 0
inverse check h2 = 0
inverse check h3 = 0
inverse check h4 = 0
inverse check h5 = 0
inverse check h6 = 0
inverse check h7 = 0
inverse check h8 = 0
inverse check h9 = 0
inverse check h10 = 0
inverse check h11 = 0
inverse check h12 = 0
inverse check h13 = 0
inverse check h14 = 0
inverse check h15 = 0

========================================================================================
PART III — COM IMAGE OF THE CARRIED 3PN SELF/STATIC SEED
========================================================================================
seed image h1 = 0
seed image h2 = 0
seed image h3 = 0
seed image h4 = 0
seed image h5 = 0
seed image h6 = 0
seed image h7 = 0
seed image h8 = 0
seed image h9 = 0
seed image h10 = 0
seed image h11 = 0
seed image h12 = 0
seed image h13 = 0
seed image h14 = 0
seed image h15 = 0

========================================================================================
FINAL FOUNDATION LEDGER
========================================================================================
Passes: 30
Fails: 0
"*)
