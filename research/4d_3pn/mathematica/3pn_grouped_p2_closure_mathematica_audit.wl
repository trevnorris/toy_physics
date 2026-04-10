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

v2 = Symbol["v2"];
u2 = v2 - d^2;

monoms = <|
  1 -> v2^4,
  2 -> v2^3 d^2,
  3 -> v2^2 d^4,
  4 -> v2 d^6,
  5 -> d^8,
  6 -> U v2^3,
  7 -> U v2^2 d^2,
  8 -> U v2 d^4,
  9 -> U d^6,
  10 -> U^2 v2^2,
  11 -> U^2 v2 d^2,
  12 -> U^2 d^4,
  13 -> U^3 v2,
  14 -> U^3 d^2,
  15 -> U^4
|>;

coeffVector[expr_] := Module[{poly = CoefficientRules[Expand[expr], {v2, d, U}]},
  Table[
    Coefficient[expr, monoms[i]],
    {i, 1, 15}
  ]
];

deltaL = <|
  6 -> nu (38 - 116 nu - 57 nu^2)/16,
  7 -> nu^2 (20 - 69 nu)/16,
  8 -> 3 nu^2 (3 - 11 nu)/16,
  9 -> 5 nu^3/16,
  10 -> nu (129 - 98 nu + 52 nu^2)/16,
  11 -> nu (-3 + 52 nu + 124 nu^2)/16,
  12 -> nu (-5 + 11 nu + 48 nu^2)/12,
  13 -> -nu (244 + 3 Pi^2 + 1272 nu + 96 nu^2)/192,
  14 -> nu (452 + 3 Pi^2 - 384 nu - 224 nu^2)/64,
  15 -> nu (-908 + 63 Pi^2)/96
|>;

C20sq = Expand[(3 d^2 - v2 - 2 U)^2/6];
C21sq = Expand[2 d^2 u2];
C22sq = Expand[u2^2/2];

T20 = Expand[U d^2 (3 u2 - U)^2/3];
T21 = Expand[U u2 (u2 - d^2 - U)^2];
T22 = Expand[U d^2 u2^2];

S20 = Expand[U^2 C20sq];
S21 = Expand[U^2 C21sq];
S22 = Expand[U^2 C22sq];

V20 = Expand[v2 S20/U];
V21 = Expand[v2 S21/U];
V22 = Expand[v2 S22/U];

banner["PART I — MINIMAL LOCAL GROUPED-P2 DEMOTION FAILURE"];

LFront = Expand[U^2 (c20 d^2 (3 u2 - U)^2/3 + c21 u2 (u2 - d^2 - U)^2 + c22 d^2 u2^2)];
LDem = Expand[LFront/U];

coeffs = Association@Table[i -> Expand[Coefficient[LDem, monoms[i]]], {i, 1, 15}];
Do[checkZero["pure kinetic slot l" <> ToString[i], coeffs[i]], {i, 1, 5}];
checkZero["slot l15", coeffs[15]];

M = Table[
  {
    Coefficient[coeffs[i], c20],
    Coefficient[coeffs[i], c21],
    Coefficient[coeffs[i], c22]
  },
  {i, 6, 14}
];
checkZero["minimal image rank - 3", MatrixRank[M] - 3];

L6 = Symbol["L6"]; L7 = Symbol["L7"]; L8 = Symbol["L8"]; L9 = Symbol["L9"];
L10 = Symbol["L10"]; L11 = Symbol["L11"]; L12 = Symbol["L12"]; L13 = Symbol["L13"]; L14 = Symbol["L14"];
relations = {
  2 L6 + 2 L7 + L8,
  -L6 - L7 + L9,
  L10 + 2 L6,
  L11 + L12 - 2 L6,
  L13 - L6,
  L14 + L11/6
};
imageSubs = Thread[{L6, L7, L8, L9, L10, L11, L12, L13, L14} -> Table[coeffs[i], {i, 6, 14}]];
Do[checkZero["minimal relation " <> ToString[i], relations[[i]] /. imageSubs], {i, 1, Length[relations]}];

targetSubs = Thread[{L6, L7, L8, L9, L10, L11, L12, L13, L14} -> Table[deltaL[i], {i, 6, 14}]];
Do[
  If[FullSimplify[Expand[relations[[i]] /. targetSubs]] === 0, Throw[$Failed]];
  passCount++;
  Print["target violation ", i, " = ", FullSimplify[Expand[relations[[i]] /. targetSubs]]];,
  {i, 1, Length[relations]}
];

banner["PART II — RICHER (T,S,V) GROUPED CLOSURE"];

names = {T20, T21, T22, S20, S21, S22, V20, V21, V22};
AMid = Table[coeffVector[names[[j]]][[i]], {i, 6, 14}, {j, 1, 9}];
checkZero["det(M_mid) + 4/27", Det[AMid] + 4/27];

targetVec = Table[deltaL[i], {i, 6, 14}];
sol = LinearSolve[AMid, targetVec];
exprTarget = Expand[Sum[sol[[i]] names[[i]], {i, 1, 9}]];
coordsTarget = coeffVector[exprTarget];
Do[checkZero["grouped middle slot l" <> ToString[i], coordsTarget[[i]] - deltaL[i]], {i, 6, 14}];
Do[checkZero["grouped pure kinetic slot l" <> ToString[i], coordsTarget[[i]]], {i, 1, 5}];

l15Pred = FullSimplify[coordsTarget[[15]]];
checkZero[
  "l15 prediction relation",
  l15Pred - (deltaL[10] + deltaL[11] + deltaL[12] + 2 (deltaL[6] + deltaL[7] + deltaL[8] + deltaL[9]))
];
checkZero[
  "remaining static gap",
  deltaL[15] - l15Pred - nu (408 nu^2 + 1232 nu - 2080 + 63 Pi^2)/96
];

banner["FINAL GROUPED-P2 CLOSURE LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — MINIMAL LOCAL GROUPED-P2 DEMOTION FAILURE
========================================================================================
pure kinetic slot l1 = 0
pure kinetic slot l2 = 0
pure kinetic slot l3 = 0
pure kinetic slot l4 = 0
pure kinetic slot l5 = 0
slot l15 = 0
minimal image rank - 3 = 0
minimal relation 1 = 0
minimal relation 2 = 0
minimal relation 3 = 0
minimal relation 4 = 0
minimal relation 5 = 0
minimal relation 6 = 0
target violation 1 = (nu*(76 - 3*nu*(61 + 95*nu)))/16
target violation 2 = (nu*(-38 + nu*(96 + 131*nu)))/16
target violation 3 = -1/16*(nu*(-205 + 330*nu + 62*nu^2))
target violation 4 = (nu*(-257 + 896*nu + 906*nu^2))/48
target violation 5 = (nu*(-700 + 12*nu*(10 + 49*nu) - 3*Pi^2))/192
target violation 6 = (nu*(-8*nu*(131 + 53*nu) + 9*(150 + Pi^2)))/192

========================================================================================
PART II — RICHER (T,S,V) GROUPED CLOSURE
========================================================================================
det(M_mid) + 4/27 = 0
grouped middle slot l6 = 0
grouped middle slot l7 = 0
grouped middle slot l8 = 0
grouped middle slot l9 = 0
grouped middle slot l10 = 0
grouped middle slot l11 = 0
grouped middle slot l12 = 0
grouped middle slot l13 = 0
grouped middle slot l14 = 0
grouped pure kinetic slot l1 = 0
grouped pure kinetic slot l2 = 0
grouped pure kinetic slot l3 = 0
grouped pure kinetic slot l4 = 0
grouped pure kinetic slot l5 = 0
l15 prediction relation = 0
remaining static gap = 0

========================================================================================
FINAL GROUPED-P2 CLOSURE LEDGER
========================================================================================
Passes: 36
Fails: 0
"*)
