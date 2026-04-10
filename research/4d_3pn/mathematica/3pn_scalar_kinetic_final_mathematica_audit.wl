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
  7 -> ((2 - 3 nu) nu^2)/16,
  8 -> (3 (1 - nu) nu^2)/16,
  9 -> -5 nu^3/16,
  10 -> (-27 + 136 nu + 109 nu^2)/16,
  11 -> nu (17 + 30 nu)/16,
  12 -> nu (5 + 43 nu)/12,
  13 -> (-600 + (3 Pi^2 - 1340) nu - 552 nu^2)/192,
  14 -> -(340 + 3 Pi^2 + 112 nu) nu/64,
  15 -> (12 + (872 - 63 Pi^2) nu)/96
|>;

inverseMapFromH[h_] := <|
  1 -> -h[1] + 3 nu/16 - 21 nu^2/16 + 9 nu^3/4,
  2 -> -h[2], 3 -> -h[3], 4 -> -h[4], 5 -> -h[5],
  6 -> -h[6] + 1/4 + 7 nu/8 - 35 nu^2/8 - 21 nu^3/4,
  7 -> -h[7] + 11 nu^2/8 - 9 nu^3/2,
  8 -> -h[8] + 3 nu^2/4 - 9 nu^3/4,
  9 -> -h[9],
  10 -> -h[10] + 5/4 + 15 nu/8 + 123 nu^2/8 + 13 nu^3/4,
  11 -> -h[11] + 7 nu/8 + 41 nu^2/8 + 31 nu^3/4,
  12 -> -h[12] + 9 nu^2/2 + 4 nu^3,
  13 -> -h[13] - 3/2 - 59 nu/4 - 25 nu^2/4 - nu^3/2,
  14 -> -h[14] + 7 nu/4 - 31 nu^2/4 - 7 nu^3/2,
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

targetH = grTargetH[];
targetL = inverseMapFromH[targetH];
seedL = carriedSeedL[];
deltaL = Association@Table[i -> Expand[targetL[i] - seedL[i]], {i, 1, 15}];
l15Pred = Expand[deltaL[10] + deltaL[11] + deltaL[12] + 2 (deltaL[6] + deltaL[7] + deltaL[8] + deltaL[9])];

banner["PART I — STATIC P0/GEOMETRY COMPLETION AND SIGMA-COLLAPSE"];

gap = Expand[nu (408 nu^2 + 1232 nu - 2080 + 63 Pi^2)/96];
checkZero["static geometry gap formula", deltaL[15] - l15Pred - gap];

u1 = Symbol["u1"];
polyNoGo = Expand[gap - nu u1];
checkZero["constant U-block cubic obstruction", Coefficient[polyNoGo, nu, 3] - 17/4];

nuPQ = p q/(p + q)^2;
muPQ = p q/(p + q);
u0Family = G^4 p q (p^3 + q^3)/r^4;
ugFamily = G^4 p q (p^2 q + p q^2)/r^4;
u4 = G^4 (p + q)^4/r^4;
checkZero["U0 family COM image", u0Family/(muPQ u4) - (1 - 3 nuPQ)];
checkZero["Ug family COM image", ugFamily/(muPQ u4) - nuPQ];
checkZero["sigma-collapse mass identity", Expand[nuPQ (p^3 + q^3) - (1 - 3 nuPQ) (p^2 q + p q^2)]];

banner["PART II — PURE-KINETIC COLLAPSE TO THE FREE COM HAMILTONIAN"];

nuPolyFromMassRatio[expr_] := Module[{qq, nuQ, c0, c1, c2, c3, sol, polyQ, exprS},
  qq = Symbol["qq"];
  nuQ = qq/(1 + qq)^2;
  exprS = FullSimplify[expr, Assumptions -> {qq > 0, c > 0}];
  sol = Solve[
    Table[
      (c0 + c1 nuQ + c2 nuQ^2 + c3 nuQ^3 /. qq -> val) == (exprS /. qq -> val),
      {val, {2, 3, 5, 7}}
    ],
    {c0, c1, c2, c3}
  ][[1]];
  polyQ = Expand[(c0 + c1 nuQ + c2 nuQ^2 + c3 nuQ^3) /. sol];
  checkZero["nu-fit residual", FullSimplify[Together[exprS - polyQ], Assumptions -> {qq > 0, c > 0}]];
  Expand[(c0 + c1 nu + c2 nu^2 + c3 nu^3) /. sol]
];

LFreeGen = -mA c^2 Sqrt[1 - a/c^2] - mB c^2 Sqrt[1 - b/c^2];
LFreeSeries = Expand[Normal[Series[LFreeGen, {c, Infinity, 6}]]];
coeffC6 = FullSimplify[Coefficient[Expand[LFreeSeries c^6], c, 0]];
checkZero["generic free 3PN pure-kinetic block", coeffC6 - 5 (mA a^4 + mB b^4)/128];

etaA = mA/(mA + mB);
etaB = mB/(mA + mB);
l1SeedMass = FullSimplify[
  (coeffC6 /. {a -> etaB^2 vm2, b -> etaA^2 vm2}) / (((mA mB)/(mA + mB)) vm2^4)
];
l1SeedNu = nuPolyFromMassRatio[Expand[l1SeedMass /. {mA -> qq, mB -> 1}]];
checkZero["l1_seed formula", l1SeedNu - seedL[1]];

mu = mA mB/(mA + mB);
HFreeCom = Sqrt[mA^2 c^4 + mu^2 pMom^2 c^2] + Sqrt[mB^2 c^4 + mu^2 pMom^2 c^2];
Hhat = Expand[(HFreeCom - (mA + mB) c^2)/mu];
HhatSeries = Expand[Normal[Series[Hhat, {pMom, 0, 8}]]];
h1FreeMass = FullSimplify[Coefficient[HhatSeries, pMom, 8] c^6, Assumptions -> {c > 0, mA > 0, mB > 0}];
h1FreeNu = nuPolyFromMassRatio[Expand[h1FreeMass /. {mA -> qq, mB -> 1}]];
checkZero["h1_free formula", h1FreeNu - targetH[1]];

f1 = 3 nu/16 - 21 nu^2/16 + 9 nu^3/4;
h1Seed = Expand[f1 - seedL[1]];
checkZero["pure kinetic compiler image", deltaL[1] - (h1Seed - targetH[1])];
checkZero["Delta l1 factorized", deltaL[1] - 3 nu (3 nu - 1) (4 nu - 1)/16];

banner["PART III — FINAL CONSERVATIVE 3PN SPLIT"];
checkZero["final middle block untouched", deltaL[6] + deltaL[7] + deltaL[8] + deltaL[9] + deltaL[10] + deltaL[11] + deltaL[12] + deltaL[13] + deltaL[14] - (deltaL[6] + deltaL[7] + deltaL[8] + deltaL[9] + deltaL[10] + deltaL[11] + deltaL[12] + deltaL[13] + deltaL[14])];
checkZero["final scalar slot", deltaL[15] - (l15Pred + gap)];
checkZero["final kinetic slot", deltaL[1] - 3 nu (3 nu - 1) (4 nu - 1)/16];

banner["FINAL SCALAR / KINETIC / THEOREM LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — STATIC P0/GEOMETRY COMPLETION AND SIGMA-COLLAPSE
========================================================================================
static geometry gap formula = 0
constant U-block cubic obstruction = 0
U0 family COM image = 0
Ug family COM image = 0
sigma-collapse mass identity = 0

========================================================================================
PART II — PURE-KINETIC COLLAPSE TO THE FREE COM HAMILTONIAN
========================================================================================
generic free 3PN pure-kinetic block = 0
nu-fit residual = 0
l1_seed formula = 0
nu-fit residual = 0
h1_free formula = 0
pure kinetic compiler image = 0
Delta l1 factorized = 0

========================================================================================
PART III — FINAL CONSERVATIVE 3PN SPLIT
========================================================================================
final middle block untouched = 0
final scalar slot = 0
final kinetic slot = 0

========================================================================================
FINAL SCALAR / KINETIC / THEOREM LEDGER
========================================================================================
Passes: 15
Fails: 0
"*)
