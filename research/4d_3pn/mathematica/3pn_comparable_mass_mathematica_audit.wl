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

swapFull[expr_] := Expand[expr /. {a -> b, b -> a, c -> c, d -> e, e -> d, p -> q, q -> p}];

canonicalSym[expr_] := Module[
  {s, coeffRules, denLCM, ints, g, sNorm, orderedTerms},
  s = Expand[expr + swapFull[expr]];
  coeffRules = CoefficientRules[s, {p, q, a, b, c, d, e}];
  denLCM = Apply[LCM, Join[{1}, Denominator /@ coeffRules[[All, 2]]]];
  ints = Numerator /@ (denLCM coeffRules[[All, 2]]);
  g = Fold[GCD, 0, Abs /@ ints];
  If[g == 0, g = 1];
  sNorm = Expand[s denLCM/g];
  orderedTerms = SortBy[CoefficientRules[sNorm, {p, q, a, b, c, d, e}], First];
  If[orderedTerms =!= {} && Last[orderedTerms][[2]] < 0, sNorm = -sNorm];
  Expand[sNorm]
];

generateBasis[massDeg_, velDeg_] := Module[
  {basis = <||>, expr, sym, tm, mp, mq, pa, pb, pc, pd, pe, key},
  Do[
    mq = massDeg - mp;
    Do[
      If[2 pa + 2 pb + 2 pc + pd + pe != velDeg, Continue[]];
      expr = p^mp q^mq a^pa b^pb c^pc d^pd e^pe;
      sym = canonicalSym[expr];
      tm = Expand[sym /. {b -> 0, c -> 0, e -> 0, p -> 0, q -> 1}];
      If[tm =!= 0, Continue[]];
      key = ToString[InputForm[sym]];
      basis[key] = sym;,
      {pa, 0, 4}, {pb, 0, 4}, {pc, 0, 4}, {pd, 0, velDeg}, {pe, 0, velDeg}
    ],
    {mp, 0, massDeg}
  ];
  SortBy[Values[basis], ToString[InputForm[#]] &]
];

banner["PART I — CUBIC-ORDER PERTURBATIVE LEGENDRE TRANSFORM"];

l0 = 1/2 m v^2;
l1 = a1 v^3 + a2 v^4;
l2 = b1 v^2 + b2 v^3;
l3 = c1 v + c2 v^2;
l = l0 + eps l1 + eps^2 l2 + eps^3 l3;

v0 = pSym/m;
vSeries = v0 + eps w1 + eps^2 w2;
pEq = Expand[(D[l, v] /. v -> vSeries) - pSym];
eq1 = Coefficient[pEq, eps, 1];
eq2 = Coefficient[pEq, eps, 2];

solW1 = w1 /. First @ Solve[eq1 == 0, w1];
solW2 = w2 /. First @ Solve[Expand[eq2 /. w1 -> solW1] == 0, w2];

hexact = Expand[
  pSym (v0 + eps solW1 + eps^2 solW2) -
  (l /. v -> (v0 + eps solW1 + eps^2 solW2))
];
hseries = Expand @ Normal @ Series[hexact, {eps, 0, 3}];

h1Exact = Expand[Coefficient[hseries, eps, 1]];
h2Exact = Expand[Coefficient[hseries, eps, 2]];
h3Exact = Expand[Coefficient[hseries, eps, 3]];

a0 = D[l1, v] /. v -> v0;
b0 = D[l2, v] /. v -> v0;
c0 = D[l1, {v, 2}] /. v -> v0;

h1Formula = Expand[-(l1 /. v -> v0)];
h2Formula = Expand[-(l2 /. v -> v0) + a0^2/(2 m)];
h3Formula = Expand[-(l3 /. v -> v0) + a0 b0/m - a0^2 c0/(2 m^2)];

checkZero["H1 exact - formula", h1Exact - h1Formula];
checkZero["H2 exact - formula", h2Exact - h2Formula];
checkZero["H3 exact - formula", h3Exact - h3Formula];

banner["PART II — NATURAL 3PN SELF/STATIC SEED"];

l3Seed = Expand[
  5 (mA vA2^4 + mB vB2^4)/128
  + 11 G mA mB (vA2^3 + vB2^3)/(16 r)
  + 47 G^2 mA mB (mB vA2^2 + mA vB2^2)/(16 r^2)
  + 13 G^3 mA mB (mB^2 vA2 + mA^2 vB2)/(8 r^3)
  - G^4 mA mB (mB^3 + mA^3)/(8 r^4)
];

Print["L3_seed = ", l3Seed];

banner["PART III — OVERCOMPLETE COMPARABLE-MASS RESIDUAL BASIS"];

g1Basis = generateBasis[0, 6];
g2Basis = generateBasis[1, 4];
g3Basis = generateBasis[2, 2];
g4Basis = generateBasis[3, 0];

checkZero["G/r sextic count - 24", Length[g1Basis] - 24];
checkZero["G^2/r^2 quartic count - 17", Length[g2Basis] - 17];
checkZero["G^3/r^3 quadratic count - 8", Length[g3Basis] - 8];
checkZero["G^4/r^4 static count - 1", Length[g4Basis] - 1];

checkZero["first Q basis term", First[g1Basis] - (a^2 b + a b^2)];
checkZero["first T basis term", First[g2Basis] - (a^2 p + b^2 q)];
checkZero["first S basis term", First[g3Basis] - (a p^2 + b q^2)];
checkZero["first U basis term", First[g4Basis] - (p^2 q + p q^2)];

banner["FINAL FOUNDATION LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — CUBIC-ORDER PERTURBATIVE LEGENDRE TRANSFORM
========================================================================================
H1 exact - formula = 0
H2 exact - formula = 0
H3 exact - formula = 0

========================================================================================
PART II — NATURAL 3PN SELF/STATIC SEED
========================================================================================
L3_seed = -1/8*(G^4*mA^4*mB)/r^4 - (G^4*mA*mB^4)/(8*r^4) + (13*G^3*mA*mB^3*vA2)/(8*r^3) + (47*G^2*mA*mB^2*vA2^2)/(16*r^2) + (11*G*mA*mB*vA2^3)/(16*r) + (5*mA*vA2^4)/128 + (13*G^3*mA^3*mB*vB2)/(8*r^3) + (47*G^2*mA^2*mB*vB2^2)/(16*r^2) + (11*G*mA*mB*vB2^3)/(16*r) + (5*mB*vB2^4)/128

========================================================================================
PART III — OVERCOMPLETE COMPARABLE-MASS RESIDUAL BASIS
========================================================================================
G/r sextic count - 24 = 0
G^2/r^2 quartic count - 17 = 0
G^3/r^3 quadratic count - 8 = 0
G^4/r^4 static count - 1 = 0
first Q basis term = 0
first T basis term = 0
first S basis term = 0
first U basis term = 0

========================================================================================
FINAL FOUNDATION LEDGER
========================================================================================
Passes: 11
Fails: 0
"*)
