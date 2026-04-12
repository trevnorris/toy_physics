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

banner["QUARTIC-ORDER PERTURBATIVE LEGENDRE COMPILER"];

L0 = (m v^2)/2;
L1 = a1 v + a2 v^2 + a3 v^3;
L2 = b1 v + b2 v^2;
L3 = c1 v + c2 v^2;
L4 = d1 v + d2 v^2;
L = L0 + eps L1 + eps^2 L2 + eps^3 L3 + eps^4 L4;

v0 = p/m;
vSeries = v0 + eps w1 + eps^2 w2 + eps^3 w3;
peq = Expand[(D[L, v] /. v -> vSeries) - p];
eq1 = Coefficient[peq, eps, 1];
eq2 = Coefficient[peq, eps, 2];
eq3 = Coefficient[peq, eps, 3];

sol1 = w1 /. First @ Solve[eq1 == 0, w1];
sol2 = w2 /. First @ Solve[(eq2 /. w1 -> sol1) == 0, w2];
sol3 = w3 /. First @ Solve[(eq3 /. {w1 -> sol1, w2 -> sol2}) == 0, w3];

Print["v1 = ", FullSimplify[sol1]];
Print["v2 = ", FullSimplify[sol2]];
Print["v3 = ", FullSimplify[sol3]];

hexact = Expand[p (v0 + eps sol1 + eps^2 sol2 + eps^3 sol3) - (L /. v -> (v0 + eps sol1 + eps^2 sol2 + eps^3 sol3))];
hseries = Expand @ Normal @ Series[hexact, {eps, 0, 4}];
h1Exact = Coefficient[hseries, eps, 1];
h2Exact = Coefficient[hseries, eps, 2];
h3Exact = Coefficient[hseries, eps, 3];
h4Exact = Coefficient[hseries, eps, 4];

A0 = D[L1, v] /. v -> v0;
B0 = D[L2, v] /. v -> v0;
D0 = D[L3, v] /. v -> v0;
C0 = D[L1, {v, 2}] /. v -> v0;
E0 = D[L2, {v, 2}] /. v -> v0;
T0 = D[L1, {v, 3}] /. v -> v0;

h1Formula = Expand[-(L1 /. v -> v0)];
h2Formula = Expand[-(L2 /. v -> v0) + A0^2/(2 m)];
h3Formula = Expand[-(L3 /. v -> v0) + A0 B0/m - A0^2 C0/(2 m^2)];
h4Formula = Expand[
  -(L4 /. v -> v0) +
  A0 D0/m +
  B0^2/(2 m) -
  B0 C0 A0/m^2 -
  A0^2 E0/(2 m^2) +
  A0^2 C0^2/(2 m^3) +
  A0^3 T0/(6 m^3)
];

checkZero["H1 exact - formula", h1Exact - h1Formula];
checkZero["H2 exact - formula", h2Exact - h2Formula];
checkZero["H3 exact - formula", h3Exact - h3Formula];
checkZero["H4 exact - formula", h4Exact - h4Formula];

banner["FINAL QUARTIC COMPILER LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
QUARTIC-ORDER PERTURBATIVE LEGENDRE COMPILER
========================================================================================
v1 = -((a1*m^2 + 2*a2*m*p + 3*a3*p^2)/m^3)
v2 = (m^3*(2*a1*a2 - b1*m) + 2*m^2*(2*a2^2 + 3*a1*a3 - b2*m)*p + 18*a2*a3*m*p^2 + 18*a3^2*p^3)/m^5
v3 = -((m^4*(a1*(4*a2^2 + 3*a1*a3) - 2*(a2*b1 + a1*b2)*m + c1*m^2) + 2*m^3*(4*a2^3 + 18*a1*a2*a3 - 3*a3*b1*m - 4*a2*b2*m + c2*m^2)*p + 18*a3*m^2*(4*a2^2 + 3*a1*a3 - b2*m)*p^2 + 180*a2*a3^2*m*p^3 + 135*a3^3*p^4)/m^7)
H1 exact - formula = 0
H2 exact - formula = 0
H3 exact - formula = 0
H4 exact - formula = 0

========================================================================================
FINAL QUARTIC COMPILER LEDGER
========================================================================================
Passes: 4
Fails: 0
"*)
