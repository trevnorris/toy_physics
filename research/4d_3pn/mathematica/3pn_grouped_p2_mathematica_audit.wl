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

banner["GROUPED REAL P2 EXPANSIONS"];

y20 = 1 + u220 omega^2 + u420 omega^4;
y21 = 1 + u221 omega^2 + u421 omega^4;
y22 = 1 + u222 omega^2 + u422 omega^4;

Print["Y20(omega) = ", y20];
Print["Y21(omega) = ", y21];
Print["Y22(omega) = ", y22];

banner["EXACT GROUPED -> AXISYMMETRIC INVERSE MAP"];

ubar2 = (u220 + 2 u221 + 2 u222)/5;
a2 = (2 u220 - u221 - u222)/10;
b2 = (u221 - u222)/2;

ubar4 = (u420 + 2 u421 + 2 u422)/5;
a4 = (2 u420 - u421 - u422)/10;
b4 = (u421 - u422)/2;

Print["ubar2 = ", ubar2];
Print["a2 = ", a2];
Print["b2 = ", b2];
Print["ubar4 = ", ubar4];
Print["a4 = ", a4];
Print["b4 = ", b4];

u220Fwd = ubar2Sym + 4 a2Sym;
u221Fwd = ubar2Sym - a2Sym + b2Sym;
u222Fwd = ubar2Sym - a2Sym - b2Sym;

u420Fwd = ubar4Sym + 4 a4Sym;
u421Fwd = ubar4Sym - a4Sym + b4Sym;
u422Fwd = ubar4Sym - a4Sym - b4Sym;

checkZero["u2^(20) recovered", (u220Fwd /. {ubar2Sym -> ubar2, a2Sym -> a2, b2Sym -> b2}) - u220];
checkZero["u2^(21) recovered", (u221Fwd /. {ubar2Sym -> ubar2, a2Sym -> a2, b2Sym -> b2}) - u221];
checkZero["u2^(22) recovered", (u222Fwd /. {ubar2Sym -> ubar2, a2Sym -> a2, b2Sym -> b2}) - u222];

checkZero["u4^(20) recovered", (u420Fwd /. {ubar4Sym -> ubar4, a4Sym -> a4, b4Sym -> b4}) - u420];
checkZero["u4^(21) recovered", (u421Fwd /. {ubar4Sym -> ubar4, a4Sym -> a4, b4Sym -> b4}) - u421];
checkZero["u4^(22) recovered", (u422Fwd /. {ubar4Sym -> ubar4, a4Sym -> a4, b4Sym -> b4}) - u422];

banner["WEIGHTED-SUM CONSTRAINTS AND ANISOTROPY NORMS"];

checkZero[
  "weighted sum constraint at O(omega^2)",
  (u220 - ubar2) + 2 (u221 - ubar2) + 2 (u222 - ubar2)
];
checkZero[
  "weighted sum constraint at O(omega^4)",
  (u420 - ubar4) + 2 (u421 - ubar4) + 2 (u422 - ubar4)
];

a2Sq = FullSimplify[4 a2^2 + (4/5) b2^2];
a4Sq = FullSimplify[4 a4^2 + (4/5) b4^2];
Print["A2^2 = ", a2Sq];
Print["A4^2 = ", a4Sq];

banner["MINIMAL ISOTROPIC BRANCH FORMULAS"];

omegaQSq = FullSimplify[1/(4 u2)];
gamma5Norm = FullSimplify[9 u2^(5/2)];
k0barTarget = FullSimplify[2 G/(45 c^5 u2^(5/2))];

Print["Omega_Q^2 = ", omegaQSq];
Print["Gamma5_norm = ", gamma5Norm];
Print["K0bar_target = ", k0barTarget];
Print["single-pole 5PN branch law: u4 = 4 u2^2"];

banner["FINAL GROUPED-P2 KICKOFF LEDGER"];
Print["First 3PN isotropy gate: a2 = 0 and b2 = 0"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
GROUPED REAL P2 EXPANSIONS
========================================================================================
Y20(omega) = 1 + omega^2*u220 + omega^4*u420
Y21(omega) = 1 + omega^2*u221 + omega^4*u421
Y22(omega) = 1 + omega^2*u222 + omega^4*u422

========================================================================================
EXACT GROUPED -> AXISYMMETRIC INVERSE MAP
========================================================================================
ubar2 = (u220 + 2*u221 + 2*u222)/5
a2 = (2*u220 - u221 - u222)/10
b2 = (u221 - u222)/2
ubar4 = (u420 + 2*u421 + 2*u422)/5
a4 = (2*u420 - u421 - u422)/10
b4 = (u421 - u422)/2
u2^(20) recovered = 0
u2^(21) recovered = 0
u2^(22) recovered = 0
u4^(20) recovered = 0
u4^(21) recovered = 0
u4^(22) recovered = 0

========================================================================================
WEIGHTED-SUM CONSTRAINTS AND ANISOTROPY NORMS
========================================================================================
weighted sum constraint at O(omega^2) = 0
weighted sum constraint at O(omega^4) = 0
A2^2 = (5*(u221 - u222)^2 + (-2*u220 + u221 + u222)^2)/25
A4^2 = (5*(u421 - u422)^2 + (-2*u420 + u421 + u422)^2)/25

========================================================================================
MINIMAL ISOTROPIC BRANCH FORMULAS
========================================================================================
Omega_Q^2 = 1/(4*u2)
Gamma5_norm = 9*u2^(5/2)
K0bar_target = (2*G)/(45*c^5*u2^(5/2))
single-pole 5PN branch law: u4 = 4 u2^2

========================================================================================
FINAL GROUPED-P2 KICKOFF LEDGER
========================================================================================
First 3PN isotropy gate: a2 = 0 and b2 = 0
Passes: 8
Fails: 0
"*)
