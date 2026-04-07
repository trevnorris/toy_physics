(* ::Package:: *)

ClearAll["Global`*"];
$HistoryLength = 0;

SectionLine[title_String] := (
  Print[""];
  Print[StringRepeat["=", 72]];
  Print[title];
  Print[StringRepeat["=", 72]];
  Print[""];
);

ShowValue[label_String, expr_] := Module[{value = FullSimplify[expr]},
  Print[label <> " = " <> ToString[InputForm[value]]];
  value
];

CheckEqual[label_String, expr_, expected_] := Module[{value = FullSimplify[expr]},
  Print[label <> " = " <> ToString[InputForm[value]]];
  If[!TrueQ[FullSimplify[value == expected]],
    Print["FAIL: expected " <> ToString[InputForm[expected]]];
    Exit[1];
  ];
  value
];

SectionLine["Adiabatic kappa_PV closure checks"];

ClearAll[n, x];

dlnFGeneral = FullSimplify[((n - 1)/2 + (-1) x + n ((1 + 2 x)/3))/(1 + x + (1 + 2 x)/3)];
CheckEqual[
  "d ln F / d ln rho (general)",
  dlnFGeneral,
  ((5 n - 3)/2 + (2 n - 3) x)/(4 + 5 x)
];

dlnFN5 = FullSimplify[dlnFGeneral /. n -> 5];
CheckEqual["d ln F / d ln rho (n=5)", dlnFN5, (11 + 7 x)/(4 + 5 x)];

xRule = Solve[dlnFN5 == 5/2, x];
ShowValue["solve d ln F / d ln rho == 5/2", xRule];
CheckEqual["x = Ef/Ew", x /. First[xRule], 2/11];

xSol = x /. First[xRule];
epvOverEw = FullSimplify[(1 + 2 xSol)/3];
fractions = FullSimplify[{1, xSol, epvOverEw}/Total[{1, xSol, epvOverEw}]];
scaledPartition = FullSimplify[11 {1, xSol, epvOverEw}];

CheckEqual["Epv/Ew", epvOverEw, 5/11];
CheckEqual["Ew:Ef:Epv", scaledPartition, {11, 2, 5}];
CheckEqual["fractions", fractions, {11/18, 2/18, 5/18}];

kappaPV = FullSimplify[(dlnFN5 /. x -> xSol) - 1];
beta1PN = FullSimplify[1 + 1/2 + kappaPV];

CheckEqual["kappa_PV", kappaPV, 3/2];
CheckEqual["beta_1PN", beta1PN, 3];

dlnAGeneral = FullSimplify[-((- (n - 1)/2) + 2 x + n (1 + 2 x))/(4 + 10 x)];
CheckEqual["d ln a / d ln rho (general)", dlnAGeneral, -((- (n - 1)/2 + 2 x + n (1 + 2 x))/(4 + 10 x))];
CheckEqual["d ln a / d ln rho (n=5, x=2/11)", dlnAGeneral /. {n -> 5, x -> xSol}, -57/64];

Print["ALL KAPPA_PV CHECKS PASSED"];

(*"
Output:

========================================================================
Adiabatic kappa_PV closure checks
========================================================================

d ln F / d ln rho (general) = (-3 - 6*x + n*(5 + 4*x))/(8 + 10*x)
d ln F / d ln rho (n=5) = (11 + 7*x)/(4 + 5*x)
solve d ln F / d ln rho == 5/2 = {{x -> 2/11}}
x = Ef/Ew = 2/11
Epv/Ew = 5/11
Ew:Ef:Epv = {11, 2, 5}
fractions = {11/18, 1/9, 5/18}
kappa_PV = 3/2
beta_1PN = 3
d ln a / d ln rho (general) = -(((1 + n)*(1 + 4*x))/(8 + 20*x))
d ln a / d ln rho (n=5, x=2/11) = -57/64
ALL KAPPA_PV CHECKS PASSED
"*)
