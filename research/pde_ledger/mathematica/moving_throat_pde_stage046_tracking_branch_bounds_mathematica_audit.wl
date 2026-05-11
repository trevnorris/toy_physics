ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectPositiveCoefficients[name_String, expr_, vars_List] := Module[{coeffs},
  coeffs = Last /@ CoefficientRules[Expand[expr], vars];
  Print[name, " coefficients = ", fmt[coeffs]];
  If[AllTrue[coeffs, Positive], pass[name], fail[name, coeffs]];
];

banner["STAGE 029 — TRACKING-BRANCH BOUNDS"];

Clear[xi, delta, r];
$Assumptions =
  Element[{xi, delta, r}, Reals] &&
  0 < xi < 1 && delta > 0 && r > 0;

gTr = FullSimplify[9*xi*(delta + xi)/(9*delta + (9 + 2*r^2)*xi), Assumptions -> $Assumptions];
fTr = FullSimplify[
  (9*delta + (9 + 2*r^2)*xi)^2*(9*delta + (9 + 2*r)*xi)^2/
    (81*(1 - xi)*(9*delta^2 + 18*delta*xi + (9 + 2*r^2)*xi^2)^2),
  Assumptions -> $Assumptions
];
gFlat = FullSimplify[gTr /. r -> 1, Assumptions -> $Assumptions];
fFlat = FullSimplify[fTr /. r -> 1, Assumptions -> $Assumptions];

Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["G_flat = ", fmt[gFlat]];
Print["F_flat = ", fmt[fFlat]];
expectZero["strong-split endpoint for G", FullSimplify[gTr /. r -> 0, Assumptions -> 0 < xi < 1 && delta > 0] - xi];
expectZero["strong-split endpoint for F", FullSimplify[fTr /. r -> 0, Assumptions -> 0 < xi < 1 && delta > 0] - 1/(1 - xi)];

dGdR = Factor[D[gTr, r]];
dFdR = Factor[D[fTr, r]];
pR = 4*r^4*xi^3 + 54*r^2*delta^2*xi + 90*r^2*delta*xi^2 + 36*r^2*xi^3 +
  162*r*delta^3 + 324*r*delta^2*xi + 162*r*delta*xi^2 + 81*delta^3 +
  243*delta^2*xi + 243*delta*xi^2 + 81*xi^3;
dGExpected = -36*r*xi^2*(delta + xi)/(2*r^2*xi + 9*delta + 9*xi)^2;
dFExpected = 4*xi*(2*r*xi + 9*delta + 9*xi)*(2*r^2*xi + 9*delta + 9*xi)*pR/
  (81*(1 - xi)*(2*r^2*xi^2 + 9*delta^2 + 18*delta*xi + 9*xi^2)^3);

Print["dG_tr/dR = ", fmt[dGdR]];
Print["dF_tr/dR = ", fmt[dFdR]];
expectZero["dG_tr/dR formula", dGdR - dGExpected];
expectZero["dF_tr/dR formula", dFdR - dFExpected];
expectPositiveCoefficients["P_R", pR, {r, delta, xi}];

deltaG = Factor[FullSimplify[gTr - gFlat, Assumptions -> $Assumptions]];
deltaF = Factor[FullSimplify[fFlat - fTr, Assumptions -> $Assumptions]];
p1 = 18*r^2*delta^2*xi + 36*r^2*delta*xi^2 + 22*r^2*xi^3 + 81*r*delta^3 +
  180*r*delta^2*xi + 99*r*delta*xi^2 + 162*delta^3 + 423*delta^2*xi +
  360*delta*xi^2 + 99*xi^3;
p2 = 18*r^3*delta^2*xi^2 + 36*r^3*delta*xi^3 + 22*r^3*xi^4 +
  81*r^2*delta^3*xi + 324*r^2*delta^2*xi^2 + 459*r^2*delta*xi^3 +
  220*r^2*xi^4 + 81*r*delta^3*xi + 243*r*delta^2*xi^2 + 261*r*delta*xi^3 +
  99*r*xi^4 + 729*delta^4 + 3078*delta^3*xi + 4959*delta^2*xi^2 +
  3600*delta*xi^3 + 990*xi^4;
gDiffExpected = 18*xi^2*(1 - r^2)*(delta + xi)/((9*delta + 11*xi)*(2*r^2*xi + 9*delta + 9*xi));
fDiffExpected = 4*xi*(1 - r)*p1*p2/
  (81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2*(2*r^2*xi^2 + 9*delta^2 + 18*delta*xi + 9*xi^2)^2);

Print["G_tr - G_flat = ", fmt[deltaG]];
Print["F_flat - F_tr = ", fmt[deltaF]];
expectZero["G_tr - G_flat formula", deltaG - gDiffExpected];
expectZero["F_flat - F_tr formula", deltaF - fDiffExpected];
expectPositiveCoefficients["P1", p1, {r, delta, xi}];
expectPositiveCoefficients["P2", p2, {r, delta, xi}];

Print[""];
Print["Stage 046 Mathematica audit passed."];

Exit[0];
