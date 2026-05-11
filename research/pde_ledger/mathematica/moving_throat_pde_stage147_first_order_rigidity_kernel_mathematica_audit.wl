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

banner["STAGE 130 — FIRST-ORDER RIGIDITY KERNEL"];

Clear[p, x, eps, gBar, sBar];
$Assumptions = Element[{p, x, eps, gBar, sBar}, Reals] && p > 0 && 0 <= x <= 1;

kap = Pi/2;
c = Cos[Pi*x/2];
kq = Cosh[kap*(1 - x)]/Cosh[kap];

gFormula = 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
sFormula = p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));

rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
gMinus = FullSimplify[rF1 - Sqrt[1 + rF1^2]/2, Assumptions -> True];
pStar = p /. FindRoot[gFormula == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
gStar = N[gFormula /. p -> pStar, 40];
sStar = N[sFormula /. p -> pStar, 40];
gPrimeStar = N[D[gFormula, p] /. p -> pStar, 40];
sPrimeStar = N[D[sFormula, p] /. p -> pStar, 40];

sigmaStar = N[pStar/(1 - sStar/4), 40];
tStar = N[Sqrt[9*sigmaStar/20], 40];

aT = N[
  -(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)),
  30
];
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 30];

banner["Stage 147 audit: first-order rigidity coefficients"];
Print["Pi_* = ", fmt[N[pStar, 30]]];
Print["Sigma_* = ", fmt[sigmaStar]];
Print["T_* = ", fmt[tStar]];
Print["A_T = ", fmt[aT]];
Print["B_T = ", fmt[bT]];
Print["|A_T|/B_T = ", fmt[N[Abs[aT]/bT, 20]]];

rQMinus = FullSimplify[((gMinus - rF1)^2)/(1 + rF1^2), Assumptions -> True];
expectZero["R_q(g_minus)-1/4", rQMinus - 1/4];

dT = FullSimplify[eps*(aT*(gBar - gMinus) + bT*(sBar - sStar))];
Print["delta T_m = ", fmt[dT]];

wCenter = FullSimplify[aT*(c - gMinus) + bT*(kq - sStar)];
Print["Centered rigidity kernel W_*(x) = ", fmt[wCenter]];

Print[""];
Print["Stage 147 Mathematica audit passed."];

Exit[0];
