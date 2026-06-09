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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

sphereAvg[expr_] := FullSimplify[
  Integrate[Integrate[expr*Sin[th], {ph, 0, 2*Pi}], {th, 0, Pi}]/(4*Pi),
  Assumptions -> $Assumptions
];

banner["STAGE 169 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE"];

Clear[x20, x21, x22, y20, y21, y22, eps, x1, y1, x0, y0, th, ph, x0s, e, xiL, xiv, xiT, iGrp, g, r];
$Assumptions =
  Element[
    {x20, x21, x22, y20, y21, y22, eps, x1, y1, x0, y0, th, ph, x0s, e, xiL, xiv, xiT, iGrp, g, r},
    Reals
  ] && x0s > 0 && g > 0 && r > 0;

x = {x20, x21, x22};
y = {y20, y21, y22};
gMat = DiagonalMatrix[{1, 2, 2}];
eBar = {1, 1, 1};

xBar = (x20 + 2*x21 + 2*x22)/5;
yBar = (y20 + 2*y21 + 2*y22)/5;
aX = (2*x20 - x21 - x22)/10;
bX = (x21 - x22)/2;
aY = (2*y20 - y21 - y22)/10;
bY = (y21 - y22)/2;

dX = x - xBar*eBar;
dY = y - yBar*eBar;
iXY = FullSimplify[(dX . gMat . dY)/5, Assumptions -> $Assumptions];
iXYExpected = FullSimplify[4*aX*aY + (4/5)*bX*bY, Assumptions -> $Assumptions];

banner["Grouped bilinear invariant"];
Print["I[x,y] = ", fmt[iXY]];
expectZero["I[x,y] - [4 a_x a_y + 4/5 b_x b_y]", iXY - iXYExpected];

subsAxis = {
  x20 -> x0 + eps*x1,
  x21 -> x0 + eps*x1/2,
  x22 -> x0 - eps*x1,
  y20 -> y0 + eps*y1,
  y21 -> y0 + eps*y1/2,
  y22 -> y0 - eps*y1
};

aXAxis = FullSimplify[aX /. subsAxis, Assumptions -> $Assumptions];
bXAxis = FullSimplify[bX /. subsAxis, Assumptions -> $Assumptions];
iXYAxis = FullSimplify[iXY /. subsAxis, Assumptions -> $Assumptions];

banner["Weak axisymmetric grouped law"];
Print["a_x(axis) = ", fmt[aXAxis]];
Print["b_x(axis) = ", fmt[bXAxis]];
expectZero["b_x - 3 a_x", bXAxis - 3*aXAxis];
expectZero["I_axis - (7/10) eps^2 x1 y1", iXYAxis - (7*eps^2*x1*y1)/10];

banner["Explicit harmonic average"];
y20Harm = Sqrt[5/(16*Pi)]*(3*Cos[th]^2 - 1);
y20Avg = sphereAvg[y20Harm];
y20SqAvg = sphereAvg[y20Harm^2];
Print["<Y20> = ", fmt[y20Avg]];
Print["<Y20^2> = ", fmt[y20SqAvg]];
expectZero["average Y20", y20Avg];
expectZero["average Y20^2 - 1/(4 pi)", y20SqAvg - 1/(4*Pi)];

logSeries = Expand[Normal[Series[Log[x0s*(1 + e*y20Harm)], {e, 0, 2}]]];
logAvg = FullSimplify[sphereAvg[logSeries], Assumptions -> $Assumptions];
linCoeff = FullSimplify[D[logAvg, e] /. e -> 0, Assumptions -> $Assumptions];
Print["<log(X0 (1 + e Y20))> through O(e^2) = ", fmt[logAvg]];
expectZero["linear coefficient in averaged log-observable", linCoeff];

banner["Stage 168 weighted transport"];
s = Sqrt[1 + r^2];
epsPerp = g*xiT*iGrp + (g + 1/(2*s))*xiv*iGrp + (2*g + 3/(4*s))*xiL*iGrp;
xiPerp = g*xiT + (g + 1/(2*s))*xiv + (2*g + 3/(4*s))*xiL;
expectZero["eps_perp - Xi_perp Igrp", epsPerp - xiPerp*iGrp];

rExact = Sqrt[4107 - 100*Pi^2]/(10*Pi);   (* canonical Family-1 radius *)
rNum = N[rExact, 30];
gNum = SetPrecision[0.758035078944663, 30];
Print["Numeric Xi_perp combination = ", fmt[N[xiPerp /. {g -> gNum, r -> rNum}, 20]]];

coeffT = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiT], 20];
coeffV = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiv], 20];
coeffL = N[Coefficient[xiPerp /. {g -> gNum, r -> rNum}, xiL], 20];
Module[{checks},
  checks = {
    {"Xi_perp coeff on xiT", coeffT, 0.758035078944663},
    {"Xi_perp coeff on xiv", coeffV, 1.00314310113848},
    {"Xi_perp coeff on xiL", coeffL, 1.88373219118005}
  };
  Do[
    With[{nm = c[[1]], got = c[[2]], want = c[[3]]},
      Print[nm, " = ", fmt[got], " (paper ", fmt[want], ")"];
      If[Abs[got - want] > 10^-12, fail[nm, got - want], pass[nm]]
    ], {c, checks}]
];

Print[""];
Print["Carry-forward formulas:"];
Print["  I[x,y] = (1/5) (x-xbar ebar)^T Ggrp (y-ybar ebar) = 4 a_x a_y + (4/5) b_x b_y"];
Print["  Pure grouped P2 anisotropy has no linear scalar feed-down on an isotropic branch."];
Print["  Weak axisymmetric branch: b = 3 a and I[x,y] = (7/10) eps^2 x^(1) y^(1)."];
Print["  Therefore grouped-lane anisotropy enters eps_L, eps_v, eps_T, eps_perp only quadratically."];

Print[""];
Print["Stage 169 Mathematica audit passed."];

Exit[0];
