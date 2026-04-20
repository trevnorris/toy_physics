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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];

Clear[p, x, gBar, sBar, eps];
$Assumptions = Element[{p, x, gBar, sBar, eps}, Reals] && p > 0 && 0 <= x <= 1;

kap = Pi/2;
sigma = p*Exp[-p*x]/(1 - Exp[-p]);
kq = Cosh[kap*(1 - x)]/Cosh[kap];

gFormula = 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
sFormula = p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));

gDirect = FullSimplify[Integrate[sigma*Cos[Pi*x/2], {x, 0, 1}], Assumptions -> p > 0];

Print["g(Pi) = ", fmt[gFormula]];
Print["S_q(Pi) = ", fmt[sFormula]];
expectZero["g(Pi) direct-formula", gDirect - gFormula];

nint[expr_] := NIntegrate[Evaluate[expr], {x, 0, 1}, WorkingPrecision -> 50, AccuracyGoal -> 20, PrecisionGoal -> 20];
Do[
  Module[{val = vv, numInt, numFormula},
    numInt = nint[(sigma*kq) /. p -> val];
    numFormula = N[sFormula /. p -> val, 40];
    Print["kernel check at Pi=", fmt[val], ": integral=", fmt[numInt], ", formula=", fmt[numFormula]];
    expectApprox["kernel formula sample " <> ToString[InputForm[val]], numInt, numFormula, 10^-12];
  ],
  {vv, {1, 3/2, 5/2}}
];

rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1];
gMinus = FullSimplify[rF1 - Sqrt[1 + rF1^2]/2, Assumptions -> True];
pStar = p /. FindRoot[gFormula == gMinus, {p, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];

Print["Pi_* = ", fmt[N[pStar, 30]]];

gStar = N[gFormula /. p -> pStar, 40];
sStar = N[sFormula /. p -> pStar, 40];
gPrimeStar = N[D[gFormula, p] /. p -> pStar, 40];
sPrimeStar = N[D[sFormula, p] /. p -> pStar, 40];

Print["g_* = ", fmt[gStar]];
Print["S_* = ", fmt[sStar]];
Print["g_*' = ", fmt[gPrimeStar]];
Print["S_*' = ", fmt[sPrimeStar]];
expectApprox["Pi_* compensation point", gStar, N[gMinus, 40], 10^-20];

dPi = FullSimplify[-eps*(gBar - gMinus)/gPrimeStar];
dS = FullSimplify[eps*(sBar - sStar) + sPrimeStar*dPi];

banner["FIRST-ORDER CANONICAL RETUNING LAW"];
Print["delta Pi = ", fmt[dPi]];
Print["delta S = ", fmt[dS]];

gEps = (1 - eps)*gMinus + eps*gBar;
sEps = (1 - eps)*sStar + eps*sBar;
expectZero["g_eps affine law", Expand[gEps - (gMinus + eps*(gBar - gMinus))]];
resSEps = Chop[Expand[sEps - (sStar + eps*(sBar - sStar))]];
Print["S_eps affine law = ", fmt[resSEps]];
If[TrueQ[resSEps === 0], pass["S_eps affine law"], fail["S_eps affine law", resSEps]];

Print[""];
Print["Stage 129 Mathematica audit passed."];

Exit[0];
