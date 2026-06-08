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

banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];

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

sDirect = FullSimplify[Integrate[sigma*kq, {x, 0, 1}], Assumptions -> p > 0];
expectZero["S_q(Pi) direct-formula", sDirect - sFormula];

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

(* Convex-family affine moment laws tested by direct quadrature of the
   assembled profile against closed-form intercepts g_* and S_*. *)
varsigmaTest = 6*x*(1 - x);
sigmaEps = (1 - eps)*(sigma /. p -> pStar) + eps*varsigmaTest;
momentQuad[expr_] := NIntegrate[
  Evaluate[expr], {x, 0, 1},
  WorkingPrecision -> 80, AccuracyGoal -> 35, PrecisionGoal -> 35
];
gStarClosed = N[gFormula /. p -> pStar, 50];
sStarClosed = N[sFormula /. p -> pStar, 50];
gBarV = momentQuad[varsigmaTest*Cos[Pi*x/2]];
sBarV = momentQuad[varsigmaTest*kq];
gSlope = N[gBarV - gStarClosed, 50];
sSlope = N[sBarV - sStarClosed, 50];
Print["g_eps convex affine moment law nonzero slope |gbar_v - g_*| = ", fmt[Abs[gSlope]]];
Print["S_eps convex affine moment law nonzero slope |Sbar_v - S_*| = ", fmt[Abs[sSlope]]];
If[!TrueQ[Abs[gSlope] > 10^-3],
  fail["g_eps convex affine moment law nonzero slope guard", gSlope]
];
If[!TrueQ[Abs[sSlope] > 10^-3],
  fail["S_eps convex affine moment law nonzero slope guard", sSlope]
];

Do[
  Module[{epsVal = sample[[1]], label = sample[[2]], sigmaEpsSample, gBarEps, sBarEps, gRes, sRes},
    sigmaEpsSample = sigmaEps /. eps -> epsVal;
    gBarEps = momentQuad[sigmaEpsSample*Cos[Pi*x/2]];
    sBarEps = momentQuad[sigmaEpsSample*kq];
    gRes = N[gBarEps - (gStarClosed + epsVal*(gBarV - gStarClosed)), 50];
    sRes = N[sBarEps - (sStarClosed + epsVal*(sBarV - sStarClosed)), 50];
    Print["g_eps convex affine moment law via direct quadrature at ", label, ": ", fmt[gRes]];
    Print["S_eps convex affine moment law via direct quadrature at ", label, ": ", fmt[sRes]];
    If[!TrueQ[NumericQ[gRes] && Abs[gRes] < 10^-25],
      fail["g_eps convex affine moment law via direct quadrature " <> label, gRes]
    ];
    If[!TrueQ[NumericQ[sRes] && Abs[sRes] < 10^-25],
      fail["S_eps convex affine moment law via direct quadrature " <> label, sRes]
    ];
  ],
  {sample, {{1/10, "eps=1/10"}, {1/2, "eps=1/2"}}}
];
pass["g_eps convex affine moment law via direct quadrature with closed-form g_* intercept and nonzero slope"];
pass["S_eps convex affine moment law via direct quadrature with closed-form S_* intercept and nonzero slope"];

Print[""];
Print["Stage 146 Mathematica audit passed."];

Exit[0];
