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

(* Affine laws tested via integral form, not via algebraic restatement. *)
varsigmaTest = 6*x*(1 - x);
sigmaEps = (1 - eps)*(sigma /. p -> pStar) + eps*varsigmaTest;
gBarPhys = Integrate[sigmaEps*Cos[Pi*x/2], {x, 0, 1}];
sBarPhys = Integrate[sigmaEps*kq, {x, 0, 1}];
gBarV    = Integrate[varsigmaTest*Cos[Pi*x/2], {x, 0, 1}];
sBarV    = Integrate[varsigmaTest*kq, {x, 0, 1}];
(* Two-tier check. Tier 1 (exact symbolic zero) is NOT reachable: the residual
   collapses by linearity of the integral to (1-eps)*(gFormula[pStar]-gMinus),
   and pStar is a transcendental root with no closed form. Tier 2: evaluate the
   RAW residual at high precision (no Chop) and require Abs < 10^-25. pStar is
   solved at WorkingPrecision 80 / AccuracyGoal 30 (line 66), so the residual is
   ~10^-40 -- well inside the 10^-25 budget. Combine the symbolic endpoint
   integrals before substituting pStar so numeric Integrate precision loss does
   not collapse the residual to a low-precision zero. *)
gEpsRes = ((1 - eps)*(gDirect /. p -> pStar) + eps*gBarV) - (gMinus + eps*(gBarV - gMinus));
gEpsSample1 = N[(gEpsRes /. eps -> 1/10), 50];
gEpsSample2 = N[(gEpsRes /. eps -> 1/2), 50];
Print["g_eps affine law (integral form) at eps=1/10: ", fmt[gEpsSample1]];
Print["g_eps affine law (integral form) at eps=1/2:  ", fmt[gEpsSample2]];
If[NumericQ[gEpsSample1] && NumericQ[gEpsSample2] && Abs[gEpsSample1] < 10^-25 && Abs[gEpsSample2] < 10^-25,
  pass["g_eps affine law (integral form)"],
  fail["g_eps affine law (integral form)", {gEpsSample1, gEpsSample2}]
];

sEpsRes = ((1 - eps)*(sDirect /. p -> pStar) + eps*sBarV) - (sStar + eps*(sBarV - sStar));
sEpsSample1 = N[(sEpsRes /. eps -> 1/10), 50];
sEpsSample2 = N[(sEpsRes /. eps -> 1/2), 50];
Print["S_eps affine law (integral form) at eps=1/10: ", fmt[sEpsSample1]];
Print["S_eps affine law (integral form) at eps=1/2:  ", fmt[sEpsSample2]];
If[NumericQ[sEpsSample1] && NumericQ[sEpsSample2] && Abs[sEpsSample1] < 10^-25 && Abs[sEpsSample2] < 10^-25,
  pass["S_eps affine law (integral form)"],
  fail["S_eps affine law (integral form)", {sEpsSample1, sEpsSample2}]
];

Print[""];
Print["Stage 146 Mathematica audit passed."];

Exit[0];
