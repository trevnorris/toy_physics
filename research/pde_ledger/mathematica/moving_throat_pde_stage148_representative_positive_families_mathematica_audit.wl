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

banner["STAGE 148 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES"];

Clear[p, lam];
$Assumptions = Element[{p, lam}, Reals] && p > 0 && 0 <= lam <= 1;

kap = Pi/2;
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
(* dSigma in terms of dPi (= -(g_target - gStar)/gPrimeStar) and dS = (s_target - sStar) *)
dSigmaOfDeltas[dPi_, dS_] := dPi/(1 - sStar/4) + pStar*dS/(4*(1 - sStar/4)^2);
dTOfDeltas[dG_, dS_] := Module[{dPi},
  dPi = -dG/gPrimeStar;
  N[(9/(40*tStar))*dSigmaOfDeltas[dPi, dS], 30]
];

gU = N[2/Pi, 30];
sU = N[2*Tanh[Pi/2]/Pi, 30];
dPiU = N[-(gU - gStar)/gPrimeStar, 30];
dTU = dTOfDeltas[gU - gStar, sU - sStar];

Print["uniform: g_u = ", fmt[gU]];
Print["uniform: S_u = ", fmt[sU]];
Print["uniform: dPi/eps = ", fmt[dPiU]];
Print["uniform: dT/eps = ", fmt[dTU]];

gD = N[Pi/4, 30];
sD = N[(1 + Sinh[Pi/2])/(2*Cosh[Pi/2]), 30];
dPiD = N[-(gD - gStar)/gPrimeStar, 30];
dTD = dTOfDeltas[gD - gStar, sD - sStar];

Print["derivative: g_d = ", fmt[gD]];
Print["derivative: S_d = ", fmt[sD]];
Print["derivative: dPi/eps = ", fmt[dPiD]];
Print["derivative: dT/eps = ", fmt[dTD]];

gLam = FullSimplify[(1 - lam)*(2/Pi) + lam*(Pi/4), Assumptions -> True];
sLam = FullSimplify[(1 - lam)*(2*Tanh[Pi/2]/Pi) + lam*((1 + Sinh[Pi/2])/(2*Cosh[Pi/2])), Assumptions -> True];
dPiLam = FullSimplify[-(gLam - gStar)/gPrimeStar];
dTLam = FullSimplify[dTOfDeltas[gLam - gStar, sLam - sStar]];

Print["dPi_lambda/eps = ", fmt[Expand[dPiLam]]];
Print["dT_lambda/eps = ", fmt[Expand[dTLam]]];

lamPiZero = FullSimplify[lam /. First[Solve[gLam == gMinus, lam, Reals]], Assumptions -> True];
lamTZero = N[lam /. First[Solve[dTLam == 0, lam, Reals]], 30];
Print["lambda_(Pi,0) = ", fmt[N[lamPiZero, 30]]];
Print["lambda_(T,0) = ", fmt[lamTZero]];
Print["1 - lambda_(Pi,0) = ", fmt[N[1 - lamPiZero, 30]]];

xiStarClosed = (-37*Sqrt[3] - 5*Pi^2 + 2*Sqrt[4107 - 100*Pi^2]) / (5*(8 - Pi^2));
Print["xi_* (Stage 126 closed form) = ", fmt[N[xiStarClosed, 30]]];
expectZero["(1-lambda_(Pi,0)) - xi_*", FullSimplify[(1 - lamPiZero) - xiStarClosed]];

Print[""];
Print["Stage 148 Mathematica audit passed."];

Exit[0];
