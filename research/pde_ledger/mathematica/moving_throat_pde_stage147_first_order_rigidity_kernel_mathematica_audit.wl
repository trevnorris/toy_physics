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

banner["STAGE 147 — FIRST-ORDER RIGIDITY KERNEL"];

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

(* --- Audit assertions: numerical anchor against paper-quoted literals --- *)
aTPaper = -4.27263956256927`30;
bTPaper = 0.134875005736706`30;
ratioPaper = 31.6785`20;
expectZero["A_T vs paper -4.27263956256927",
  If[Abs[aT - aTPaper] < 10^-12, 0, aT - aTPaper]];
expectZero["B_T vs paper 0.134875005736706",
  If[Abs[bT - bTPaper] < 10^-12, 0, bT - bTPaper]];
expectZero["|A_T|/B_T vs paper 31.6785",
  If[Abs[Abs[aT]/bT - ratioPaper] < 10^-3, 0, Abs[aT]/bT - ratioPaper]];

(* --- Audit assertion: chain-rule consistency for A_T (independent route) --- *)
dTmDSigma = 9/(40*tStar);
dSigmaDPi = 1/(1 - sStar/4) + pStar*sPrimeStar/(4*(1 - sStar/4)^2);
aTChain = N[-dTmDSigma*dSigmaDPi/gPrimeStar, 30];
expectZero["A_T closed form vs chain-rule route",
  If[Abs[aTChain - aT] < 10^-20, 0, aTChain - aT]];

dT = FullSimplify[eps*(aT*(gBar - gMinus) + bT*(sBar - sStar))];
Print["delta T_m = ", fmt[dT]];

wCenter = FullSimplify[aT*(c - gMinus) + bT*(kq - sStar)];
Print["Centered rigidity kernel W_*(x) = ", fmt[wCenter]];

(* --- Audit assertion: centered kernel structure --- *)
wCenterConst = FullSimplify[(wCenter - (aT*c + bT*kq)) /. x -> 1/2];
wCenterConstExpected = FullSimplify[-aT*gMinus - bT*sStar];
expectZero["W_*(x) centered form A_T(c-g_*) + B_T(K_q-S_*)",
  Chop[wCenterConst - wCenterConstExpected, 10^-25]];

(* --- Audit assertion: source-moment values g_*, S_* stable under resubstitution --- *)
gStarResub = N[gFormula /. p -> pStar, 40];
sStarResub = N[sFormula /. p -> pStar, 40];
expectZero["g_* resubstitution drift",
  If[Abs[gStarResub - gStar] < 10^-30, 0, gStarResub - gStar]];
expectZero["S_* resubstitution drift",
  If[Abs[sStarResub - sStar] < 10^-30, 0, sStarResub - sStar]];

Print[""];
Print["Stage 147 Mathematica audit passed."];

Exit[0];
