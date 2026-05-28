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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 130 — EXACT MOUTH-BIAS MAP AND FAMILY-1 COMPENSATION POINT"];

Clear[z, lM, piM];
$Assumptions = Element[{z, lM, piM}, Reals] && lM > 0 && piM > 0;

gMinus = N[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi), 80];

sigma = piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]));
f = Cos[Pi*z/(2*lM)];

gPi = FullSimplify[Integrate[sigma*f, {z, 0, lM}], Assumptions -> $Assumptions];
gPiBoxed = 2*piM*(2*piM*Exp[piM] + Pi) / ((4*piM^2 + Pi^2) * (Exp[piM] - 1));
expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed];
eZ = FullSimplify[Integrate[sigma*z, {z, 0, lM}], Assumptions -> $Assumptions];
eFZ = FullSimplify[Integrate[sigma*f*z, {z, 0, lM}], Assumptions -> $Assumptions];
covId = FullSimplify[D[gPi, piM] + (eFZ - gPi*eZ)/lM, Assumptions -> $Assumptions];

Print["g_Pi = ", fmt[gPi]];
Print["Covariance identity residual = ", fmt[covId]];
expectZero["covariance identity", covId];

Clear[pi0, piInf];
g0 = FullSimplify[Limit[gPi /. piM -> pi0, pi0 -> 0, Direction -> "FromAbove"], Assumptions -> pi0 > 0];
gInf = FullSimplify[Limit[gPi /. piM -> piInf, piInf -> Infinity], Assumptions -> piInf > 0];
Print["limit Pi->0+ = ", fmt[g0]];
Print["limit Pi->oo = ", fmt[gInf]];
expectZero["uniform-source limit", g0 - 2/Pi];
expectZero["point-source limit", gInf - 1];

(* Strict monotonicity sweep: dg/dpiM > 0 for piM > 0 (notes boxed result) *)
dgPi = D[gPi, piM];
Module[{vals = {1/10, 1/2, 1, 15088/10000, 3, 10}, dv},
  Do[
    dv = N[dgPi /. piM -> v, 40];
    Print["dg/dpiM at piM=", fmt[v], " = ", fmt[dv]];
    If[TrueQ[dv > 0], pass["dg/dpiM > 0 at piM=" <> ToString[v]],
      fail["dg/dpiM > 0 at piM=" <> ToString[v], dv]],
    {v, vals}]
];

piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
Print["Pi_* = ", fmt[N[piStar, 30]]];
Print["x_* = 1/Pi_* = ", fmt[N[1/piStar, 30]]];
Print["g(Pi_*) = ", fmt[N[gPi /. piM -> piStar, 30]]];
Print["g'(Pi_*) = ", fmt[N[D[gPi, piM] /. piM -> piStar, 30]]];
expectApprox["Family-1 compensation point", gPi /. piM -> piStar, gMinus, 10^-20];

Print[""];
Print["Stage 130 Mathematica audit passed."];

Exit[0];
