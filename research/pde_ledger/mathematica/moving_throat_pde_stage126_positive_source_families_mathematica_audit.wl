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

banner["STAGE 109 — EXPLICIT POSITIVE SOURCE FAMILIES"];

Clear[z, lM, xi];
$Assumptions = Element[{z, lM, xi}, Reals] && lM > 0;

k = Pi/(2*lM);
sigmaMatch = k*Cos[k*z];
normMatch = FullSimplify[Integrate[sigmaMatch, {z, 0, lM}], Assumptions -> $Assumptions];
gMatch = FullSimplify[Integrate[sigmaMatch*Cos[k*z], {z, 0, lM}], Assumptions -> $Assumptions];

Print["sigma_match normalization = ", fmt[normMatch]];
Print["g_match = ", fmt[gMatch]];
expectZero["self-matched normalization", normMatch - 1];
expectZero["self-matched bias", gMatch - Pi/4];

rDisc = Sqrt[4107 - 100*Pi^2];
gMinus = FullSimplify[(2*rDisc - 37*Sqrt[3])/(20*Pi), Assumptions -> $Assumptions];
deltaG = FullSimplify[Pi/4 - gMinus, Assumptions -> $Assumptions];
tractionRatio = FullSimplify[(Pi/4)/gMinus, Assumptions -> $Assumptions];

Print["g_-^F1 = ", fmt[gMinus]];
Print["delta_g_match = ", fmt[deltaG]];
Print["T_- / T_match = ", fmt[tractionRatio]];
Print["g_match numeric = ", fmt[N[Pi/4, 20]]];
Print["g_- numeric     = ", fmt[N[gMinus, 20]]];
Print["traction ratio  = ", fmt[N[tractionRatio, 20]]];

sigmaXi = (1 - xi)*k*Cos[k*z] + xi/lM;
normXi = FullSimplify[Integrate[sigmaXi, {z, 0, lM}], Assumptions -> $Assumptions];
gXi = FullSimplify[Integrate[sigmaXi*Cos[k*z], {z, 0, lM}], Assumptions -> $Assumptions];

Print["sigma_xi normalization = ", fmt[normXi]];
Print["g_xi = ", fmt[gXi]];
expectZero["convex-family normalization", normXi - 1];

xiStar = FullSimplify[xi /. First[Solve[gXi == gMinus, xi, Reals]], Assumptions -> $Assumptions];
Print["xi_* = ", fmt[xiStar]];
Print["xi_* numeric = ", fmt[N[xiStar, 20]]];
expectZero["g_xi(xi_*) - g_-", (gXi /. xi -> xiStar) - gMinus];

Print["2/pi numeric = ", fmt[N[2/Pi, 20]]];
Print["pi/4 numeric = ", fmt[N[Pi/4, 20]]];
intervalCheck = N[2/Pi, 40] < N[gMinus, 40] < N[Pi/4, 40];
Print["Check 2/pi < g_- < pi/4 -> ", intervalCheck];
If[!TrueQ[intervalCheck], fail["interval check", {N[2/Pi, 20], N[gMinus, 20], N[Pi/4, 20]}]];

Print[""];
Print["Stage 126 Mathematica audit passed."];

Exit[0];
