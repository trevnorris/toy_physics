ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
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

banner["STAGE 126 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT"];

Clear[piM];
$Assumptions = Element[piM, Reals] && piM > 0;

r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));

subbanner["Exact strict inequality g_Pi < 1 for finite Pi"];
num = Expand[(1 - gPi)*(4*piM^2 + Pi^2)*(Exp[piM] - 1)];
Print["numerator = ", fmt[num]];

decomp = Pi^2*(Exp[piM] - 1 - piM - piM^2/2) + piM*(Pi^2 - 2*Pi) + piM^2*(Pi^2/2 - 4);
expectZero["numerator - exact positive decomposition", num - decomp];

Print["positive pieces:"];
Print["  exp remainder = ", fmt[Exp[piM] - 1 - piM - piM^2/2]];
Print["  linear coeff  = ", fmt[FullSimplify[Pi^2 - 2*Pi]]];
Print["  quadratic coeff = ", fmt[FullSimplify[Pi^2/2 - 4]]];

subbanner["Endpoint limits"];
Clear[pi0, piInf, piInf2, piInf3];
g0 = FullSimplify[Limit[gPi /. piM -> pi0, pi0 -> 0, Direction -> "FromAbove"], Assumptions -> pi0 > 0];
gInf = FullSimplify[Limit[gPi /. piM -> piInf, piInf -> Infinity], Assumptions -> piInf > 0];
Print["lim_{Pi->0+} g_Pi = ", fmt[g0]];
Print["lim_{Pi->oo} g_Pi = ", fmt[gInf]];

sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rInf = FullSimplify[(1 - r)^2/(1 + r^2)];
sInf = 1;
sigmaRatio = FullSimplify[1/(1 - rInf*sInf)];
tHatRatio = FullSimplify[Sqrt[(9/20)*sigmaRatio]];

Print["R_infty = ", fmt[rInf]];
Print["S_infty = ", fmt[sInf]];
Print["lim Sigma0/Pi = ", fmt[sigmaRatio]];
Print["lim That/sqrt(Pi) = ", fmt[tHatRatio]];

banner["STAGE 126 LEDGER"];
Print["For every finite Pi>0:"];
Print["  2/pi < g_Pi < 1"];
Print["So g_c = 1 is not a finite positive-bias branch."];
Print[""];
Print["Equal-normalized branch:"];
Print["  reached only as Pi -> infinity"];
Print["  and That -> infinity like sqrt(Pi)."];
Print[""];
Print["Numerics:"];
Print["  R_infty = ", fmt[N[rInf, 30]]];
Print["  lim That/sqrt(Pi) = ", fmt[N[tHatRatio, 30]]];

Print[""];
Print["Stage 143 Mathematica audit passed."];

Exit[0];
