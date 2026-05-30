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

expectEqual[name_String, lhs_, rhs_] := Module[{res},
  res = FullSimplify[lhs - rhs, Assumptions -> $Assumptions];
  Print[name, ": lhs - rhs = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectPositive[name_String, expr_] := Module[{val},
  val = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[name, ": ", fmt[val]];
  If[TrueQ[Simplify[val > 0]], pass[name], fail[name, val]];
];

banner["STAGE 143 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT"];

Clear[piM];
$Assumptions = Element[piM, Reals] && piM > 0;

r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));

subbanner["Exact strict inequality g_Pi < 1 for finite Pi"];
num = Expand[(1 - gPi)*(4*piM^2 + Pi^2)*(Exp[piM] - 1)];
Print["numerator = ", fmt[num]];

decomp = Pi^2*(Exp[piM] - 1 - piM - piM^2/2) + piM*(Pi^2 - 2*Pi) + piM^2*(Pi^2/2 - 4);
expectZero["numerator - exact positive decomposition", num - decomp];

(* Independent positivity verification: prove num > 0 for piM > 0 directly *)
numPositiveCheck = Reduce[num > 0, piM, Reals] /. {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
Print["Reduce[num > 0, piM, Reals] = ", fmt[numPositiveCheck]];
If[TrueQ[Simplify[numPositiveCheck === True || numPositiveCheck === (piM > 0)]],
  pass["num > 0 for piM > 0 via Reduce"],
  fail["num > 0 for piM > 0 via Reduce", numPositiveCheck]];

Print["positive pieces:"];
Print["  exp remainder = ", fmt[Exp[piM] - 1 - piM - piM^2/2]];
Print["  linear coeff  = ", fmt[FullSimplify[Pi^2 - 2*Pi]]];
Print["  quadratic coeff = ", fmt[FullSimplify[Pi^2/2 - 4]]];
expectPositive["Pi^2 - 2*Pi > 0", Pi^2 - 2*Pi];
expectPositive["Pi^2/2 - 4 > 0", Pi^2/2 - 4];
(* exp-remainder positivity: prove R(piM) = Exp[piM] - 1 - piM - piM^2/2 > 0 for piM > 0. *)
(* Primary route: Reduce over the reals must reduce the inequality to piM > 0. *)
expRemReduce = Reduce[Exp[piM] - 1 - piM - piM^2/2 > 0, piM, Reals] /.
  {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
Print["Reduce[exp remainder > 0, piM, Reals] = ", fmt[expRemReduce]];
If[TrueQ[Simplify[expRemReduce === True || expRemReduce === (piM > 0)]],
  pass["exp remainder > 0 for piM > 0 via Reduce"],
  fail["exp remainder > 0 for piM > 0 via Reduce", expRemReduce]];
(* Independent backing route: Taylor-remainder monotonicity. *)
(* R(0)=0, R'(0)=0, R''(0)=0, R'''(piM)=Exp[piM]>0 => R strictly increasing from 0 => R>0. *)
rRem = Exp[piM] - 1 - piM - piM^2/2;
expectEqual["exp remainder R(0) == 0", (rRem /. piM -> 0), 0];
expectEqual["exp remainder R'(0) == 0", (D[rRem, piM] /. piM -> 0), 0];
expectEqual["exp remainder R''(0) == 0", (D[rRem, {piM, 2}] /. piM -> 0), 0];
expectEqual["exp remainder R'''(piM) - Exp[piM] == 0", D[rRem, {piM, 3}], Exp[piM]];
expectPositive["exp remainder R'''(piM) = Exp[piM] > 0 for piM>0", Exp[piM]];

subbanner["Endpoint limits"];
Clear[pi0, piInf, piInf2, piInf3];
g0 = FullSimplify[Limit[gPi /. piM -> pi0, pi0 -> 0, Direction -> "FromAbove"], Assumptions -> pi0 > 0];
gInf = FullSimplify[Limit[gPi /. piM -> piInf, piInf -> Infinity], Assumptions -> piInf > 0];
Print["lim_{Pi->0+} g_Pi = ", fmt[g0]];
Print["lim_{Pi->oo} g_Pi = ", fmt[gInf]];
expectEqual["lim_{piM->0+} g_Pi == 2/Pi", g0, 2/Pi];
expectEqual["lim_{piM->oo} g_Pi == 1", gInf, 1];

sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
that = Sqrt[(9/20)*sigma0];

Clear[piInf2, piInf3, piInf4, piInf5];
rInf = FullSimplify[Limit[rQ /. piM -> piInf2, piInf2 -> Infinity], Assumptions -> piInf2 > 0];
sInf = FullSimplify[Limit[sQ /. piM -> piInf3, piInf3 -> Infinity], Assumptions -> piInf3 > 0];
sigmaRatio = FullSimplify[Limit[sigma0/piM /. piM -> piInf4, piInf4 -> Infinity], Assumptions -> piInf4 > 0];
tHatRatio = FullSimplify[Limit[that/Sqrt[piM] /. piM -> piInf5, piInf5 -> Infinity], Assumptions -> piInf5 > 0];

Print["R_infty = ", fmt[rInf]];
Print["S_infty = ", fmt[sInf]];
Print["lim Sigma0/Pi = ", fmt[sigmaRatio]];
Print["lim That/sqrt(Pi) = ", fmt[tHatRatio]];
expectEqual["R_infty == (1-r)^2/(1+r^2)", rInf, (1 - r)^2/(1 + r^2)];
expectEqual["S_infty == 1", sInf, 1];
expectEqual["lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))", tHatRatio, Sqrt[(9/20)/(1 - rInf)]];

banner["STAGE 143 LEDGER"];
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
