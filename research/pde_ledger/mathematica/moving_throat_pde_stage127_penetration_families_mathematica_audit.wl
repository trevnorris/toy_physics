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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES"];

Clear[x];
$Assumptions = Element[x, Reals] && x > 0;

rDisc = Sqrt[4107 - 100*Pi^2];
gMinus = N[(2*rDisc - 37*Sqrt[3])/(20*Pi), 80];

gSlab = FullSimplify[2*Sin[Pi*x/2]/(Pi*x), Assumptions -> $Assumptions];
gExp = FullSimplify[2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x])), Assumptions -> $Assumptions];

Print["g_slab(x) = ", fmt[gSlab]];
Print["g_exp(x) = ", fmt[gExp]];

xStarSlab = x /. FindRoot[gSlab == gMinus, {x, 0.8}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 200];
xStarExp = x /. FindRoot[gExp == gMinus, {x, 0.66}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 200];

Print["g_-^F1 = ", fmt[N[gMinus, 30]]];
Print["x_*^slab = ", fmt[N[xStarSlab, 30]]];
Print["x_*^exp  = ", fmt[N[xStarExp, 30]]];

expectApprox["slab compensation root", gSlab /. x -> xStarSlab, gMinus, 10^-20];
expectApprox["exponential compensation root", gExp /. x -> xStarExp, gMinus, 10^-20];

Print[""];
Print["Stage 127 Mathematica audit passed."];

Exit[0];
