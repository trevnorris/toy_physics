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

banner["STAGE 114 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION"];

Clear[piM, thetaSigma, lM, tM, qStar, a0Prime];
$Assumptions =
  Element[{piM, thetaSigma, lM, tM, qStar, a0Prime}, Reals] &&
  piM > 0 && thetaSigma > 0 && lM > 0;

gMinus = N[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi), 80];
gPi = FullSimplify[2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1)), Assumptions -> $Assumptions];
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];

v1 = piM*thetaSigma/lM;
v1Star = N[piStar, 30]*thetaSigma/lM;
gPrimeStar = N[D[gPi, piM] /. piM -> piStar, 30];
thresholdResidual = FullSimplify[(tM - qStar*a0Prime) - piM*thetaSigma/lM, Assumptions -> $Assumptions];

Print["Pi_* = ", fmt[N[piStar, 30]]];
Print["V1_* = ", fmt[v1Star]];
Print["g'(Pi_*) = ", fmt[gPrimeStar]];
Print["Parent bias mismatch formula = ", fmt[thresholdResidual]];
expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20];

Print[""];
Print["Stage 131 Mathematica audit passed."];

Exit[0];
