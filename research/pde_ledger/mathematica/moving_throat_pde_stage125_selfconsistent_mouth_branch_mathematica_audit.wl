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

expectApprox[name_String, value_, target_, tol_] := Module[{delta},
  delta = N[value - target, 40];
  Print[name, " residual = ", fmt[delta]];
  If[TrueQ[Abs[delta] < tol], pass[name], fail[name, delta]];
];

banner["STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW"];

Clear[piM];
$Assumptions = Element[piM, Reals] && piM > 0;

r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
tHat = Sqrt[(9/20)*sigma0];

subbanner["Core-to-mouth reduction"];
Print["r_F1 = ", fmt[r]];
Print["g_Pi = ", fmt[gPi]];
Print["S_q(Pi) = ", fmt[sQ]];
Print["R_q(Pi) = ", fmt[rQ]];
Print["Sigma0(Pi) = ", fmt[sigma0]];
Print["That(Pi) = ", fmt[tHat]];

gMinus = FullSimplify[r - Sqrt[1 + r^2]/2, Assumptions -> $Assumptions];
rQMinus = FullSimplify[(gMinus - r)^2/(1 + r^2), Assumptions -> $Assumptions];
expectZero["R_q(g_minus)-1/4", rQMinus - 1/4];

subbanner["Canonical compensation point"];
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30, MaxIterations -> 100];
gStar = N[gPi /. piM -> piStar, 30];
rQStar = N[rQ /. piM -> piStar, 30];
sQStar = N[sQ /. piM -> piStar, 30];
sigmaStar = N[sigma0 /. piM -> piStar, 30];
tHatStar = N[tHat /. piM -> piStar, 30];

Print["g_minus^F1 = ", fmt[N[gMinus, 30]]];
Print["Pi_*       = ", fmt[N[piStar, 30]]];
Print["g(Pi_*)    = ", fmt[gStar]];
Print["R_q(Pi_*)  = ", fmt[rQStar]];
Print["S_q(Pi_*)  = ", fmt[sQStar]];
Print["Sigma0(Pi_*) = ", fmt[sigmaStar]];
Print["That(Pi_*)   = ", fmt[tHatStar]];
expectApprox["Pi_* compensation solve", gStar, N[gMinus, 30], 10^-12];

banner["STAGE 125 LEDGER"];
Print["Self-consistent Family-1 mouth branch:"];
Print["  Pi = Sigma0 * [1 - R_q(Pi) S_q(Pi)]"];
Print["  Sigma0(Pi) = Pi / (1 - R_q(Pi) S_q(Pi))"];
Print["  That(Pi)   = sqrt(9 Sigma0(Pi) / 20)"];
Print[""];
Print["Canonical point:"];
Print["  Pi_*       = ", fmt[N[piStar, 30]]];
Print["  Sigma0_*   = ", fmt[sigmaStar]];
Print["  That_*     = ", fmt[tHatStar]];

Print[""];
Print["Stage 125 Mathematica audit passed."];

Exit[0];
