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

banner["STAGE 079 — FAMILY-1 QUADRUPOLE-DEMAND / PE MAP"];

Clear[pe, y];
$Assumptions = Element[{pe, y}, Reals] && pe > 0;

kappaF1 = 12321/5;
etaF1 = 37;
yF1 = y /. FindRoot[y*Tan[y] == etaF1, {y, 1.53}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yF1^2), 50];

Print["y_F1 = ", fmt[N[yF1, 30]]];
Print["A_F1 = ", fmt[aF1]];
Print["kappa_F1 = ", fmt[kappaF1]];
Print["eta_F1 = ", fmt[etaF1]];
expectApprox["Robin residual", N[yF1*Tan[yF1] - etaF1, 50], 0, 10^-28];

omega = FullSimplify[
  Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1)),
  Assumptions -> $Assumptions
];
omega0 = FullSimplify[Limit[omega, pe -> 0, Direction -> "FromAbove"]];
omegaInf = FullSimplify[Limit[omega, pe -> Infinity]];

Print["Omega(Pe) = ", fmt[omega]];
Print["Omega(0+) = ", fmt[omega0]];
Print["Omega(inf) = ", fmt[omegaInf]];
expectZero["Omega(0+) - 1", omega0 - 1];
expectZero["Omega(inf) - Pi/2", omegaInf - Pi/2];

zeta0 = N[aF1*omega0^2, 50];
zetaInf = N[aF1*omegaInf^2, 50];

Print["zeta_F1(Pe) = A_F1 Omega(Pe)^2"];
Print["zeta_F1(0+) = ", fmt[zeta0]];
Print["zeta_F1(inf) = ", fmt[zetaInf]];
zeta0Sym = FullSimplify[Limit[aF1*omega^2, pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
zetaInfSym = FullSimplify[Limit[aF1*omega^2, pe -> Infinity], Assumptions -> $Assumptions];
expectApprox["zeta_F1(0+) - A_F1", N[zeta0Sym - aF1, 50], 0, 10^-40];
expectApprox["zeta_F1(inf) - A_F1 Pi^2/4", N[zetaInfSym - aF1*Pi^2/4, 50], 0, 10^-40];

zetaSeries = FullSimplify[Normal[Series[aF1*omega^2, {pe, 0, 1}]], Assumptions -> $Assumptions];
expectedSeries = FullSimplify[aF1*(1 + ((4 - Pi)/Pi)*pe), Assumptions -> $Assumptions];
seriesDiff = N[Expand[zetaSeries - expectedSeries], 50];
expectApprox["small-Pe constant coefficient check", Coefficient[seriesDiff, pe, 0], 0, 10^-28];
expectApprox["small-Pe linear coefficient check", Coefficient[seriesDiff, pe, 1], 0, 10^-28];
omegaPrime0 = FullSimplify[Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
Print["Omega'(0+) = ", fmt[omegaPrime0]];
expectApprox["Omega'(0+) - (4-Pi)/(2 Pi)", N[omegaPrime0 - (4 - Pi)/(2*Pi), 50], 0, 10^-40];

Print[""];
Print["Stage 079 Mathematica audit passed."];

Exit[0];
