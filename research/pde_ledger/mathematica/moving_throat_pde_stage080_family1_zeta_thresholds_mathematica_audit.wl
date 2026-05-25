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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 063 — FAMILY-1 ZETA THRESHOLDS"];

Clear[lambdaMu, pe];
$Assumptions = Element[{lambdaMu, pe}, Reals] && lambdaMu > 0 && pe > 0;

kappaF1 = 12321/5;
etaF1 = 37;
yRoot = y /. FindRoot[y*Tan[y] == etaF1, {y, 1.53}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yRoot^2), 50];
omegaFn[p_] := Pi*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
zetaFn[p_] := aF1*omegaFn[p]^2;
zetaMax = N[aF1*Pi^2/4, 50];

peSuffChi = ToExpression["96.5285247264386`25"]*lambdaMu^2;
peFailChi = ToExpression["11220.5441626259`25"]*lambdaMu^2;
peSuffJ = ToExpression["22.0062226330754`25"]*lambdaMu^2;
peFailJ = ToExpression["2558.01892349205`25"]*lambdaMu^2;

zetaSuffChi = zetaFn[peSuffChi];
zetaFailChi = zetaFn[peFailChi];
zetaSuffJ = zetaFn[peSuffJ];
zetaFailJ = zetaFn[peFailJ];

omegaIndep[p_] := Pi*p*(2*p*Exp[p] + Pi) / ((4*p^2 + Pi^2) * (Exp[p] - 1));
aF1Indep = N[(12321/5 + Pi^2/4) / (12321/5 + yRoot^2), 50];
zetaTargetSuffChi = N[aF1Indep * omegaIndep[ToExpression["96.5285247264386`30"]]^2, 40];
zetaTargetFailChi = N[aF1Indep * omegaIndep[ToExpression["11220.5441626259`30"]]^2, 40];
zetaTargetSuffJ   = N[aF1Indep * omegaIndep[ToExpression["22.0062226330754`30"]]^2, 40];
zetaTargetFailJ   = N[aF1Indep * omegaIndep[ToExpression["2558.01892349205`30"]]^2, 40];

Print["zeta_max^(F1) = ", fmt[zetaMax]];
Print["zeta_suff^(chi)(lambda_mu) = zeta_F1(96.5285247264386 lambda_mu^2)"];
Print["zeta_fail^(chi)(lambda_mu) = zeta_F1(11220.5441626259 lambda_mu^2)"];
Print["zeta_suff^(J)(lambda_mu) = zeta_F1(22.0062226330754 lambda_mu^2)"];
Print["zeta_fail^(J)(lambda_mu) = zeta_F1(2558.01892349205 lambda_mu^2)"];

zetaSuffChi1 = N[zetaSuffChi /. lambdaMu -> 1, 40];
zetaFailChi1 = N[zetaFailChi /. lambdaMu -> 1, 40];
zetaSuffJ1 = N[zetaSuffJ /. lambdaMu -> 1, 40];
zetaFailJ1 = N[zetaFailJ /. lambdaMu -> 1, 40];

Print["zeta_suff^(chi)(1) = ", fmt[zetaSuffChi1]];
Print["zeta_fail^(chi)(1) = ", fmt[zetaFailChi1]];
Print["zeta_suff^(J)(1) = ", fmt[zetaSuffJ1]];
Print["zeta_fail^(J)(1) = ", fmt[zetaFailJ1]];

expectApprox["zeta_suff^(chi)(1) numeric check", zetaSuffChi1, zetaTargetSuffChi, 10^-14];
expectApprox["zeta_fail^(chi)(1) numeric check", zetaFailChi1, zetaTargetFailChi, 10^-14];
expectApprox["zeta_suff^(J)(1) numeric check", zetaSuffJ1, zetaTargetSuffJ, 10^-14];
expectApprox["zeta_fail^(J)(1) numeric check", zetaFailJ1, zetaTargetFailJ, 10^-14];

limSuffChi = N[Limit[zetaSuffChi, lambdaMu -> Infinity], 40];
limFailChi = N[Limit[zetaFailChi, lambdaMu -> Infinity], 40];
limSuffJ = N[Limit[zetaSuffJ, lambdaMu -> Infinity], 40];
limFailJ = N[Limit[zetaFailJ, lambdaMu -> Infinity], 40];

Print["limit zeta_suff^(chi) = ", fmt[limSuffChi]];
Print["limit zeta_fail^(chi) = ", fmt[limFailChi]];
Print["limit zeta_suff^(J) = ", fmt[limSuffJ]];
Print["limit zeta_fail^(J) = ", fmt[limFailJ]];

expectApprox["limit zeta_suff^(chi) -> zeta_max", limSuffChi, zetaMax, 10^-14];
expectApprox["limit zeta_fail^(chi) -> zeta_max", limFailChi, zetaMax, 10^-14];
expectApprox["limit zeta_suff^(J) -> zeta_max", limSuffJ, zetaMax, 10^-14];
expectApprox["limit zeta_fail^(J) -> zeta_max", limFailJ, zetaMax, 10^-14];
expectTrue["zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max", zetaSuffChi1 < zetaFailChi1 < zetaMax];
expectTrue["zeta_suff^(J)(1) < zeta_fail^(J)(1) < zeta_max", zetaSuffJ1 < zetaFailJ1 < zetaMax];
expectTrue["zeta_suff^(J)(1) <= zeta_suff^(chi)(1)", zetaSuffJ1 <= zetaSuffChi1];
expectTrue["zeta_fail^(J)(1) <= zeta_fail^(chi)(1)", zetaFailJ1 <= zetaFailChi1];

Print[""];
Print["Stage 080 Mathematica audit passed."];

Exit[0];
